#ifdef CUDA

#include "PEList.h"
#include "HostCUDA.h"

/// @brief Each TreePiece on a given PE checks in as its tree walk completes
///        Once all TreePieces are done, launch a gravity kernel on the GPU
/// @param treePiece A reference to the TreePiece that checked in
void PEList::finishWalk(TreePiece *treePiece) {
    vtpLocal.push_back(treePiece);

    // On first call, find the total number of active pieces on this PE.
    // The charm++ location manager gives us this count in cTreePieces
    if(vtpLocal.size() == 1) {
        CkLocMgr *locMgr = treeProxy.ckLocMgr();
        locMgr->iterate(cTreePieces);
    }

    // check if we have everyone
    if(vtpLocal.size() < cTreePieces.count)
        return;

    // bucketMarkers[i+1] is needed to determine # of IL entries per bucket
    if(finalBucketMarker != -1)
	bucketMarkers.push_back(finalBucketMarker);

    finishCb = new CkCallback(CkIndex_PEList::finishWalkCb(), CkMyPe(), thisProxy);

    // If the DataManager device pointer is NULL, the GPU data transfer is
    // still in progress and we need to delay the kernel launch
    bool dataReady = (!bRemote &&
        dMProxy.ckLocalBranch()->bLocalDataTransferred.load()) ||
        (bRemote && dMProxy.ckLocalBranch()->bRemoteDataTransferred.load());
    if (!dataReady)
        bKernelDelayed = 1;
    else
        launchKernel();
}

/// @brief Called from DataManager after remote transfer finishes. Launch our kernel if it was delayed
void PEList::tryLaunchDelayedKernel() {
    if (bKernelDelayed) {
        launchKernel();
    }
}

void PEList::finishWalkCb() {
     dMProxy.ckLocalBranch()->transferParticleVarsBack();
     reset();
}

/// @brief Launch the corresponding CUDA kernel, depending what type of request this was
void PEList::launchKernel() {
    request = new CudaRequest;

    // Ensure required data is present on the GPU
    CkAssert(dMProxy.ckLocalBranch()->d_localMoments != nullptr);
    CkAssert(dMProxy.ckLocalBranch()->d_localParts != nullptr);
    CkAssert(dMProxy.ckLocalBranch()->d_localVars != nullptr);
    // The following checks can fail if remote prefetch does not bring in any data.
    // if (bRemote) {
    //    CkAssert(dMProxy.ckLocalBranch()->d_remoteParts != nullptr);
    //    CkAssert(dMProxy.ckLocalBranch()->d_remoteMoments != nullptr);
    // }

    request->d_localMoments = dMProxy.ckLocalBranch()->d_localMoments;
    request->d_localParts = dMProxy.ckLocalBranch()->d_localParts;
    request->d_localVars = dMProxy.ckLocalBranch()->d_localVars;
    request->d_remoteParts = dMProxy.ckLocalBranch()->d_remoteParts;
    request->d_remoteMoments = dMProxy.ckLocalBranch()->d_remoteMoments;
    request->stream = stream;

    request->numBucketsPlusOne = bucketSizes.size()+1;

    request->node = bNode;
    request->remote = bRemote;

    request->fperiod = fperiod;

    const char* funcTag = "PEList::finish";
    if (iList.size() > 0) {
      hostMalloc((void**)&iListHost, iList.size()*sizeof(ILCell), funcTag);
      memcpy(iListHost, iList.data(),  iList.size()*sizeof(ILCell));
      hostMalloc((void**)&bucketMarkersHost, bucketMarkers.size()*sizeof(int), funcTag);
      memcpy(bucketMarkersHost, bucketMarkers.data(), bucketMarkers.size()*sizeof(int));
      hostMalloc((void**)&bucketStartsHost, bucketStarts.size()*sizeof(int), funcTag);
      memcpy(bucketStartsHost, bucketStarts.data(), bucketStarts.size()*sizeof(int));
      hostMalloc((void**)&bucketSizesHost, bucketSizes.size()*sizeof(int), funcTag);
      memcpy(bucketSizesHost, bucketSizes.data(), bucketSizes.size()*sizeof(int));
    }

    if (missedParts.size() > 0) {
      hostMalloc((void**)&missedPartsHost, missedParts.size()*sizeof(CompactPartData), funcTag);
      memcpy(missedPartsHost, missedParts.data(), missedParts.size()*sizeof(CompactPartData));
    }
    if (missedNodes.size() > 0) {
      hostMalloc((void**)&missedNodesHost, missedNodes.size()*sizeof(CudaMultipoleMoments), funcTag);
      memcpy(missedNodesHost, missedNodes.data(), missedNodes.size()*sizeof(CudaMultipoleMoments));
    }

    request->list = iListHost;
    request->missedParts = missedPartsHost;
    request->missedNodes = missedNodesHost;
    request->sMissed = bNode ? missedNodes.size()*sizeof(CudaMultipoleMoments) : missedParts.size()*sizeof(CompactPartData);
    request->bucketMarkers = bucketMarkersHost;
    request->bucketStarts = bucketStartsHost;
    request->bucketSizes = bucketSizesHost;
    request->numInteractions = iList.size();

    request->cb = finishCb;

    void (*transferFunc)(CudaRequest*);
    if (bNode) {
        transferFunc = bRemote ? PEListNodeListDataTransferRemote : PEListNodeListDataTransferLocal;
        if (bResume) {
                transferFunc = PEListNodeListDataTransferRemoteResume;
        }
    } else {
        transferFunc = bRemote ? PEListPartListDataTransferRemote : PEListPartListDataTransferLocal;
        if (bResume) {
                transferFunc = PEListPartListDataTransferRemoteResume;
        }
    }

    transferFunc(request);
}

/// @brief Collect the interaction list results each a Compute operation completes
/// @param treePiece The TreePiece that sent the operation
/// @data A CudaRequest object containing the interaction list data
void PEList::sendList(TreePiece *treePiece, CudaRequest* data) {
    int numBucketsPlusOne = data->numBucketsPlusOne;
    int numBuckets = numBucketsPlusOne-1;

    // bucketMarkers need an offset because we are concatenating the interaction lists
    for (int i = 0; i < numBuckets; i++) {
	bucketMarkers.push_back(data->bucketMarkers[i] + iList.size());
    }
    finalBucketMarker = data->bucketMarkers[numBuckets] + iList.size();

    for (int i = 0; i < numBuckets; i++) {
	bucketStarts.push_back(data->bucketStarts[i]);
	bucketSizes.push_back(data->bucketSizes[i]);
    }

    // If we have missed parts/nodes, the indices in the interaction list
    // need to be shifted because the remote data is being concatenated
    if (data->missedParts) {
	// Note that many TreePieces will have the same missed particles
	// We are copying a lot of duplicate data to the GPU here
	int numMissedParts = data->sMissed/sizeof(CompactPartData);
	int missedOffset = missedParts.size();
	for (int i = 0; i < data->numInteractions; i++) {
	    ((ILCell *)data->list)[i].index += missedOffset;
	    iList.push_back(((ILCell *)data->list)[i]);
	}
	for (int i = 0; i < numMissedParts; i++) {
	    missedParts.push_back(((CompactPartData *)data->missedParts)[i]);
	}
    } else if (data->missedNodes) {
	int numMissedNodes = data->sMissed/sizeof(CudaMultipoleMoments);
	int missedOffset = missedNodes.size();
	for (int i = 0; i < data->numInteractions; i++) {
	    ((ILCell *)data->list)[i].index += missedOffset;
	    iList.push_back(((ILCell *)data->list)[i]);
	}
	for (int i = 0; i < numMissedNodes; i++) {
	    missedNodes.push_back(((CudaMultipoleMoments *)data->missedNodes)[i]);
	}
    } else {
	for (int i = 0; i < data->numInteractions; i++) {
	    iList.push_back(((ILCell *)data->list)[i]);
	}
    }

    fperiod = data->fperiod;

    // Call finishBucket for all buckets involved in this interaction
    treePiece->cudaFinishAffectedBuckets(data->affectedBuckets, numBuckets, bRemote);

    // deallocate the memory used by the incoming cudaRequest
    free(data->list);
    free(data->bucketMarkers);
    free(data->bucketStarts);
    free(data->bucketSizes);
    delete[] data->affectedBuckets;
}

/// @brief Re-initalize data arrays and clean up callback objects at the end of the step
void PEList::reset() {
    const char* funcTag = "PEList::reset";
    if (iList.size() > 0) {
      hostFree(iListHost, funcTag);
      hostFree(bucketMarkersHost, funcTag);
      hostFree(bucketStartsHost, funcTag);
      hostFree(bucketSizesHost, funcTag);
    }
    if (missedParts.size() > 0) {
      hostFree(missedPartsHost, funcTag);
    }
    if (missedNodes.size() > 0) {
      hostFree(missedNodesHost, funcTag);
    }
    iList.clear();
    missedParts.clear();
    missedNodes.clear();
    bucketMarkers.clear();
    bucketStarts.clear();
    bucketSizes.clear();

    cTreePieces.reset();
    vtpLocal.clear();
    bRemoteReady = 0;
    bKernelDelayed = 0;
    finalBucketMarker = -1;
    delete finishCb;
    delete request;
    request = nullptr;
}

#endif
