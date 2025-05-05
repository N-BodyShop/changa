#include "PEList.h"
#include "HostCUDA.h"

void PEList::setType(int _bNode, int _bRemote, int _bResume) {
        bNode = _bNode;
        bRemote = _bRemote;
        bResume = _bResume;
    }

void PEList::finishWalk(TreePiece *treePiece, CkCallback &cb) {
    vtpLocal.push_back(treePiece);

    // On first call, find the total number of active pieces on this PE.
    // The charm++ location manager gives us this count in cTreePieces
    if(vtpLocal.length() == 1) {
        CkLocMgr *locMgr = treeProxy.ckLocMgr();
        locMgr->iterate(cTreePieces);
    }
    // check if we have everyone
    //CkPrintf("(%d %d %d) %d of %d\n", bNode, bRemote, bResume, vtpLocal.length(), cTreePieces.count);
    if(vtpLocal.length() < cTreePieces.count)
        return;

    //CkPrintf("(%d) finishWalk %d %d %d\n", CmiMyPe(), bNode, bRemote, bResume);

    // bucketMarkers[i+1] is needed to determine # of IL entries per bucket
    if(vtpLocal.length() == cTreePieces.count && finalBucketMarker != -1)
	bucketMarkers.push_back(finalBucketMarker);

    CudaRequest *request = new CudaRequest;

    request->d_localMoments = d_localMoments;
    request->d_localParts = d_localParts;
    request->d_localVars = d_localVars;
    request->d_remoteParts = d_remoteParts;
    request->d_remoteMoments = d_remoteMoments;
    request->stream = stream;

    request->numBucketsPlusOne = bucketSizes.size()+1;

    request->node = bNode;
    request->remote = bRemote;

    request->fperiod = fperiod;

    request->list = iList.data();
    request->missedParts = missedParts.data();
    request->missedNodes = missedNodes.data();
    request->sMissed = missedParts.size() > 0 ? missedParts.size()*sizeof(CompactPartData) : missedNodes.size()*sizeof(CudaMultipoleMoments);
    request->bucketMarkers = bucketMarkers.data();
    request->bucketStarts = bucketStarts.data();
    request->bucketSizes = bucketSizes.data();
    request->numInteractions = iList.size();

    CProxy_PEList proxy;
    void (*transferFunc)(CudaRequest*);
    if (bNode) {
	proxy = bRemote ? peNodeRemoteListProxy : peNodeLocalListProxy;
	transferFunc = bRemote ? PEListNodeListDataTransferRemote : PEListNodeListDataTransferLocal;
	if (bResume) {
		proxy = peNodeRemoteResumeListProxy;
		transferFunc = PEListNodeListDataTransferRemoteResume;
	}
    } else {
	proxy = bRemote ? pePartRemoteListProxy : pePartLocalListProxy;
	transferFunc = bRemote ? PEListPartListDataTransferRemote : PEListPartListDataTransferLocal;
	if (bResume) {
		proxy = pePartRemoteResumeListProxy;
		transferFunc = PEListPartListDataTransferRemoteResume;
	}
    }

    request->callback = &cb;
    transferFunc(request);
}

void PEList::sendList(TreePiece *treePiece, CudaRequest* data) {
    // list, bucketMarkers, bucketStarts and bucketSizes need to be collected into one array
    // note that bucketMarkers will need to be updated
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

    d_localMoments = data->d_localMoments;
    d_localParts = data->d_localParts;
    d_localVars = data->d_localVars;
    d_remoteParts = data->d_remoteParts;
    d_remoteMoments = data->d_remoteMoments;
    fperiod = data->fperiod;

    treePiece->cudaFinishAffectedBuckets(data->affectedBuckets, numBuckets, bRemote);

    // deallocate the memory used by the incoming cudaRequest
    freePinnedHostMemory(data->list);
    freePinnedHostMemory(data->bucketMarkers);
    freePinnedHostMemory(data->bucketStarts);
    freePinnedHostMemory(data->bucketSizes);
    delete[] data->affectedBuckets;
}

void PEList::reset() {
    iList.clear();
    missedParts.clear();
    missedNodes.clear();
    bucketMarkers.clear();
    bucketStarts.clear();
    bucketSizes.clear();

    cTreePieces.reset();
    vtpLocal.length() = 0;
    finalBucketMarker = -1;
}
