#ifdef CUDACOOL
#include "PECool.h"
#include "HostCUDA.h"

/// @brief Each TreePiece on this PE checks in once it completes updateuDot
///        Once all TreePieces are done, copy the particle cooling data to
///        the GPU and launch StiffStep
/// @param treePiece A reference to the treePiece that checked in
void PECool::finish(TreePiece *treePiece) {
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

    // PE with no particles
    if (coolData.empty()) {
	finishIntegrateCb();
	return;
    }
    int numParts = coolData.size();
    int nv = stiff[0].nv;
    cudaChk(cudaMallocAsync(&d_CoolData, numParts*sizeof(clDerivsData), stream));
    cudaChk(cudaMallocAsync(&d_Stiff, numParts*sizeof(STIFF), stream));
    cudaChk(cudaMallocAsync(&d_dtg, numParts*sizeof(double), stream));

    cudaChk(cudaMallocAsync(&d_y, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_ymin, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_y0, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_y1, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_q, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_d, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_rtau, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_ys, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_qs, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_rtaus, numParts*nv*sizeof(double), stream));
    cudaChk(cudaMallocAsync(&d_scrarray, numParts*nv*sizeof(double), stream));

    // Deep copy each IntegratorContext
    for (int i = 0; i < numParts; i++) {
        coolData[i].IntegratorContext = &d_Stiff[i];

        stiff[i].ymin = &d_ymin[i*nv];
        stiff[i].y0 = &d_y0[i*nv];
        stiff[i].y1 = &d_y1[i*nv];
        stiff[i].q = &d_q[i*nv];
        stiff[i].d = &d_d[i*nv];
        stiff[i].rtau = &d_rtau[i*nv];
        stiff[i].ys = &d_ys[i*nv];
        stiff[i].qs = &d_qs[i*nv];
        stiff[i].rtaus = &d_rtaus[i*nv];
        stiff[i].scrarray = &d_scrarray[i*nv];
	stiff[i].Data = &d_CoolData[i];
    }

    cudaChk(cudaMemcpyAsync(d_CoolData, coolData.data(), numParts*sizeof(clDerivsData), cudaMemcpyHostToDevice, stream));
    cudaChk(cudaMemcpyAsync(d_Stiff, stiff.data(), sizeof(STIFF)*numParts, cudaMemcpyHostToDevice, stream));
    cudaChk(cudaMemcpyAsync(d_ymin, yMin.data(), numParts*nv*sizeof(double), cudaMemcpyHostToDevice, stream));
    cudaChk(cudaMemcpyAsync(d_y, yInt.data(), numParts*nv*sizeof(double), cudaMemcpyHostToDevice, stream));
    cudaChk(cudaMemcpyAsync(d_dtg, dtg.data(), dtg.size()*sizeof(double), cudaMemcpyHostToDevice, stream));

    double t = 0.0;
    PeODESolver(d_Stiff, d_y, d_dtg, t, numParts, stream);

    cudaChk(cudaMemcpyAsync(yInt.data(), d_y, numParts*nv*sizeof(double), cudaMemcpyDeviceToHost, stream));
    cudaChk(cudaMemcpyAsync(coolData.data(), d_CoolData, numParts*sizeof(clDerivsData), cudaMemcpyDeviceToHost, stream));

    finishCb = new CkCallback(CkIndex_PECool::finishIntegrateCb(), CkMyPe(), thisProxy);
    hapiAddCallback(stream, finishCb);
}

/// @brief Callback that triggers after the StiffStep kernel finishes. Calls finishuDot
///        every TreePiece on this PE
void PECool::finishIntegrateCb() {
    for (int i = 0; i < vtpLocal.size(); i++) {
       vtpLocal[i]->finishuDot();
    }
}

/// @brief Used by updateuDot to accumulate cooling data for individual particles
/// @param data A struct containing the data to eventually be shipped to the GPU
/// @return Current number of particles accumulated on this PE. Used to map TreePiece
///         arrays onto the PE arrays
int PECool::sendData(CoolRequest data) {
  d_Cool = data.d_Cool;

  coolData.push_back(*data.coolData);
  stiff.push_back(*(data.coolData->IntegratorContext));
  int nv = data.coolData->IntegratorContext->nv;
  for (int i = 0; i < nv; i++) {
    yMin.push_back(data.coolData->IntegratorContext->ymin[i]);
    yInt.push_back(data.y[i]);
  }

  dtg.push_back(data.dtg);

  return coolData.size()-1;
}

/// @brief Wipes the accumulated PE-level data for the next cooling calculation
void PECool::reset() {
    // This function is called by every TreePiece
    // Dont reset until the final TreePiece checks in
    treePiecesDone++;
    if (treePiecesDone == vtpLocal.size()) {
      treePiecesDone = 0;
      if (finishCb != nullptr) {
        delete finishCb;
        finishCb = nullptr;
      }
      if (!coolData.empty()) {
        cudaChk(cudaFree(d_CoolData));
        cudaChk(cudaFree(d_Stiff));
        cudaChk(cudaFree(d_dtg));

        cudaChk(cudaFree(d_y));
        cudaChk(cudaFree(d_ymin));
        cudaChk(cudaFree(d_y0));
        cudaChk(cudaFree(d_y1));
        cudaChk(cudaFree(d_q));
        cudaChk(cudaFree(d_d));
        cudaChk(cudaFree(d_rtau));
        cudaChk(cudaFree(d_ys));
        cudaChk(cudaFree(d_qs));
        cudaChk(cudaFree(d_rtaus));
        cudaChk(cudaFree(d_scrarray));
      }
      coolData.clear();
      stiff.clear();
      yMin.clear();
      yInt.clear();
      dtg.clear();

      cTreePieces.reset();
      vtpLocal.clear();
    }
}

#endif
