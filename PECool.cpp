#include "PECool.h"
#include "HostCUDA.h"

void PECool::finish(TreePiece *treePiece) {
    vtpLocal.push_back(treePiece);

    // On first call, find the total number of active pieces on this PE.
    // The charm++ location manager gives us this count in cTreePieces
    if(vtpLocal.length() == 1) {
        CkLocMgr *locMgr = treeProxy.ckLocMgr();
        locMgr->iterate(cTreePieces);
    }

    // check if we have everyone
    if(vtpLocal.length() < cTreePieces.count)
        return;

    // PE with no particles
    if (coolData.empty()) {
        treeProxy.finishIntegrateCb();
	return;
    }
    int numParts = coolData.size();
    int nv = stiff[0].nv;
    cudaChk(cudaMalloc(&d_CoolData, numParts*sizeof(clDerivsData)));
    cudaChk(cudaMalloc(&d_Stiff, numParts*sizeof(STIFF)));
    cudaChk(cudaMalloc(&d_dtg, numParts*sizeof(double)));

    cudaChk(cudaMalloc(&d_y, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_ymin, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_y0, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_y1, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_q, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_d, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_rtau, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_ys, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_qs, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_rtaus, numParts*nv*sizeof(double)));
    cudaChk(cudaMalloc(&d_scrarray, numParts*nv*sizeof(double)));

    for (int i = 0; i < coolData.size(); i++) {
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

    cudaChk(cudaMemcpy(d_CoolData, coolData.data(), numParts*sizeof(clDerivsData), cudaMemcpyHostToDevice));
    cudaChk(cudaMemcpy(d_Stiff, stiff.data(), sizeof(STIFF)*numParts, cudaMemcpyHostToDevice));
    cudaChk(cudaMemcpy(d_ymin, yMin.data(), numParts*nv*sizeof(double), cudaMemcpyHostToDevice));
    cudaChk(cudaMemcpy(d_y, yInt.data(), numParts*nv*sizeof(double), cudaMemcpyHostToDevice));
    cudaChk(cudaMemcpy(d_dtg, dtg.data(), dtg.size()*sizeof(double), cudaMemcpyHostToDevice));

    double t = 0.0;
    PeODESolver(d_Stiff, d_y, d_dtg, t, numParts, stream);

    cudaChk(cudaMemcpy(yInt.data(), d_y, numParts*nv*sizeof(double), cudaMemcpyDeviceToHost));
    cudaChk(cudaMemcpy(coolData.data(), d_CoolData, numParts*sizeof(clDerivsData), cudaMemcpyDeviceToHost));

    // TODO hapi callback
    treeProxy.finishIntegrateCb();
}

int PECool::sendData(CoolRequest data) {
  d_Cool = data.d_Cool;

  coolData.push_back(*data.coolData);
  stiff.push_back(*data.coolData->IntegratorContext);
  int nv = data.coolData->IntegratorContext->nv;
  for (int i = 0; i < nv; i++) {
    yMin.push_back(data.coolData->IntegratorContext->ymin[i]);
    yInt.push_back(data.y[i]);
  }

  dtg.push_back(data.dtg);

  return coolData.size()-1;
}

void PECool::reset() {
    coolData.clear();
    stiff.clear();
    yMin.clear();
    yInt.clear();
    dtg.clear();

    cTreePieces.reset();
    vtpLocal.length() = 0;
}
