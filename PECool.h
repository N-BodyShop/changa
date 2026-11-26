#ifndef PE_COOL_H
#define PE_COOL_H

#ifdef CUDACOOL

#include "ParallelGravity.h"

/// @brief PE-level group to coordinate cooling calculations for gas particles
/// on the GPU
///
/// TreePieces call updateuDot, which loops over individual particles. This loop
/// calls sendData, which accumulates IntegratorContext and clDerivsData for each
/// particle. One all TreePieces have checked in, the finish routine copies this
/// data to device memory and launches a GPU version of StiffStep. A callback
/// function finishIntegrateCb then updates the internal energies of the particles
/// on each TreePiece accordingly.
class PECool : public CBase_PECool
{
  /// TreePieces on this PE
  CkVec<TreePiece*> vtpLocal;
  /// Count of TreePieces with particles on this PE
  NonEmptyGasTreePieceCounter cTreePieces;

  /// Counter to ensure all TreePieces have finished
  /// before clearing PE-level data
  int treePiecesDone;

  // vectors for host data
  vector<clDerivsData> coolData;
  vector<STIFF> stiff;
  vector<double> yMin;
  vector<double> yInt;
  vector<double> dtg;

  // Device pointers for data
  COOL *d_Cool;
  clDerivsData *d_CoolData;
  STIFF *d_Stiff;
  double *d_dtg;
  double *d_ymin;
  double *d_y;
  double *d_y0;
  double *d_y1;
  double *d_q;
  double *d_d;
  double *d_rtau;
  double *d_ys;
  double *d_qs;
  double *d_rtaus;
  double *d_scrarray;

  cudaStream_t stream;
  CkCallback *finishCb;

public:
  PECool() { finishCb=nullptr; treePiecesDone=0; cudaStreamCreate(&stream); }
  PECool(CkMigrateMessage *m) : CBase_PECool(m) {}
  ~PECool() { cudaStreamDestroy(stream); }
  void pup(PUP::er &p) {}

  clDerivsData* getCoolData() { return coolData.data(); }
  double* getYInt() { return yInt.data(); }

  void finish(TreePiece *treePiece);
  void finishIntegrateCb();
  int sendData(CoolRequest data);
  void reset();
};

#endif
#endif
