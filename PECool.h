#ifndef PE_COOL_H
#define PE_COOL_H

#include "ParallelGravity.h"

class PECool : public CBase_PECool
{
  /// TreePieces on this PE
  CkVec<TreePiece*> vtpLocal;
  /// Count of TreePieces with particles on this PE
  NonEmptyTreePieceCounter cTreePieces;

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

public:
  PECool() { cudaStreamCreate(&stream); }
  PECool(CkMigrateMessage *m) : CBase_PECool(m) {}
  ~PECool() { cudaStreamDestroy(stream); }
  void pup(PUP::er &p) {}

  clDerivsData* getCoolData() { return coolData.data(); }
  double* getYInt() { return yInt.data(); }

  void finish(TreePiece *treePiece);
  int sendData(CoolRequest data);
  void reset();
};

#endif
