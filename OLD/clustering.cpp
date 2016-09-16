#include "SFC.h"

class KeySorted {
public:
  SFC::Key key;
  int index;
  inline bool operator<(const KeySorted& k) const {
    return key < k.key;
  }
};

void clustering(int n, double *weights, CkVec<TaggedVector3D> &centers, int procs, int *to) {
  int saved_peanoKey = peanoKey;
  peanoKey = 1;

  // create the SFC keys
  KeySorted *keys = new KeySorted[n];
  for (int i=0; i<n; ++i) {
    float centerArray[3];
    centers[i].vec.array_form(centerArray);
    
    keys[i].key = SFC::makeKey(centerArray);
    keys[i].index = centers[i].tag;
  }
  sort(keys, keys+n);

  // divide the keys in equal weighted segments
  double sum = 0;
  for (int i=0; i<n; ++i) sum += weights[i];
  CkAssert(sum > 0.0);
  double increment = sum / procs;
  CkPrintf("sum: %f, increment: %f\n", sum, increment);
  double target = increment;
  double partial = 0;
  int cluster = 0;
  for (int i=0; i<n; ++i) {
    partial += weights[keys[i].index];
    to[keys[i].index] = cluster;
    if (partial >= target) {
      if (2*(partial-target) > weights[keys[i].index]) to[keys[i].index] ++;
      cluster ++;
      target += increment;
    }
    CkAssert(to[keys[i].index] < procs);
  }

  peanoKey = saved_peanoKey;
  delete [] keys;
}
