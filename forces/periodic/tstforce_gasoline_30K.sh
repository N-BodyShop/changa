#!/bin/bash
echo Testing periodic forces in gasoline
echo Expect RMS force errors of .004
echo Expect Max force errors of .11
../gasoline lambs_30K.param >& DIAG
../../array/subarr lambs_30K.accg direct.acc > diff.acc
../../array/magvec < diff.acc > magdiff.arr
../../array/magvec < direct.acc > mag.acc
../../array/divarr magdiff.arr mag.acc > rdiff.acc
echo RMS relative force error:
../../array/rmsarr < rdiff.acc
echo Maximum relative force error:
../../array/maxarr < rdiff.acc
