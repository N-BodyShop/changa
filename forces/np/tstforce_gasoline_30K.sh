#!/bin/sh
echo Testing forces in gasoline
echo Expect RMS errors of .002
echo Expect Max errors of .09
../gasoline lambs_30K.param >& DIAG
../../array/subarr lambs_30K.accg direct.acc > diff.acc
../../array/magvec < diff.acc > magdiff.arr
../../array/magvec < direct.acc > mag.acc
../../array/divarr magdiff.arr mag.acc > rdiff.acc
echo RMS relative force error:
../../array/rmsarr < rdiff.acc
echo Maximum relative force error:
../../array/maxarr < rdiff.acc
