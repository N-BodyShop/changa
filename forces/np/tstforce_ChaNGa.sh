#!/bin/sh
# get rid of old output
#rm lambs_30K.000000.acc2
echo Testing forces in ChaNGa
echo For HEXADECAPOLE:
echo Expect RMS errors of .0008
echo Expect Max errors of .03
# run program
#../../charmrun +p 2 ../../ChaNGa -v 1 lambs_30K.param >& DIAG
../../array/subarr lambs_30K.000000.acc2 direct.acc > diff.acc
../../array/magvec < diff.acc > magdiff.arr
../../array/magvec < direct.acc > mag.acc
../../array/divarr magdiff.arr mag.acc > rdiff.acc
echo RMS relative force error:
../../array/rmsarr < rdiff.acc
echo Maximum relative force error:
../../array/maxarr < rdiff.acc
