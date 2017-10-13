#!/bin/sh
# Use this for gasoline version after 6/06: which no longer supports
# unpacked ASCII vectors 
rm -f rdiff.acc lambs_30K.accg
echo Testing forces in gasoline
echo Expect RMS errors of .002
echo Expect Max errors of .09
../../charmrun ../../gasoline +p 4 lambs_30K.param >& DIAG
../../array/vecpackunpk < lambs_30K.accg > lambs_30K.up.accg
../../array/subarr lambs_30K.up.accg direct.acc > diff.acc
../../array/magvec < diff.acc > magdiff.arr
../../array/magvec < direct.acc > mag.acc
../../array/divarr magdiff.arr mag.acc > rdiff.acc
echo RMS relative force error:
../../array/rmsarr < rdiff.acc
echo Maximum relative force error:
../../array/maxarr < rdiff.acc
