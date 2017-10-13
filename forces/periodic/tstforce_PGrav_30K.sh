# get rid of old output
rm ../lambs.00200_subsamp_30K.acc2
# run program
../../charmrun +p 2 ../../ParallelGravity +f -v 1 -D 1 -c 1 -p 4 -z lambs_30K.pgrav.param >& DIAG

../../array/sortvec ../lambs.00200_subsamp_30K.acc2 ../lambs.00200_subsamp_30K.iord > sorted.acc
../../array/subarr sorted.acc direct.acc > diff.acc
../../array/magvec < diff.acc > magdiff.arr
../../array/magvec < direct.acc > mag.acc
../../array/divarr magdiff.arr mag.acc > rdiff.acc
echo RMS relative force error:
../../array/rmsarr < rdiff.acc
echo Maximum relative force error:
../../array/maxarr < rdiff.acc
