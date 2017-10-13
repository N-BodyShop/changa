#!/bin/sh
# Non-periodic version of "test4", a cluster.
TESTDIR=../..

echo "# theta a_rms a_max time" > acc_cha_test4np.out
for theta in .3 .4 .5 .6 .7 .8 .9 1.0 ; do

#$TESTDIR/ChaNGa -consph -theta $theta test4np_changa.param >& DIAG
$TESTDIR/ChaNGa -theta $theta test4np_changa.param >& DIAG
$TESTDIR/array/subarr test4.000000.acc2 ../np/test4np.direct.acc > diff.acc
$TESTDIR/array/magvec < diff.acc > magdiff.arr
$TESTDIR/array/magvec < ../np/test4np.direct.acc > mag.acc
$TESTDIR/array/divarr magdiff.arr mag.acc > rdiff.acc
RMS=`$TESTDIR/array/rmsarr < rdiff.acc`
MAX=`$TESTDIR/array/maxarr < rdiff.acc`
TIME=`grep Calc DIAG | awk '{print $8;}'`
# TIME=`grep Calc DIAG | awk '{print $12;}'`
# FLOPS=0
# NPARTS=`grep particle-particle DIAG | awk '{print $2/30000.;}'`
# NCELLS=`grep particle-node DIAG | awk '{print $2/30000.;}'`
echo $theta $RMS $MAX $TIME >> acc_cha_test4np.out

done
