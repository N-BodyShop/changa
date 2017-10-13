#!/bin/sh
TESTDIR=../..

echo "# theta a_rms a_max time" > acc_cha_lambs.out
for theta in .3 .4 .5 .6 .7 .8 .9 1.0 ; do

NREPS=1
if [ $theta = .3 ] ; then
        NREPS=3
fi
if [ $theta = .4 ] ; then
        NREPS=2
fi
if [ $theta = .5 ] ; then
        NREPS=2
fi

./ChaNGa -nrep $NREPS -theta $theta lambs_30K_changa.param >& DIAG
$TESTDIR/array/subarr lambs_30K.000000.acc2 ../periodic/direct.acc > diff.acc
$TESTDIR/array/magvec < diff.acc > magdiff.arr
$TESTDIR/array/magvec < ../periodic/direct.acc > mag.acc
$TESTDIR/array/divarr magdiff.arr mag.acc > rdiff.acc
RMS=`$TESTDIR/array/rmsarr < rdiff.acc`
MAX=`$TESTDIR/array/maxarr < rdiff.acc`
# TIME=`grep Calc DIAG | awk '{print $8;}'`
TIME=`grep Calc DIAG | awk '{print $12;}'`
# FLOPS=0
# NPARTS=`grep particle-particle DIAG | awk '{print $2/30000.;}'`
# NCELLS=`grep particle-node DIAG | awk '{print $2/30000.;}'`
echo $theta $RMS $MAX $TIME >> acc_cha_lambs.out

done
