#!/bin/sh
#
# edit these
#
PKDGRAVDIR=..
SKIDDIR=$HOME/peak/src/skid
SODIR=$HOME/peak/src/so
#
# end of edits
#
echo "Beginning pkdgrav"
../charmrun +p 4 $PKDGRAVDIR/gasoline cube300.param > cube300.out
echo "pkdgrav finished"
echo "Moving to SKID"
$SKIDDIR/skid -tau 0.005  -std -O 0.3 -H 2.894405 -s 12 -m 16 -o cube300.00128.skid -stats <cube300.00128
echo "Moving to SO"
$SODIR/so -i cube300.00128.skid.gtp  -o cube300.00128.so -all -std -grp -gtp -stat cube300.00128.skid.stat -m 16 -O 0.3 -L  -p 1 -u 3.672e18 300 < cube300.00128
echo "SO complete"
echo "Moving to Mass Function comparison"
idl idl.mf.batch
