#!/bin/sh
#
# edit these
#
SKIDDIR=XXX
SODIR=XXX
#
# end of edits
#
echo "SKIDDING outputs"
$SKIDDIR/skid -tau 0.005  -std -O 0.3 -H 2.894405 -s 12 -m 16 -o cube300.000128.skid -stats <cube300.000128
echo "Moving to SO"
$SODIR/so -i cube300.000128.skid.gtp  -o cube300.000128.so -all -std -grp -gtp -stat cube300.000128.skid.stat -m 16 -O 0.3 -L  -p 1 -u 3.672e18 300 < cube300.000128
echo "SO complete"
echo "Moving to Mass Function comparison"
python3 mf2.py cube300.000.so.sovcirc
