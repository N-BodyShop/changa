#!/bin/sh
#
#
CHARMRUN=../charmrun
CHANGA=../ChaNGa
#
echo "Beginning ChaNGa"
$CHARMRUN $CHANGA +LBPeriod 0.01 onestar.param > onestar.out
echo "ChaNGa finished"
