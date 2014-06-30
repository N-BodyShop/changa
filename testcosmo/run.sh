#!/bin/sh
#
#
CHARMRUN=../charmrun
CHANGA=../ChaNGa
#
echo "Beginning ChaNGa"
$CHARMRUN +p8 $CHANGA cube300.param > cube300.out
echo "ChaNGa finished"
