#!/bin/sh
#
#
CHARMRUN=../charmrun
CHANGA=../ChaNGa
#
echo "Beginning ChaNGa"
$CHARMRUN $CHANGA onestar.param > onestar.out
echo "ChaNGa finished"
