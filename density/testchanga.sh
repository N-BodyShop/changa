#!/bin/sh
rm -f lambs.00200_subsamp_30K.den
charmrun +p 1 ChaNGa -z -v 1 density.param
#
# difference with the standard
# 
echo EXPECT about 7e-7
../array/subarr lambs.00200_subsamp_30K.den bench.den > diff.den
../array/divarr diff.den bench.den | ../array/absarr | ../array/maxarr
