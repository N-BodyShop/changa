#!/bin/sh
rm -f pkdden.den
# gasoline density.param
charmrun +p 16 gasoline density.param
#
# difference with the standard
# 
echo EXPECT about 7e-7
../array/subarr pkdden.den bench.den > diff.den
../array/divarr diff.den bench.den | ../array/maxarr
