grep 'feedback 0:' onestar.out | awk '{print $3,$4,$5}' > SNII.out
grep 'feedback 1:' onestar.out | awk '{print $3,$4,$5}' > SNIa.out
grep 'feedback 2:' onestar.out | awk '{print $3,$4,$5}' > winds.out

grep dDeltaStarForm onestar.log > eff.out
grep dKpcUnit onestar.param >> eff.out
grep dMsolUnit onestar.param >> eff.out

idl idl.metals.batch

