idl idl.imf.batch

echo IGNORE WARNING ABOUT COOLING EOS TURNED OFF

pkdgravdir/gasoline onestar.param > onestar.out

rm onestar.000*
rm onestar.001*
rm onestar.002*

grep 'feedback 0:' onestar.out | awk '{print $3,$4,$5}' > SNII.out
grep 'feedback 1:' onestar.out | awk '{print $3,$4,$5}' > SNIa.out
grep 'feedback 2:' onestar.out | awk '{print $3,$4,$5}' > winds.out

grep effectively onestar.log > eff.out

idl idl.metals.batch

