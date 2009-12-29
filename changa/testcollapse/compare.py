def column(filename, n) :
    """Read a column from a file into a list."""
    f = file(filename)
    s = f.readline()
    x = list()
    while s:
        if not s.startswith('#') :
            l = s.split()
            x.append(float(l[n]))
        s = f.readline()
    return x

r = column('gas.prof', 0)
rho = column('gas.prof', 12)
pres = column('gas.prof', 14)
vr = column('gas.prof', 5)
rc = column('cha.prof', 0)
rhoc = column('cha.prof', 12)
presc = column('cha.prof', 14)
vrc = column('cha.prof', 5)
s = [pres[i]/rho[i]**(5./3.) for i in range(len(pres))]
sc = [presc[i]/rhoc[i]**(5./3.) for i in range(len(presc))]
# compare inner entropy profile
rds = [abs((sc[i] - s[i])/s[i]) for i in range(len(sc)) if r[i] < .1]
print 'Expect relative differences of less than 0.04'
print 'maximum relative entropy difference in inner region: ', max(rds)
rdvr = [abs((vrc[i] - vr[i])/vr[i]) for i in range(len(vrc)) if r[i] > .25]
print 'maximum relative velocity difference in outer region: ', max(rdvr)

