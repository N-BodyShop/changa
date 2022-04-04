import yt
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
G = 6.67e-8

tOrig, sepOrig = np.loadtxt("orbit.data", unpack=True)

tArray = []
sepArray = []
for i in range(1,21) : 
    filename = "cee.out.{0:06d}".format(i)
    dataset = yt.load(filename)
    ad = dataset.all_data()
    dm1, dm2 = ad[("DarkMatter", "Coordinates")].v
    sep = np.linalg.norm( dm1-dm2)/7e10
    t = dataset.current_time.v*G**-0.5/3600/24
    tArray.append(t)
    sepArray.append(sep)

pl.plot(tOrig, sepOrig, ls="dotted", lw=2, label="baseline")
pl.plot(tArray, sepArray, ls="solid", lw=2, label="current run")
pl.legend(loc="best")
pl.savefig("orbit.pdf")
