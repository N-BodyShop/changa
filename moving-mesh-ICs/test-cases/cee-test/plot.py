import yt

for i in range(1,21) : 
    filename = "cee.out.{0:06d}".format(i)
    dataset = yt.load(filename)
    yt.ProjectionPlot( dataset, "z", ("gas", "density"), width=500*7e10).save("frame{0:04d}.png".format(i))

