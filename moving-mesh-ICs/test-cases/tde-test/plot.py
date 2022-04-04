import yt

for i in range(1,21) : 
    filename = "tde.out.{0:06d}".format(i)
    dataset = yt.load(filename)
    yt.ProjectionPlot( dataset, "y", ("gas", "density"), width=2000*7e10).save("frame{0:04d}.png".format(i))

