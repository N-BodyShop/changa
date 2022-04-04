import matplotlib
matplotlib.use("Agg")
import yt
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
import argparse

fig = pl.figure()

bound = 2e12
tdyn = 6.7e-8**-.5
Rsun = 7e10

grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 2),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")

def plotfile(fig,grid,i,filename, title=None) :
        ds = yt.load( filename, n_ref=8)
	p = yt.SlicePlot( ds, "z", ("gas", "density"), center="c")
	#p.set_zlim(("gas", "density"), 1e-25, 1e-28)
	if( title != None) : 
		p.annotate_text((3*Rsun,4.*Rsun,0.), title, coord_system='data',text_args={'color':'white'})
	plot = p.plots[("gas", "density")]
	plot.figure = fig
	plot.axes=grid[i].axes
	plot.cax=grid.cbar_axes[i]
	p._setup_plots()
	

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('start', metavar='N', type=int,
                    help='an integer for the accumulator')
parser.add_argument('end', metavar='N', type=int, 
                    help='an integer for the accumulator')
parser.add_argument('step', metavar='N', type=int, 
                    help='an integer for the accumulator')
parser.add_argument('--framestep', default=1, metavar='N', type=int, 
                    help='an integer for the accumulator')

args = parser.parse_args()

for i in range(args.start, args.end, args.step) :
        plotfile(fig,grid,1, "mm/star.out.{0:06d}".format(i), title="MM ".format(i))
        #plotfile(fig,grid,0, "sm/star.out.{0:06d}".format(i), title="Entropy ".format(i))
        plotfile(fig,grid,0, "sph/star.out.{0:06d}".format(i), title="SPH ")

        pl.savefig("movie_frame{0:04d}.png".format(i/args.framestep))

