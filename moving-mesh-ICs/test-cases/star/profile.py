import matplotlib
matplotlib.use("Agg")
import yt
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
import argparse

fig = None
grid = None
#fig = pl.figure()

bound = 2e12
tdyn = 6.7e-8**-.5
Rsun = 7e10

#grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
#                nrows_ncols = (1, 3),
#                axes_pad = 0.05,
#                label_mode = "L",
#                share_all = True,
#                cbar_location="right",
#                cbar_mode="single",
#                cbar_size="3%",
#                cbar_pad="0%")

def plotfile(filename) :
        ds = yt.load( filename, n_ref=8)
        sphere = ds.sphere( [0.,0.,0.], (2e11, "cm"))
        profile = yt.create_profile( sphere, "radius", ("gas", "density"), weight_field="cell_volume")
	r = profile.x
	rho = profile[("gas", "density")]
	time = ds.current_time*tdyn
	
	return time, r, rho
	#if( title != None) : 
#		p.annotate_text((7.5*Rsun, 1e1), title, coord_system='data',text_args={'color':'black'})
#	plot = p.plots[("gas", "density")]
#	plot.figure = fig
#	plot.axes=grid[i].axes
#	plot.cax=grid.cbar_axes[i]
#	p._setup_plots()
	

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
	pl.clf()
	labels = []
	radii = []
	rhos = []
        t, r, rho = plotfile("mm/star.out.{0:06d}".format(i))
	radii.append(r)
	rhos.append(rho)
        labels.append( "MM {0}".format(t/60))
        #t, r, rho = plotfile("sm/star.out.{0:06d}".format(i))
	#radii.append(r)
	#rhos.append(rho)
        #labels.append( "Entropy {0}".format(t/60))
        #plotfile(fig,grid,0, "sm/test_out.{0:06d}".format(i), title="SM ".format(i))
        #plotfile(fig,grid,1, "sph/star.out.{0:06d}".format(i), title="SPH ")

        t, r, rho = plotfile("sph/star.out.{0:06d}".format(i))
	radii.append(r)
	rhos.append(rho)
        labels.append( "SPH")
	for label, r, rho in zip( labels, radii, rhos) : 
		pl.loglog( r, rho, label=label, lw=2)
	pl.legend(loc="best",fontsize=16)
	pl.xlim(3e9,3e11)
	pl.ylim(1e-6,1e2)
	pl.xlabel("$r$ [cm]", fontsize=20)
	pl.xlabel("$\\rho$ [g/cc]", fontsize=20)
	pl.xticks(fontsize=18)
	pl.yticks(fontsize=18)
        pl.savefig("movie_frame{0:04d}.png".format(i/args.framestep))
	

