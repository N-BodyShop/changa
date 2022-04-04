import mesa_reader as mr
import matplotlib
matplotlib.use("Agg")
import math
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from tqdm import tqdm
import io
from scipy.interpolate import interp1d

import sys
sys.path.append('..')

from readWriteTipsy import *

import argparse

mSun = 1.99e33
rSun = 7e10
rho_max = 1.
rotSun = 2.0*math.pi/25.0/24.0/3600.0/math.sqrt(6.674e-8)

def ScaleFactor( rParticle, rStar, rMax=1000.) : 
    scaleR = 2**np.minimum(np.maximum((rParticle/rStar - 1.)*6,0),np.log2(rMax))
    scaleR = np.maximum( scaleR, 1.)
    scaleR = np.minimum( scaleR, rMax)
    return scaleR

# def ScaleFactor( radius, rStar, maxRatio=100 ) :
# 	scale = (radius / rStar)**3
# 	scale = np.maximum( scale, 1.)
# 	scale = np.minimum( scale, maxRatio)
# 	return scale

def filterMonotonicity( rho) :
	boolArray = np.ones(rho.size, dtype=bool)
	return boolArray
	mono = True
	for i in range(rho.size-1) :
		if( mono and rho[i] > rho[i+1]) :
			boolArray[i] = True
		else :
			mono = False
			boolArray[i] = False

	boolArray[-1] = boolArray[-2]
	return boolArray

def makeStarMesh(profileFile, massParticleInGram, posArray, sphICs=False) :
	mcore, rcore, mModel, rModel, rhoModel, eModel, XModel = buildStar( profileFile, True)

	mcore = mcore*mSun
	rcore = rcore*rSun
	mModel = mModel*mSun
	rModel = rModel*rSun
	mass=mModel.max()
	rho_max = rhoModel.max()

	massParticleInCode = massParticleInGram
	numParticles = int(mass/massParticleInGram)
	rArray = np.linalg.norm( posArray, axis=1)
	rCircle = rArray.max()/math.sqrt(3.)
	boolArray = np.argsort(rArray)
	rArray = rArray[boolArray]
	posArray = posArray[boolArray]
	eArray = np.ones(rArray.size)
	vArray = np.ones([rArray.size, 3])
	XArray = np.zeros(rArray.size)
	rMax = None
	mArray = np.ones(rArray.size)
	lastModelP = 0

	if( not sphICs) :
		# now rescale
		rMax = rModel.max()
		rArray = rArray*rMax/rCircle*(rArray.size/numParticles)**0.33333
		interpolate = interp1d(rModel,rhoModel,kind = 'linear')
		interpolateE = interp1d(rModel,eModel,kind = 'linear')
		interpolateX = interp1d(rModel,XModel,kind = 'linear')
		lastModelP = rArray[rArray <= rMax].size
		mArray[:lastModelP] = interpolate( rArray[:lastModelP])
		eArray[:lastModelP] = interpolateE( rArray[:lastModelP])
		XArray[:lastModelP] = interpolateX( rArray[:lastModelP])
		posArray[:lastModelP,:] = posArray[:lastModelP,:]/np.linalg.norm(posArray[:lastModelP,:],axis=1)[:,np.newaxis]*rArray[:lastModelP,np.newaxis]
		mArrayOut = mArray[:lastModelP] * mass / mArray[:lastModelP].sum()
	else : 
		mArray = np.arange(1,rArray.size)*massParticleInGram
		interpolate = interp1d(mModel,rModel,kind = 'linear')
		interpolateE = interp1d(mModel,eModel,kind = 'linear')
		interpolateX = interp1d(mModel,XModel,kind = 'linear')

		lastModelP = 0

		for i in range(numParticles) :
			if( mArray[i] > mModel[-1] ) :
				break
			lastModelP = i
			r = interpolate( mArray[i])
			eArray[i] = interpolateE( mArray[i])
			XArray[i] = interpolateX( mArray[i])
			posArray[i,:] = posArray[i,:]/rArray[i] * r

		lastModelP = lastModelP+1
		mArrayOut = massParticleInCode*np.ones(lastModelP)

	return posArray[:lastModelP], vArray[:lastModelP], eArray[:lastModelP], XArray[:lastModelP], mArrayOut, rho_max, mcore, rcore, mass

def buildStar( profileFile="profile15.data", makePlot=False) :
	data = mr.MesaData(profileFile)

	logT   = data.logT
	mass   = data.mass
	logR   = data.logR
	logRho = data.logRho
	logP   = data.logP
	Hfrac  = data.x_mass_fraction_H
	dlogTdlogP = (logT[1:] - logT[:-1])/(logP[1:] - logP[:-1])
	dlogTdlogP = np.append( np.array([dlogTdlogP[0]]), dlogTdlogP)
	R = 1e1**logR

	mtot = mass.max()
	print( 'mtot = ' + str(mtot))
	maxR = R.max()

	rho = 1e1**logRho
	if( makePlot ) :
		plt.clf()
		plt.semilogy(R, 1e1**logRho, lw=2, color="blue")

	meanRho = mtot*mSun/(4.*math.pi/3.*(maxR*rSun)**3)
	print( 'mean density = ' + str(meanRho))
	cutoffRho = 500.0*meanRho
	# cutoffRho = 0.05
	print( 'cutoff density = ' + str(cutoffRho))

	boolArray = rho<cutoffRho
	nboolArray = np.logical_not(boolArray)
	massCentral = 0.
	logRCentral = -9.
	if( any(nboolArray)) :
		massCentral = mass[nboolArray].max()
		logRCentral = logR[nboolArray].max()

	mass = mass[boolArray]
	logR = logR[boolArray]
	logT = logT[boolArray]
	logRho = logRho[boolArray]
	logP = logP[boolArray]
	dlogTdlogP = dlogTdlogP[boolArray]
	Hfrac = Hfrac[boolArray]
	R = R[boolArray]

	R = 1e1**logR
	pipe = subprocess.Popen("./hydro", stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
	#pipe = subprocess.Popen("./hydro", stdin=subprocess.PIPE)

	input_string = []
	input_string.append( "{0:12.10E} {1:12.10E}\n".format(massCentral, 1e1**logRCentral))
	input_string.append("{0}\n".format(R.size))
	for m, lgr, lgrho, lgT, lgP, dlTdlP, Xfrac in list(zip( mass, logR, logRho, logT, logP, dlogTdlogP, Hfrac))[::-1] :
		input_string.append( "{0:12.10E} {1:12.10E} {2:12.10E} {3:12.10E} {4:12.10E} {5:12.10E} {6:12.10E}\n".format( m, lgr, lgrho, lgT, lgP, dlTdlP, Xfrac))
	output_string = pipe.communicate( "".join(input_string))[0]
	#pipe.communicate( "".join(input_string))[0]
	#print output_string

	#mcore, rcore = np.genfromtxt(io.BytesIO(output_string),max_rows=1,unpack=True)
	mcore, rcore = np.genfromtxt(io.StringIO(output_string),max_rows=1,unpack=True)

	print( "core mass: {0}, core radius = {1}".format(mcore,rcore))
	m,r,rho,T,p,Xout = np.genfromtxt(io.StringIO(output_string),skip_header=2,unpack=True)
	#m,r,rho,T,p,Xout = np.genfromtxt(io.BytesIO(output_string),skip_header=2,unpack=True)	
	print( "maxR: {0}, Rmax = {1}".format(maxR,r.max()))
	boolArray = filterMonotonicity(rho)
   # for dens,temp,g in zip( rho, T,gamma):
   #   print dens, temp, g
	if( makePlot ) :
		plt.semilogy(r, rho, ls="dashed", lw=3, color="green")
		plt.savefig( "test.pdf")

	return mcore, rcore, m[boolArray], r[boolArray], rho[boolArray], T[boolArray], Xout[boolArray]

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('profilefile', default="profile15.dat", help="<profilename>")
parser.add_argument('--no_atm', action='store_true', help='no atmosphere surrounding object')
parser.add_argument('--sph_ic', action='store_true', help='sph initial conditions (equal mass particles)')
parser.add_argument('--add_companion', nargs=1, type=float, help='add dark matter companion')
parser.add_argument('--rotation', nargs=1, type=float, help='rotation of primary as fraction of solar rotation')
parser.add_argument('--corotation', nargs=1, type=float, help='rotation of primary as fraction of companions orbit')
parser.add_argument('--comp_radius', nargs=1, type=float, help='companion-core separation in solar radii')
parser.add_argument('--bh', action='store_true', help='companion has a small softening length')
args = parser.parse_args()

rot = None
if( args.rotation != None) :
	rot = args.rotation[0]
	omega = rot * rotSun

corot = None
if( args.corotation != None) :
	corot = args.corotation[0]

mComp = None
if( args.add_companion != None) :
	mComp = args.add_companion[0]
print( 'mComp = ' + str(mComp))
N = 32
posArray, velArray = makePosVel( N)

# now center
massParticleInGram = mSun/16e6

xmax = 2500.0*7.0e10
rhoExt = 1.0e-14
tempExt = 1.0e5

ymax = xmax
zmax = xmax

p1, v1, e1, X1, m1, rho_max, mcore, rcore, mStar = makeStarMesh( args.profilefile, massParticleInGram, posArray, args.sph_ic)
mStar = mStar + mcore

pos  = p1
vel  = v1
ener = e1
mass = m1
Hfrac = X1

numberStarParticles = e1.size

r1Max = np.linalg.norm(p1[-1])
print( "Number of star particles: {0}, Rmax = {1}".format(numberStarParticles,r1Max))
r1 = np.array([0.,0.,0.])

# check xmax
if( xmax < r1Max) :
	print( 'warning')

if( args.comp_radius != None) :
        rComp = args.comp_radius[0] * rSun
else :
        rComp = r1Max

if( mComp != None and mComp > 0) :
	mComp = mComp*mSun
	xComp = rComp
	yComp = 0.
	zComp = 0.
	mtot = mStar + mComp
	vrel  = math.sqrt( mtot/rComp)
	v1 = -mComp/mtot * vrel
	v2 = mStar/mtot*vrel
	print( "Mstar = {0}, mComp = {1}, v1 = {2}, v2 = {3}, rComp = {4}".format(mStar,mComp, v1, v2, rComp))

omega = 0.0
r = np.linalg.norm(pos, axis=1 )

if( corot != None) :
	rotComp = vrel/rComp
	omega = corot * rotComp

vel[:,0] = vel[:,0] - omega * pos[:,1]
vel[:,1] = vel[:,1] + omega * pos[:,0]

# do the external gas
if( not args.no_atm) :
	# find the mean separation
	dr1 = np.linalg.norm( p1[0:-1,:] - p1[-1,:], axis=1)
	dr1.sort()
	h = dr1[0:16].mean()
	'''
	rad = np.linalg.norm( p1, axis=1 )
	shellBoolArray = rad > ( r1Max * 0.95 )
	numShell = shellBoolArray.sum()
	print 'numShell = ', numShell
	p1Shell = p1[shellBoolArray,:]
	h1 = np.ones(numShell)
	for i in range(0,numShell):
		dr1 = np.linalg.norm( p1 - p1Shell[i,:], axis=1 )
		dr1.sort()
		h1[i] = dr1[1:17].mean()
	h = h1.max()
	'''
	print( 'h = ', h)

	massParticleExtInGram = rhoExt*h*h*h
	massParticleExtInCode = massParticleExtInGram

	# find the total mass
	totalExtMass = rhoExt * 8. * xmax*xmax*xmax
	numberExtParticles = totalExtMass/massParticleExtInCode
	print( "posArray: ", posArray.size/3, " numberExtParticles: ", numberExtParticles)
	#find the scaling
	scaling = (posArray.size/3/numberExtParticles)**0.3333/(0.5*N)

	ExtPosArray = posArray*scaling*xmax

	radius = np.linalg.norm(ExtPosArray, axis=1 )
	rescaleFactor = ScaleFactor( radius, r1Max)
	ExtPosArray[:,0] = ExtPosArray[:,0] * rescaleFactor
	ExtPosArray[:,1] = ExtPosArray[:,1] * rescaleFactor
	ExtPosArray[:,2] = ExtPosArray[:,2] * rescaleFactor

	# cut the cube first
	boolArray = ExtPosArray[:,0] < xmax
	boolArray = np.logical_and( boolArray, ExtPosArray[:,0] > -xmax)

	boolArray = np.logical_and( boolArray, ExtPosArray[:,1] < ymax)
	boolArray = np.logical_and( boolArray, ExtPosArray[:,1] > -ymax)

	boolArray = np.logical_and( boolArray, ExtPosArray[:,2] < zmax)
	boolArray = np.logical_and( boolArray, ExtPosArray[:,2] > -zmax)

	# now exclude the regions around star
	boolArray = np.logical_and( boolArray, np.linalg.norm(ExtPosArray - r1, axis=1) > r1Max)

	ExtPosArray = ExtPosArray[boolArray]
	rescaleFactor = rescaleFactor[boolArray]

	radius = np.linalg.norm(ExtPosArray, axis=1 )

	# create the other arrays
	NumExtParticles = int(ExtPosArray.size/3)
	ExtVelArray = np.zeros([NumExtParticles,3])
	ExtMArray = massParticleExtInCode*np.ones(NumExtParticles)
	ExtEArray = tempExt*np.ones(NumExtParticles)
	ExtXArray = 0.7*np.ones( NumExtParticles)

	ExtVelArray[:,0] = - r1Max*r1Max/radius/radius*ExtPosArray[:,1]*omega
	ExtVelArray[:,1] = r1Max*r1Max/radius/radius*ExtPosArray[:,0]*omega

	print( 'Number of external particles: ' + str(NumExtParticles)  + " " + str(numberExtParticles))

	ExtMArray = ExtMArray * rescaleFactor**3 # ScaleFactor( radius, r1Max )**3

	totalExtMassPcles = ExtMArray.sum()
	rhoExtActual = totalExtMassPcles / (8. * xmax*xmax*xmax - 4./3.*3.1416*r1Max*r1Max*r1Max)
	#print "rhoExtActual, ", rhoExtActual
	ExtMArray = ExtMArray * rhoExt / rhoExtActual

	# now append
	pos = np.append(pos, ExtPosArray, axis=0)
	vel = np.append(vel, ExtVelArray, axis=0)
	ener = np.append(ener, ExtEArray)
	mass = np.append(mass, ExtMArray)
	Hfrac = np.append( Hfrac, ExtXArray)

x =  pos[:,0]
y =  pos[:,1]
z =  pos[:,2]
vx = vel[:,0]
vy = vel[:,1]
vz = vel[:,2]

print( "Number of total particles: {0}".format(mass.size))

xboolArray = np.logical_and( x <= xmax, x > -xmax)
yboolArray = np.logical_and( y <= ymax, y > -ymax)
zboolArray = np.logical_and( z <= zmax, z > -zmax)

boolArray = np.logical_and( np.logical_and(xboolArray, yboolArray), zboolArray)
x = x[boolArray]
y = y[boolArray]
z = z[boolArray]
vx = vx[boolArray]
vy = vy[boolArray]
vz = vz[boolArray]

mass = mass[boolArray]
#rhos = rhoArray[boolArray]
temp = ener[boolArray]
metals = Hfrac[boolArray]
#pArray = pArray[boolArray]
print( "rho_max = ", rho_max)
print( "mcore = ", mcore)
xcore = np.zeros(1)
ycore = np.zeros(1)
zcore = np.zeros(1)
vxcore = np.zeros(1)
vycore = np.zeros(1)
vzcore = np.zeros(1)
mdm = np.array([])
rdm = np.array([])
if( mcore > 0.) :
	mdm = np.array([mcore])
	rdm = np.array([rcore])

if( mComp != None and mComp > 0) :
	xcore = np.array([0., xComp])
	ycore = np.array([0., yComp])
	zcore = np.array([0., zComp])
	vxcore = np.zeros(2)
	vzcore = np.zeros(2)
	vycore = np.array([v1, v2])
	vy = vy + v1
	mdm = np.append( mdm, mComp)
	if args.bh :
		rdm = np.append( rdm, 2.0*h )
	else :
		rdm = np.append( rdm, mComp/mSun*rSun )

if( mcore > 0. or mComp is not None) :
	write_tipsy( "star", x, y, z, vx, vy, vz, mass, temp, rho_max, gasmetals=metals,
		xdm=xcore, ydm=ycore, zdm=zcore, vxdm=vxcore, vydm=vycore, vzdm=vzcore, mdm=mdm, dmeps=rdm)
else :
	write_tipsy( "star", x, y, z, vx, vy, vz, mass, temp, rho_max, gasmetals=metals)
