import math
import numpy as np
from scipy.interpolate import interp1d

from readWriteTipsy import *

maxRatio = 100
rscaling = 3
N = 24

def ScaleFactor( radius, rStar) :
  global maxRatio, rscaling
  scale = (radius / rStar)**rscaling
  scale = np.maximum( scale, 1.)
  scale = np.minimum( scale, maxRatio)
  return scale

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

def makeStarMesh(mcore, rcore, mModel, rModel, rhoModel, eModel, XModel, vModel, massParticleInGram, posArray, rMax=None, sphICs=False, maxGrid=None) :
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
  vrArray = np.zeros(rArray.size)
  vArray = np.zeros([rArray.size, 3])
  XArray = np.zeros(rArray.size)
  mArray = np.ones(rArray.size)
  lastModelP = 0

  if( vModel is None) : 
    vModel = np.zeros(rModel.size)

  if( not sphICs) :
    # now rescale
    rescale = True
    scale = 1
    if( rMax is None) :
      rMax = rModel.max()
      rescale=False
    rArray = rArray*rMax/rCircle*(rArray.size/numParticles)**0.33333
    if( rescale) :
      scale = ScaleFactor( rArray, rMax)
      rArray *= scale
      scale3 = scale**3

    interpolate = interp1d(rModel,rhoModel,kind = 'linear')
    interpolateE = interp1d(rModel,eModel,kind = 'linear')
    interpolateX = interp1d(rModel,XModel,kind = 'linear')
    interpolateV = interp1d(rModel,vModel,kind = 'linear')

    #for i in range(numParticles) :
      #if( i%10000 == 0) : 
      #  print( "Working on particle number {0} r = {1:.2e} scale={2:.1e}".format(i, rArray[i], scale[i]))
      #if( (rArray[i] > rMax and not rescale) or rArray[i] >= rModel.max() or (not maxGrid is None and rArray[i] > maxGrid)) :
      #  print("breaking at ", rArray[i])
      #  break
          # vol = 4.*math.pi/3.*rMax**3/numParticles
          # mArray[i] = interpolate( rArray[i])*vol
    rlimit = rModel.max()
    if( not rescale) :
      rlimit = rMax
    if( not maxGrid is None) : 
      rlimit = min(rlimit,maxGrid)
    boolArray = rArray < rlimit
    
    mArray[boolArray] = interpolate( rArray[boolArray])
    eArray[boolArray] = interpolateE( rArray[boolArray])
    XArray[boolArray] = interpolateX( rArray[boolArray])
    posArray[boolArray,:] = posArray[boolArray,:]* (rArray[boolArray]/np.linalg.norm(posArray[boolArray,:],axis=-1))[:,np.newaxis]
    vrArray[boolArray] = interpolateV( rArray[boolArray])
    vArray[boolArray,:] = posArray[boolArray,:]*(vrArray[boolArray]/np.linalg.norm(posArray[boolArray,:],axis=-1))[:,np.newaxis]
    lastModelP = rArray[boolArray].size
    # adjust for scaling 
    mArrayOut = mArray[boolArray]*scale3[boolArray] * mass / (mArray[boolArray]*scale3[:lastModelP]).sum()*0.5
    print("rlimit = {0:.2e} {1} {2:.2e}".format(rlimit, rArray[boolArray][-1], np.linalg.norm(posArray[boolArray][-1])))
  else : 
    mArray = np.arange(1,rArray.size)*massParticleInGram
    interpolate = interp1d(mModel,rModel,kind = 'linear')
    interpolateE = interp1d(mModel,eModel,kind = 'linear')
    interpolateX = interp1d(mModel,XModel,kind = 'linear')
    interpolateV = interp1d(mModel,vModel,kind = 'linear')
    lastModelP = 0
    for i in range(numParticles) :
      if( mArray[i] > mModel[-1] ) :
        break
      lastModelP = i
      r = interpolate( mArray[i])
      eArray[i] = interpolateE( mArray[i])
      XArray[i] = interpolateX( mArray[i])
      posArray[i,:] = posArray[i,:]/rArray[i] * r
      vrArray[i] = interpolateV( rArray[i])
    lastModelP = lastModelP+1
    mArrayOut = massParticleInCode*np.ones(lastModelP)
  return posArray[:lastModelP], vArray[:lastModelP], eArray[:lastModelP], XArray[:lastModelP], mArrayOut, rho_max, mcore, rcore, mass



def makeStar(mcore, rcore, mModel, rModel, rhoModel, eModel, XModel, vModel, xmax, ymax, zmax, massParticleInGram, rhoExt, tempExt, rMax=None, outfile="star") :
    global N
    posArray, velArray = makePosVel( N)
    pos, vel, ener, Hfrac, mass, rho_max, mcore, rcore, mStar = makeStarMesh( mcore, rcore, mModel, rModel, rhoModel, eModel, XModel, vModel, massParticleInGram, posArray, rMax=rMax, sphICs=False, maxGrid=1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2))
    mStar = mStar + mcore

    numberStarParticles = ener.size
    print(pos.shape)
    rStar = np.linalg.norm(pos[-1,:])
    print( "Number of star particles: {0}, Rmax = {1:.2e}".format(numberStarParticles,rStar))
    r1 = np.array([0.,0.,0.])
    # check xmax
    if( xmax < rStar) :
        print( 'warning')


    # do the external gas
    if( not rhoExt is None) :
        # find the mean separation
        dr1 = np.linalg.norm( pos[0:-1,:] - pos[-1,:], axis=1)
        dr1.sort()
        h = dr1[0:16].mean()
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
        rescaleFactor = ScaleFactor( radius, rStar)
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
        boolArray = np.logical_and( boolArray, np.linalg.norm(ExtPosArray - r1, axis=1) > rStar)

        ExtPosArray = ExtPosArray[boolArray]
        rescaleFactor = rescaleFactor[boolArray]

        radius = np.linalg.norm(ExtPosArray, axis=1 )

        # create the other arrays
        NumExtParticles = int(ExtPosArray.size/3)
        ExtVelArray = np.zeros([NumExtParticles,3])
        ExtMArray = massParticleExtInCode*np.ones(NumExtParticles)
        ExtEArray = tempExt*np.ones(NumExtParticles)
        ExtXArray = 0.7*np.ones( NumExtParticles)

        #ExtVelArray[:,0] = - rStar*rStar/radius/radius*ExtPosArray[:,1]*omega
        #ExtVelArray[:,1] = rStar*rStar/radius/radius*ExtPosArray[:,0]*omega

        print( 'Number of external particles: ' + str(NumExtParticles)  + " " + str(numberExtParticles))

        ExtMArray = ExtMArray * rescaleFactor**3 # ScaleFactor( radius, rStar )**3

        totalExtMassPcles = ExtMArray.sum()
        rhoExtActual = totalExtMassPcles / (8. * xmax*xmax*xmax - 4./3.*3.1416*rStar*rStar*rStar)
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
    print( "rho_max = {0:.3e} ".format(rho_max))
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

    if( mcore > 0.) :
	    write_tipsy( outfile, x, y, z, vx, vy, vz, mass, temp, rho_max, gasmetals=metals,
		    xdm=xcore, ydm=ycore, zdm=zcore, vxdm=vxcore, vydm=vycore, vzdm=vzcore, mdm=mdm, dmeps=rdm)
    else :
	    write_tipsy( outfile, x, y, z, vx, vy, vz, mass, temp, rho_max, gasmetals=metals)
