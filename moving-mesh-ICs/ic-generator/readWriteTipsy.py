import numpy as np
import struct
import os

def swap_endian( i) :
    return struct.unpack("<I", struct.pack(">I",i))[0]
def swap_endianf( i) :
	return struct.unpack("<f", struct.pack(">f",i))[0]

def readSquare():
    N = 16
    pos = np.zeros([N*N*N,3])
    for i in range(N) :
        for j in range(N) :
            for k in range(N) :
                pos[i*N*N + j*N + k,0] = (1.0*i+0.5)/N
                pos[i*N*N + j*N + k,1] = (1.0*j+0.5)/N
                pos[i*N*N + j*N + k,2] = (1.0*k+0.5)/N
    vel = np.zeros([N*N*N,3])
    return pos, vel

def readPosVel():
    f = open("d2300.glass", "r")

    nf77,v0,v1,v2,nf77e = np.fromfile(f,dtype=np.int32,count=5)
    array = np.fromfile(f,dtype=np.int32,count=402)
    nobj=swap_endian(array[302]) # get the number of gas particles

    # read type junk
    nf77 = np.fromfile(f,dtype=np.int32,count=1)
    types = np.fromfile(f,dtype=np.float32,count=nobj)
    nf77e = np.fromfile(f,dtype=np.int32,count=1)

    # read mass junk
    nf77 = np.fromfile(f,dtype=np.int32,count=1)
    mass = np.fromfile(f,dtype=np.float32,count=nobj)
    nf77e = np.fromfile(f,dtype=np.int32,count=1)

    # read pos junk
    nf77 = np.fromfile(f,dtype=np.int32,count=1)
    pos = np.fromfile(f,dtype=np.float32,count=nobj*3)
    nf77e = np.fromfile(f,dtype=np.int32,count=1)

    # read velocity
    nf77 = np.fromfile(f,dtype=np.int32,count=1)
    vel = np.fromfile(f,dtype=np.float32,count=nobj*3)
    nf77e = np.fromfile(f,dtype=np.int32,count=1)

    for i in range(vel.size) :
        vel[i] = swap_endianf(vel[i])
        pos[i] = swap_endianf(pos[i])

    return pos.reshape([nobj,3]), vel.reshape([nobj,3])

def makePosVelOld() :
    pos = np.zeros( [4096,3])
    N = 16
    ipart = 0
    for i in range(N) :
    	for j in range(N) :
            for k in range(N) :
                noise = (np.random.rand(3) - 0.5) * 0.025/N
                x = (1.*i + 0.5)/N + noise[0]
                y = (1.*j + 0.5)/N + noise[1]
                z = (1.*k + 0.5)/N + noise[2]
                pos[ipart,0] = x
                pos[ipart,1] = y
                pos[ipart,2] = z
                ipart = ipart + 1
    vel = pos*0.

    return pos, vel

def makePosVel( NX, NY=None, NZ=None, useGlass=True) :
    if( NY == None) :
        NY=NX
    if( NZ == None) :
        NZ=NX
    pos = None
    vel = None
    if( useGlass) :
        pos, vel = readPosVel()
    else :
        pos, vel = makePosVelOld()

    xc = np.arange(0., NX, 1.)
    yc = np.arange(0., NY, 1.)
    zc = np.arange(0., NZ, 1.)
    numParticles = int(pos.size/3)
    posArray = np.zeros([int(NX*NY*NZ*numParticles),3])
    velArray = np.zeros([int(NX*NY*NZ*numParticles),3])

    totNumParticles = posArray.size/3
    offset = np.zeros(3)

    for i in range(NX) :
        offset[0] = 1.*i - 0.5*NX;

        for j in range(NY) :
            offset[1] = 1.*j - 0.5*NY;

            for k in range(NZ) :
                offset[2] = 1.*k - 0.5*NZ;
                pStart = int(i*NY*NZ+j*NZ + k)*numParticles
                pEnd = pStart+numParticles
                posArray[pStart:pEnd,:] = pos[0:numParticles,:] + offset[:]
                velArray[pStart:pEnd,:] = vel[0:numParticles,:]*0.

    return posArray, velArray

def write_tipsy( filename, xs, ys, zs, vx, vy, vz, ms, temp, rho_max,
    xdm=None, ydm=None, zdm=None, vxdm=None, vydm=None, vzdm=None, mdm=None, dmeps=None,
    gasmetals=None, STANDARD=True):#, gsmooths) :
    f = None
    try:
        f = open(filename+".std", 'wb')
    except:
        print("WTIPSY ERROR: Can't open file")
        return 1

    print("Writing {0}.std".format(filename))
    time = 0.
    ngas = xs.size
    print ('Number of gas particles for simulation box is {0}\n'.format( ngas))
    nbodies = ngas # isn't this just the number of stars?
    nstars = 0
    ndims = 3
    ntime = 0 # not sure what this represents

    darkphi = np.array([])
    if( xdm is None) :
        xdm = np.array([])
        ydm = np.array([])
        zdm = np.array([])
        vxdm = np.array([])
        vydm = np.array([])
        vzdm = np.array([])
        mdm = np.array([])
        dmeps = np.array([])
    else :
        darkphi = np.zeros(xdm.size)
    ndm = xdm.size
    nbodies = nbodies + ndm
    gasmass = ms
    darkmass = mdm
    starmass = np.array([])

    gasx = xs
    darkx = xdm
    starx = np.array([])
    gasy = ys
    darky = ydm
    stary = np.array([])
    gasz = zs
    darkz = zdm
    starz = np.array([])

    gasvx = vx
    darkvx = vxdm
    starvx = np.array([])
    gasvy = vy
    darkvy = vydm
    starvy = np.array([])
    gasvz = vz
    darkvz = vzdm
    starvz = np.array([])

    #Other crap.
    darkeps = dmeps
    stareps = np.array([])
    gasdens = np.ones(ngas) # this is calculated later
    if(gasmetals is None) :
        gasmetals = np.zeros(ngas)
    starmetals = np.array([])
#    temp = 8*nbodies + 3*ngas + nstars
    startform = np.array([])
#    temp = 8*nbodies + 3*ngas + 2*nstars
    gasphi = np.zeros(ngas)
    starphi = np.array([])
    gastemp = temp

    #Get OLD smoothing length
    gashsmooth = np.ones([ngas])
    if( ngas > 0) : 
       gasr = np.sqrt(gasx**2 + gasy**2 + gasz**2)
       gashsmoothOLD = 3.0*np.ones(ngas)*max(gasr)/(nbodies**(1.0/3.0))
       #Get correct smoothing length.
       gashsmooth = 1.0*(gasmass[0]/rho_max)**(1.0/3.0)*np.ones(ngas)    #min interparticle separation shoul    d be reasonable
      
    newstring = repr(nbodies) + ' ' + repr(ngas) + ' ' + repr(nstars) + '\n'

    f.write(struct.pack(">diiiii", time, nbodies, ndims, ngas, ndm, nstars))
    if STANDARD:
        f.write(struct.pack("xxxx"))
    for i in range(ngas):
        f.write(struct.pack(">ffffffffffff", gasmass[i], gasx[i], gasy[i], gasz[i], gasvx[i], gasvy[i], 
            gasvz[i], gasdens[i], gastemp[i], gashsmooth[i], gasmetals[i], gasphi[i]))
    for i in range(ndm):
        f.write(struct.pack(">fffffffff", darkmass[i], darkx[i], darky[i], darkz[i], darkvx[i], darkvy[i], 
            darkvz[i], darkeps[i], darkphi[i]))
    for i in range(nstars):
        f.write(struct.pack(">fffffffffff", starmass[i], starx[i], stary[i], starz[i], starvx[i], starvy[i], 
            starvz[i], starmetals[i], startform[i], stareps[i], starphi[i]))
    f.close()


def read_tipsy(filename, VERBOSE=False):
    """rtipsy Reads tipsy files detecting the format: 
    big endian, little endian, padded (standard) or non-padded header 
    Usage: 
            rtipsy(filename, VERBOSE=False)
    Input parameters: 
    filename  filename string
    VERBOSE  print messages (optional)
    Return values:
    (header,g,d,s)
    header    tipsy header struct
    g,d,s     gas, dark and star structures
    Please read rtipsy.py for the structure definitions
    Example: 
    h,g,d,s = rtipsy('/home/wadsley/usr5/mihos/mihos.std')
    print, h['ndark']
    plt.plot(d['x'], d['y'], 'k,')"""
    try:
        f = open(filename, 'rb')
    except:
        print( "RTIPSY ERROR: Can't open file")
        return 1
    fs = len(f.read())
    f.seek(0)
    #Read in the header
    t, n, ndim, ng, nd, ns = struct.unpack("<diiiii", f.read(28))
    endianswap = False
    #Check Endianness
    if (ndim < 1 or ndim > 3):
        endianswap = True
        f.seek(0)
        t, n, ndim, ng, nd, ns = struct.unpack(">diiiii", f.read(28))
        if VERBOSE:
            print( "SWAP_ENDIAN")
    if VERBOSE:
        print( "Read time,n,ngas,ndark,nstar: ", t, n, ng, nd, ns)
    #Catch for 4 byte padding
    if (fs == 32+48*ng+36*nd+44*ns):
        f.read(4)
    #File is borked if this is true
    elif (fs != 28+48*ng+36*nd+44*ns):
        print( "RTIPSY ERROR: Header and file size inconsistent")
        print( "Estimates: Header bytes:  28 or 32 (either is OK)")
        print( "     ngas: ",ng," bytes:",48*ng)
        print( "    ndark: ",nd," bytes:",36*nd)
        print( "    nstar: ",ns," bytes:",44*ns)
        print( "Actual File bytes:",fs,"  not one of:",28+48*ng+36*nd+44*ns,32+48*ng+36*nd+44*ns)
        f.close()
        return 1

    catg = {'mass':np.zeros(ng), 'pos':np.zeros((ng,3)), 'vel':np.zeros((ng,3)), 'dens':np.zeros(ng),
            'tempg':np.zeros(ng), 'h':np.zeros(ng), 'zmetal':np.zeros(ng),	'phi':np.zeros(ng)}
    catd = {'mass':np.zeros(nd), 'pos':np.zeros((nd,3)), 'vel':np.zeros((nd,3)),
            'eps':np.zeros(nd), 'phi':np.zeros(nd)}
    cats = {'mass':np.zeros(ns), 'pos':np.zeros((ns,3)), 'vel':np.zeros((ns,3)),
            'metals':np.zeros(ns), 'tform':np.zeros(ns), 'eps':np.zeros(ns), 'phi':np.zeros(ns)}
    for cat in ['g','d','s']:
        j = 0
        for qty in ['x','y','z']:
            locals()['cat'+cat][qty] = locals()['cat'+cat]['pos'][:,j]
            locals()['cat'+cat]['v'+qty] = locals()['cat'+cat]['vel'][:,j]
            j += 1

    if (ng > 0):
        for i in range(ng):
            if endianswap:
                mass, x, y, z, vx, vy, vz, dens, tempg, h, zmetal, phi = struct.unpack(">ffffffffffff", f.read(48))
            else:
                mass, x, y, z, vx, vy, vz, dens, tempg, h, zmetal, phi = struct.unpack("<ffffffffffff", f.read(48))
            catg['mass'][i] = mass
            catg['x'][i] = x
            catg['y'][i] = y
            catg['z'][i] = z
            catg['vx'][i] = vx
            catg['vy'][i] = vy
            catg['vz'][i] = vz
            catg['dens'][i] = dens
            catg['tempg'][i] = tempg
            catg['h'][i] = h
            catg['zmetal'][i] = zmetal
            catg['phi'][i] = phi
    if (nd > 0):
        for i in range(nd):
            if endianswap:
                mass, x, y, z, vx, vy, vz, eps, phi = struct.unpack(">fffffffff", f.read(36))
            else:
                mass, x, y, z, vx, vy, vz, eps, phi = struct.unpack("<fffffffff", f.read(36))
            catd['mass'][i] = mass
            catd['x'][i] = x
            catd['y'][i] = y
            catd['z'][i] = z
            catd['vx'][i] = vx
            catd['vy'][i] = vy
            catd['vz'][i] = vz
            catd['eps'][i] = eps
            catd['phi'][i] = phi
    if (ns > 0):
        for i in range(ns):
            if endianswap:
                mass, x, y, z, vx, vy, vz, metals, tform, eps, phi = struct.unpack(">fffffffffff", f.read(44))
            else:
                mass, x, y, z, vx, vy, vz, metals, tform, eps, phi = struct.unpack("<fffffffffff", f.read(44))
            cats['mass'][i] = mass
            cats['x'][i] = x
            cats['y'][i] = y
            cats['z'][i] = z
            cats['vx'][i] = vx
            cats['vy'][i] = vy
            cats['vz'][i] = vz
            cats['metals'][i] = metals
            cats['tform'][i] = tform
            cats['eps'][i] = eps
            cats['phi'][i] = phi
    header = {'time':t, 'n':n, 'ndim':ndim, 'ngas':ng, 'ndark':nd, 'nstar':ns}
    return (header,catg,catd,cats)

