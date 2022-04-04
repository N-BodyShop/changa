!
! Initialized the helmholtz equation of state
!
      subroutine initializeHelmholtz
      
      include 'implno.dek'
      include 'vector_eos.dek'

      call read_helm_table

      return
      end

!
! get the values of the equation of state
! 
      subroutine getHelmholtzEOS( temp, dens, abar, zbar, pressure, energy, soundspeed )

      include 'implno.dek'
      include 'vector_eos.dek'

      double precision temp, dens, abar, zbar, pressure, energy, soundspeed
 
      temp_row(1) = temp ; den_row(1)  = dens ; abar_row(1) = abar ; zbar_row(1) = zbar
      jlo_eos = 1 ; jhi_eos = 1

      call helmeos
      pressure = ptot_row(1)
      energy = etot_row(1)
      soundspeed = cs_row(1)

!pion = pion_row(1)

!temp = pion/(avogadro's number*kb/abar)

      return
      end


!      program teos
!      implicit none
! tests the eos routine
! 
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

!      integer          ionmax
!      parameter        (ionmax=3)
!      double precision xmass(ionmax),aion(ionmax),zion(ionmax),temp,den,abar,zbar
!      external         initializeHelmholtz, getHelmholtzEOS
!      double precision pressure

! set the mass fractions, z's and a's of the composition
! hydrogen, heliu, and carbon
!      xmass(1) = 0.75d0 ; aion(1)  = 1.0d0  ; zion(1)  = 1.0d0
!      xmass(2) = 0.23d0 ; aion(2)  = 4.0d0  ; zion(2)  = 2.0d0
!      xmass(3) = 0.02d0 ; aion(3)  = 12.0d0 ; zion(3)  = 6.0d0

! average atomic weight and charge
!      abar   = 1.0d0/sum(xmass(1:ionmax)/aion(1:ionmax))
!      zbar   = abar * sum(xmass(1:ionmax) * zion(1:ionmax)/aion(1:ionmax))

! set the input vector. pipeline is only 1 element long
!      call initializeHelmholtz
!      call getHelmholtzEOS( 1D4, 1D6, abar, zbar, pressure) 
      
!      write (*,*) pressure

!      end   


      subroutine getEOSInversion( temp, dens, abar, zbar, pressure, energy, soundspeed, success)

      include 'implno.dek'
      include 'vector_eos.dek'

      double precision temp, dens, abar, zbar, pressure, energy, soundspeed
      logical success
 
      temp_row(1) = temp ; den_row(1)  = dens ; abar_row(1) = abar ; zbar_row(1) = zbar
      etot_row(1) = energy
      jlo_eos = 1 ; jhi_eos = 1

      call invert_helm_ed

      pressure = ptot_row(1)
      temp = temp_row(1)
      soundspeed = cs_row(1)
      success = .not. eosfail

      return
      end

      subroutine invert_helm_ed
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the specific internal energy density, density, and composition
! find everything else

! it is assumed that etot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision tmpold,tmp,f,df,tmpnew,eostol,eostola,fpmin
      double precision eostolfail
      parameter        (eostola = 1.0d-8, &
                        fpmin  = 1.0d-14)

!if error is greater than eostolfail, stop the simulation;
!otherwise just report it.
      if(den_row(1) .lt. 1.0d10) then
	eostol = 1.0d-8
	eostolfail = 1.0d-6
      elseif(den_row(1) .ge. 1.0d10 .and. &
		den_row(1) .lt. 1.0d12) then
	eostol = 1.0d-7
	eostolfail = 1.0d-6
      elseif(den_row(1) .ge. 1.0d12) then
	eostol = 1.0d-6
	eostolfail = 1.0d-5
      endif


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = etot_row(j)	!original energy
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe

      call helmeos

      do j = jlo_eos, jhi_eos

       f     = etot_row(j)/eoswrk03(j) - 1.0d0
       df    = det_row(j)/eoswrk03(j)

!Newton-Raphson: T1 = T0 - (E(T0)-E_exp)/E'(T0); we want to get
!the temperature for E_exp, which is why we subtract E_exp.
!Timmes decides to do this via ratios, hence f and df
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       tmp    = temp_row(j)
       tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
       eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new temperature, keep it within the table limits
       temp_row(j)  = min(1.0d14,max(tmpnew,1.0d-11))


      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos

        f     = etot_row(j)/eoswrk03(j) - 1.0d0
        df    = det_row(j)/eoswrk03(j)

        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,1.0d3))


! end of netwon loop
       end do


! we did not converge if we land here
      write(0,*)
      write(0,*) 'newton-raphson failed in routine invert_helm_ed'
      write(0,*) 'pipeline element',j
      write(0,*) 'ewant  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(0,*) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(0,*) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j)
      write(0,*) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(0,*) 'density = ',den_row(j)
      write(0,*) 'ehave =',etot_row(j)
      write(0,*)
      if(eoswrk01(j) .gt. eostolfail) then
	stop 'could not find a temperature in routine invert_helm_ed'
      endif


! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos

      return
      end




      subroutine getPDInversion( temp, dens, abar, zbar, pressure, energy)

      include 'implno.dek'
      include 'vector_eos.dek'

      double precision temp, dens, abar, zbar, pressure, energy, soundspeed
 
      den_row(1) = dens ; ptot_row(1) = pressure ; abar_row(1) = abar ; zbar_row(1) = zbar
      jlo_eos = 1 ; jhi_eos = 1

      temp_row(1) = temp

      call invert_helm_pd

      temp = temp_row(1)
      energy = etot_row(1)
      soundspeed = cs_row(1)

      return
      end

      subroutine invert_helm_pd
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


! given the pressure, density, and composition
! find everything else

! it is assumed that ptot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! this version is quiet on all errors


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision temp,f,df,tempnew,eostol,fpmin
      parameter        (eostol = 1.0d-4, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpt_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       temp    = temp_row(j)
       tempnew = min(max(0.5d0*temp,temp - eoswrk02(j)),2.0d0*temp)

! compute the error
       eoswrk01(j)  = abs((tempnew - temp)/temp)

! store the new density, keep it within the table limits
       temp_row(j)  = min(1.0d12,max(tempnew,1.0d3))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpt_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
       temp    = temp_row(j)
       tempnew = min(max(0.5d0*temp,temp - eoswrk02(j)),2.0d0*temp)

! compute the error
       eoswrk01(j)  = abs((tempnew - temp)/temp)

! store the new density, keep it within the table limits
       temp_row(j)  = min(1.0d12,max(tempnew,1.0d3))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_pd'
      write(6,*) 'pipeline element',j
      write(6,01) 'pres  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                 '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a density in routine invert_helm_pd'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the temperature

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos

      return
      end









! here is the tabular helmholtz free energy eos:
!
! routine read_helm_table reads an electron helm free energy table
! routine read_helm_iontable reads an ion free energy table
! routine helmeos computes the pressure, energy and entropy via tables
! routine helmeos3 uses fermionic ions
! routine helmeos2 adds third derivatives
! routine helmeos_orig as helmeos with the table read inside the routine



      subroutine read_helm_table
      include 'implno.dek'
      include 'helm_table_storage.dek'

! this routine reads the helmholtz eos file, and
! must be called once before the helmeos routine is invoked.

! declare local variables
      integer          i,j
      double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i, &
                       dd,dd2,ddi,dd2i,dd3i


! open the file (use softlinks to input the desired table)

       open(unit=19,file='helm_table.dat',status='old')

! for standard table limits
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -12.0d0
       dhi   = 15.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp

! read the helmholtz free energy and its derivatives
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                  fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo


! read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

! read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

! read the number density table
       do j=1,jmax
        do i=1,imax
         read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

! close the file and write a summary message
      close(unit=19)


! construct the temperature and density deltas and their inverses
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt3i        = dt2i*dti
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
        dt3i_sav(j) = dt3i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd3i        = dd2i*ddi
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
        dd3i_sav(i) = dd3i
       enddo



!      write(6,*)
!      write(6,*) 'finished reading eos table'
!      write(6,04) 'imax=',imax,' jmax=',jmax
!04    format(1x,4(a,i4))
!      write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
!      write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
!03    format(1x,4(a,1pe11.3))
!      write(6,*)

      return
      end





      subroutine helmeos
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      include 'helm_table_storage.dek'


! given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar and zbar, this routine returns most of the other
! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
! their derivatives with respect to temperature, density, abar, and zbar.
! other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.
!
! this routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. interpolation in a table of the helmholtz free energy
! is used to return the electron-positron thermodynamic quantities.
! all other derivatives are analytic.
!
! references: cox & giuli chapter 24 ; timmes & swesty apj 1999


! declare
      integer          i,j
      double precision temp,den,abar,zbar,ytot1,ye, &
                       x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                       dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                       dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                       deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                       dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                       sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                       dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                       gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                       detadt,detadd,xnefer,dxnedt,dxnedd,s

      double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz, &
                       egas,degasdd,degasdt,degasda,degasdz, &
                       sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz, &
                       cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas, &
                       chit_gas,chid_gas,nabad_gas,sound_gas


      double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
                        forth   = 4.0d0/3.0d0, &
                        forpi   = 4.0d0 * pi, &
                        kergavo = kerg * avo, &
                        ikavo   = 1.0d0/kergavo, &
                        asoli3  = asol/3.0d0, &
                        light2  = clight * clight)

! for the abar derivatives
      double precision dpradda,deradda,dsradda, &
                       dpionda,deionda,dsionda, &
                       dpepda,deepda,dsepda, &
                       dpresda,denerda,dentrda, &
                       detada,dxneda

! for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz, &
                       dpiondz,deiondz,dsiondz, &
                       dpepdz,deepdz,dsepdz, &
                       dpresdz,denerdz,dentrdz, &
                       detadz,dxnedz

! for the interpolations
      integer          iat,jat
      double precision free,df_d,df_t,df_dd,df_tt,df_dt
      double precision xt,xd,mxt,mxd, &
                       si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                       si0d,si1d,si2d,si0md,si1md,si2md, &
                       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                       ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                       z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
                       dpsi2,ddpsi2,din,h5,fi(36), &
                       xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
                       w0t,w1t,w2t,w0mt,w1mt,w2mt, &
                       w0d,w1d,w2d,w0md,w1md,w2md


! for the uniform background coulomb correction
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
                       plasg,plasgdd,plasgdt,plasgda,plasgdz, &
                       ecoul,decouldd,decouldt,decoulda,decouldz, &
                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
                       a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
      parameter        (a1    = -0.898004d0, &
                        b1    =  0.96786d0, &
                        c1    =  0.220703d0, &
                        d1    = -0.86097d0, &
                        e1    =  2.5269d0, &
                        a2    =  0.29561d0, &
                        b2    =  1.9885d0, &
                        c2    =  0.288675d0, &
                        third =  1.0d0/3.0d0, &
                        esqu  =  qe * qe)


! quintic hermite polynomial statement functions
! psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


! psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


! psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


! biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
             fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



! cubic hermite polynomial statement functions
! psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)


! psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


! bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
             fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
           + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
           + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
           + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
           + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



! popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))



! start of pipeline loop, normal execution starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

!       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
!       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = max(1.0d-16,ytot1 * zbar)



! initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


! radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0

       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0


! ion section:
        xni     = avo * ytot1 * den
        dxnidd  = avo * ytot1
        dxnida  = -xni * ytot1

        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = dxnida * kt
        dpiondz = 0.0d0

        eion    = 1.5d0 * pion*deni
        deiondd = (1.5d0 * dpiondd - eion)*deni
        deiondt = 1.5d0 * dpiondt*deni
        deionda = 1.5d0 * dpionda*deni
        deiondz = 0.0d0


! sackur-tetrode equation for the ion entropy of
! a single ideal gas characterized by abar
        x       = abar*abar*sqrt(abar) * deni/avo
        s       = sioncon * temp
        z       = x * s * sqrt(s)
        y       = log(z)

!        y       = 1.0d0/(abar*kt)
!        yy      = y * sqrt(y)
!        z       = xni * sifac * yy
!        etaion  = log(z)


        sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
        dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                   - kergavo * deni * ytot1
        dsiondt = (dpiondt*deni + deiondt)*tempi - &
                  (pion*deni + eion) * tempi*tempi &
                  + 1.5d0 * kergavo * tempi*ytot1
        x       = avo*kerg/abar
        dsionda = (dpionda*deni + deionda)*tempi &
                  + kergavo*ytot1*ytot1* (2.5d0 - y)
        dsiondz = 0.0d0



! electron-positron section:


! assume complete ionization
        xnem    = xni * zbar


! enter the table with ye*den
        din = ye*den


! bomb proof the input
        if (temp .gt. t(jmax)) then
         write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
         write(6,*) 'temp too hot, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (temp .lt. t(1)) then
         write(6,01) 'temp=',temp,' t(1)=',t(1)
         write(6,*) 'temp too cold, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .gt. d(imax)) then
         write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
         write(6,*) 'ye*den too big, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .lt. d(1)) then
         write(6,01) 'ye*den=',din,' d(1)=',d(1)
         write(6,*) 'ye*den too small, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if

! hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jmax-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,imax-1))


! access the table locations only once
        fi(1)  = f(iat,jat)
        fi(2)  = f(iat+1,jat)
        fi(3)  = f(iat,jat+1)
        fi(4)  = f(iat+1,jat+1)
        fi(5)  = ft(iat,jat)
        fi(6)  = ft(iat+1,jat)
        fi(7)  = ft(iat,jat+1)
        fi(8)  = ft(iat+1,jat+1)
        fi(9)  = ftt(iat,jat)
        fi(10) = ftt(iat+1,jat)
        fi(11) = ftt(iat,jat+1)
        fi(12) = ftt(iat+1,jat+1)
        fi(13) = fd(iat,jat)
        fi(14) = fd(iat+1,jat)
        fi(15) = fd(iat,jat+1)
        fi(16) = fd(iat+1,jat+1)
        fi(17) = fdd(iat,jat)
        fi(18) = fdd(iat+1,jat)
        fi(19) = fdd(iat,jat+1)
        fi(20) = fdd(iat+1,jat+1)
        fi(21) = fdt(iat,jat)
        fi(22) = fdt(iat+1,jat)
        fi(23) = fdt(iat,jat+1)
        fi(24) = fdt(iat+1,jat+1)
        fi(25) = fddt(iat,jat)
        fi(26) = fddt(iat+1,jat)
        fi(27) = fddt(iat,jat+1)
        fi(28) = fddt(iat+1,jat+1)
        fi(29) = fdtt(iat,jat)
        fi(30) = fdtt(iat+1,jat)
        fi(31) = fdtt(iat,jat+1)
        fi(32) = fdtt(iat+1,jat+1)
        fi(33) = fddtt(iat,jat)
        fi(34) = fddtt(iat+1,jat)
        fi(35) = fddtt(iat,jat+1)
        fi(36) = fddtt(iat+1,jat+1)


! various differences
        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
        mxt = 1.0d0 - xt
        mxd = 1.0d0 - xd

! the six density and six temperature basis functions
        si0t =   psi0(xt)
        si1t =   psi1(xt)*dt_sav(jat)
        si2t =   psi2(xt)*dt2_sav(jat)

        si0mt =  psi0(mxt)
        si1mt = -psi1(mxt)*dt_sav(jat)
        si2mt =  psi2(mxt)*dt2_sav(jat)

        si0d =   psi0(xd)
        si1d =   psi1(xd)*dd_sav(iat)
        si2d =   psi2(xd)*dd2_sav(iat)

        si0md =  psi0(mxd)
        si1md = -psi1(mxd)*dd_sav(iat)
        si2md =  psi2(mxd)*dd2_sav(iat)

! derivatives of the weight functions
        dsi0t =   dpsi0(xt)*dti_sav(jat)
        dsi1t =   dpsi1(xt)
        dsi2t =   dpsi2(xt)*dt_sav(jat)

        dsi0mt = -dpsi0(mxt)*dti_sav(jat)
        dsi1mt =  dpsi1(mxt)
        dsi2mt = -dpsi2(mxt)*dt_sav(jat)

        dsi0d =   dpsi0(xd)*ddi_sav(iat)
        dsi1d =   dpsi1(xd)
        dsi2d =   dpsi2(xd)*dd_sav(iat)

        dsi0md = -dpsi0(mxd)*ddi_sav(iat)
        dsi1md =  dpsi1(mxd)
        dsi2md = -dpsi2(mxd)*dd_sav(iat)

! second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
        ddsi2t =   ddpsi2(xt)

        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

!        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
!        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
!        ddsi2d =   ddpsi2(xd)

!        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
!        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
!        ddsi2md =  ddpsi2(mxd)


! the free energy
        free  = h5(iat,jat, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density
        df_d  = h5(iat,jat, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


! derivative with respect to temperature
        df_t = h5(iat,jat, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density**2
!        df_dd = h5(iat,jat,
!     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
!     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

! derivative with respect to temperature**2
        df_tt = h5(iat,jat, &
              ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to temperature and density
        df_dt = h5(iat,jat, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



! now get the pressure derivative with density, chemical potential, and
! electron positron number densities
! get the interpolation weight functions
        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt)*dt_sav(jat)

        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt)*dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd)*dd_sav(iat)

        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd)*dd_sav(iat)


! derivatives of weight functions
        dsi0t  = xdpsi0(xt)*dti_sav(jat)
        dsi1t  = xdpsi1(xt)

        dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
        dsi1mt = xdpsi1(mxt)

        dsi0d  = xdpsi0(xd)*ddi_sav(iat)
        dsi1d  = xdpsi1(xd)

        dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
        dsi1md = xdpsi1(mxd)


! look in the pressure derivative only once
        fi(1)  = dpdf(iat,jat)
        fi(2)  = dpdf(iat+1,jat)
        fi(3)  = dpdf(iat,jat+1)
        fi(4)  = dpdf(iat+1,jat+1)
        fi(5)  = dpdft(iat,jat)
        fi(6)  = dpdft(iat+1,jat)
        fi(7)  = dpdft(iat,jat+1)
        fi(8)  = dpdft(iat+1,jat+1)
        fi(9)  = dpdfd(iat,jat)
        fi(10) = dpdfd(iat+1,jat)
        fi(11) = dpdfd(iat,jat+1)
        fi(12) = dpdfd(iat+1,jat+1)
        fi(13) = dpdfdt(iat,jat)
        fi(14) = dpdfdt(iat+1,jat)
        fi(15) = dpdfdt(iat,jat+1)
        fi(16) = dpdfdt(iat+1,jat+1)

! pressure derivative with density
        dpepdd  = h3(iat,jat, &
                       si0t,   si1t,   si0mt,   si1mt, &
                       si0d,   si1d,   si0md,   si1md)
        dpepdd  = max(ye * dpepdd,1.0d-30)



! look in the electron chemical potential table only once
        fi(1)  = ef(iat,jat)
        fi(2)  = ef(iat+1,jat)
        fi(3)  = ef(iat,jat+1)
        fi(4)  = ef(iat+1,jat+1)
        fi(5)  = eft(iat,jat)
        fi(6)  = eft(iat+1,jat)
        fi(7)  = eft(iat,jat+1)
        fi(8)  = eft(iat+1,jat+1)
        fi(9)  = efd(iat,jat)
        fi(10) = efd(iat+1,jat)
        fi(11) = efd(iat,jat+1)
        fi(12) = efd(iat+1,jat+1)
        fi(13) = efdt(iat,jat)
        fi(14) = efdt(iat+1,jat)
        fi(15) = efdt(iat,jat+1)
        fi(16) = efdt(iat+1,jat+1)


! electron chemical potential etaele
        etaele  = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)


! derivative with respect to density
        x       = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
        detadd  = ye * x

! derivative with respect to temperature
        detadt  = h3(iat,jat, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
       detada = -x * din * ytot1
       detadz =  x * den * ytot1



! look in the number density table only once
        fi(1)  = xf(iat,jat)
        fi(2)  = xf(iat+1,jat)
        fi(3)  = xf(iat,jat+1)
        fi(4)  = xf(iat+1,jat+1)
        fi(5)  = xft(iat,jat)
        fi(6)  = xft(iat+1,jat)
        fi(7)  = xft(iat,jat+1)
        fi(8)  = xft(iat+1,jat+1)
        fi(9)  = xfd(iat,jat)
        fi(10) = xfd(iat+1,jat)
        fi(11) = xfd(iat,jat+1)
        fi(12) = xfd(iat+1,jat+1)
        fi(13) = xfdt(iat,jat)
        fi(14) = xfdt(iat+1,jat)
        fi(15) = xfdt(iat,jat+1)
        fi(16) = xfdt(iat+1,jat+1)

! electron + positron number densities
       xnefer   = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to density
       x        = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
       x = max(x,1.0d-30)
       dxnedd   = ye * x

! derivative with respect to temperature
       dxnedt   = h3(iat,jat, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
       dxneda = -x * din * ytot1
       dxnedz =  x  * den * ytot1


! the desired electron-positron thermodynamic quantities

! dpepdd at high temperatures and low densities is below the
! floating point limit of the subtraction of two large terms.
! since dpresdd doesn't enter the maxwell relations at all, use the
! bicubic interpolation done above instead of the formally correct expression
        x       = din * din
        pele    = x * df_d
        dpepdt  = x * df_dt
!        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
        s       = dpepdd/ye - 2.0d0 * din * df_d
        dpepda  = -ytot1 * (2.0d0 * pele + s * din)
        dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


        x       = ye * ye
        sele    = -df_t * ye
        dsepdt  = -df_tt * ye
        dsepdd  = -df_dt * x
        dsepda  = ytot1 * (ye * df_dt * din - sele)
        dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


        eele    = ye*free + temp * sele
        deepdt  = temp * dsepdt
        deepdd  = x * df_d + temp * dsepdd
        deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
        deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




! coulomb section:

! uniform background corrections only
! from yakovlev & shalybkov 1989
! lami is the average ion seperation
! plasg is the plasma coupling parameter

        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd
        dsda     = z * dxnida

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s
        lamida   = z * dsda/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami
        plasgdd  = z * lamidd
        plasgda  = z * lamida
        plasgdt  = -plasg*ktinv * kerg
        plasgdz  = 2.0d0 * plasg/zbar


! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0) then
         x        = plasg**(0.25d0)
         y        = avo * ytot1 * kerg
         ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
         pcoul    = third * den * ecoul
         scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
                    + d1 * (log(plasg) - 1.0d0) - e1)

         y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
         decouldd = y * plasgdd
         decouldt = y * plasgdt + ecoul/temp
         decoulda = y * plasgda - ecoul/abar
         decouldz = y * plasgdz

         y        = third * den
         dpcouldd = third * ecoul + y*decouldd
         dpcouldt = y * decouldt
         dpcoulda = y * decoulda
         dpcouldz = y * decouldz


         y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
         dscouldd = y * plasgdd
         dscouldt = y * plasgdt
         dscoulda = y * plasgda - scoul/abar
         dscouldz = y * plasgdz


! yakovlev & shalybkov 1989 equations 102, 103, 104
        else if (plasg .lt. 1.0) then
         x        = plasg*sqrt(plasg)
         y        = plasg**b2
         z        = c2 * x - third * a2 * y
         pcoul    = -pion * z
         ecoul    = 3.0d0 * pcoul/den
         scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

         s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
         dpcouldd = -dpiondd*z - pion*s*plasgdd
         dpcouldt = -dpiondt*z - pion*s*plasgdt
         dpcoulda = -dpionda*z - pion*s*plasgda
         dpcouldz = -dpiondz*z - pion*s*plasgdz

         s        = 3.0d0/den
         decouldd = s * dpcouldd - ecoul/den
         decouldt = s * dpcouldt
         decoulda = s * dpcoulda
         decouldz = s * dpcouldz

         s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
         dscouldd = s * plasgdd
         dscouldt = s * plasgdt
         dscoulda = s * plasgda - scoul/abar
         dscouldz = s * plasgdz
        end if


! bomb proof
        x   = prad + pion + pele + pcoul
        y   = erad + eion + eele + ecoul
        z   = srad + sion + sele + scoul

!        write(6,*) x,y,z
!        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then
        if (x .le. 0.0 .or. y .le. 0.0) then
!        if (x .le. 0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if


! sum all the gas components
       pgas    = pion + pele + pcoul
       egas    = eion + eele + ecoul
       sgas    = sion + sele + scoul

       dpgasdd = dpiondd + dpepdd + dpcouldd
       dpgasdt = dpiondt + dpepdt + dpcouldt
       dpgasda = dpionda + dpepda + dpcoulda
       dpgasdz = dpiondz + dpepdz + dpcouldz

       degasdd = deiondd + deepdd + decouldd
       degasdt = deiondt + deepdt + decouldt
       degasda = deionda + deepda + decoulda
       degasdz = deiondz + deepdz + decouldz

       dsgasdd = dsiondd + dsepdd + dscouldd
       dsgasdt = dsiondt + dsepdt + dscouldt
       dsgasda = dsionda + dsepda + dscoulda
       dsgasdz = dsiondz + dsepdz + dscouldz




! add in radiation to get the total
       pres    = prad + pgas
       ener    = erad + egas
       entr    = srad + sgas

       dpresdd = dpraddd + dpgasdd
       dpresdt = dpraddt + dpgasdt
       dpresda = dpradda + dpgasda
       dpresdz = dpraddz + dpgasdz

       denerdd = deraddd + degasdd
       denerdt = deraddt + degasdt
       denerda = deradda + degasda
       denerdz = deraddz + degasdz

       dentrdd = dsraddd + dsgasdd
       dentrdt = dsraddt + dsgasdt
       dentrda = dsradda + dsgasda
       dentrdz = dsraddz + dsgasdz


! for the gas
! the temperature and density exponents (c&g 9.81 9.82)
! the specific heat at constant volume (c&g 9.92)
! the third adiabatic exponent (c&g 9.93)
! the first adiabatic exponent (c&g 9.97)
! the second adiabatic exponent (c&g 9.105)
! the specific heat at constant pressure (c&g 9.98)
! and relativistic formula for the sound speed (c&g 14.29)

       zz        = pgas*deni
       zzi       = den/pgas
       chit_gas  = temp/pgas * dpgasdt
       chid_gas  = dpgasdd*zzi
       cv_gas    = degasdt
       x         = zz * chit_gas/(temp * cv_gas)
       gam3_gas  = x + 1.0d0
       gam1_gas  = chit_gas*x + chid_gas
       nabad_gas = x/gam1_gas
       gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
       cp_gas    = cv_gas * gam1_gas/chid_gas
       z         = 1.0d0 + (egas + light2)*zzi
       sound_gas = clight * sqrt(gam1_gas/z)



! for the totals
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)



! maxwell relations; each is zero if the consistency is perfect
       x   = den * den

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*x/dpresdt - 1.0d0


! store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda
        dez_row(j)    = denerdz

        stot_row(j)   = entr
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda
        dsz_row(j)    = dentrdz


        pgas_row(j)   = pgas
        dpgast_row(j) = dpgasdt
        dpgasd_row(j) = dpgasdd
        dpgasa_row(j) = dpgasda
        dpgasz_row(j) = dpgasdz

        egas_row(j)   = egas
        degast_row(j) = degasdt
        degasd_row(j) = degasdd
        degasa_row(j) = degasda
        degasz_row(j) = degasdz

        sgas_row(j)   = sgas
        dsgast_row(j) = dsgasdt
        dsgasd_row(j) = dsgasdd
        dsgasa_row(j) = dsgasda
        dsgasz_row(j) = dsgasdz


        prad_row(j)   = prad
        dpradt_row(j) = dpraddt
        dpradd_row(j) = dpraddd
        dprada_row(j) = dpradda
        dpradz_row(j) = dpraddz

        erad_row(j)   = erad
        deradt_row(j) = deraddt
        deradd_row(j) = deraddd
        derada_row(j) = deradda
        deradz_row(j) = deraddz

        srad_row(j)   = srad
        dsradt_row(j) = dsraddt
        dsradd_row(j) = dsraddd
        dsrada_row(j) = dsradda
        dsradz_row(j) = dsraddz


        pion_row(j)   = pion
        dpiont_row(j) = dpiondt
        dpiond_row(j) = dpiondd
        dpiona_row(j) = dpionda
        dpionz_row(j) = dpiondz

        eion_row(j)   = eion
        deiont_row(j) = deiondt
        deiond_row(j) = deiondd
        deiona_row(j) = deionda
        deionz_row(j) = deiondz

        sion_row(j)   = sion
        dsiont_row(j) = dsiondt
        dsiond_row(j) = dsiondd
        dsiona_row(j) = dsionda
        dsionz_row(j) = dsiondz

        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda
        deepz_row(j)  = deepdz

        sele_row(j)   = sele
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = dsepda
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xnem
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxnedt
        dxned_row(j)  = dxnedd
        dxnea_row(j)  = dxneda
        dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0
        zeff_row(j)   = zbar

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = 0.0d0

        pcou_row(j)   = pcoul
        dpcout_row(j) = dpcouldt
        dpcoud_row(j) = dpcouldd
        dpcoua_row(j) = dpcoulda
        dpcouz_row(j) = dpcouldz

        ecou_row(j)   = ecoul
        decout_row(j) = decouldt
        decoud_row(j) = decouldd
        decoua_row(j) = decoulda
        decouz_row(j) = decouldz

        scou_row(j)   = scoul
        dscout_row(j) = dscouldt
        dscoud_row(j) = dscouldd
        dscoua_row(j) = dscoulda
        dscouz_row(j) = dscouldz

        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_gas_row(j)    = cv_gas
        cp_gas_row(j)    = cp_gas
        gam1_gas_row(j)  = gam1_gas
        gam2_gas_row(j)  = gam2_gas
        gam3_gas_row(j)  = gam3_gas
        nabad_gas_row(j) = nabad_gas
        cs_gas_row(j)    = sound_gas

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        nabad_row(j)  = nabad
        cs_row(j)     = sound

! end of pipeline loop
      enddo
      return
      end






      subroutine pretty_eos_out(whose)
      include 'implno.dek'
      include 'vector_eos.dek'

! writes a pretty output for the eos tester


! declare the pass
      character*(*) whose


! local variables
      integer          i,j
      double precision ye,xcess,avo,kerg,xka
      parameter        (avo     = 6.0221417930d23, &
                        kerg    = 1.380650424d-16, &
                        xka = kerg*avo)


! popular formats
01    format(1x,t2,a,t11,a,t27,a,t43,a,t59,a,t75,a,t91,a,t107,a)
02    format(1x,t2,a,1p7e16.8)
03    format(1x,t2,a7,1pe12.4,t22,a7,1pe12.4, &
               t42,a7,1pe12.4,t62,a7,1pe12.4)
04    format(1x,t2,a,t11,'total',t24,'ion',t34,'e- + e+', &
             t58,'radiation',t70,'coulomb')
05    format(1x,t2,a,1p3e12.4,t56,1p2e12.4)
06    format(1x,t2,a,a,1pe12.4, &
                t30,a,a,1pe12.4, &
                t58,a,a,1pe12.4)



! loop over the pipeline
      do j=jlo_eos,jhi_eos


! the input
      write(6,03) 'temp  =',temp_row(j),'den   =',den_row(j), &
                  'abar  =',abar_row(j),'zbar  =',zbar_row(j)

      ye = zbar_row(1)/abar_row(1)
      xcess = 1.0d0 - 2.0d0*ye
      write(6,03) 'ye    =',ye,'xcess =',xcess
      write(6,*) ' '


! and the output

       write(6,01)  whose,'value','d/dd','d/dt','d/da','d/dz'

       write(6,02) 'p tot=',ptot_row(j), &
                    dpd_row(j),dpt_row(j),dpa_row(j),dpz_row(j)
       write(6,02) 'p gas=',pgas_row(j), &
                 dpgasd_row(j),dpgast_row(j),dpgasa_row(j),dpgasz_row(j)
       write(6,02) 'p rad=',prad_row(j), &
                 dpradd_row(j),dpradt_row(j),dprada_row(j),dpradz_row(j)
       write(6,02) 'p ion=',pion_row(j), &
                dpiond_row(j),dpiont_row(j),dpiona_row(j),dpionz_row(j)
       write(6,02) 'p  e-=',pele_row(j), &
                dpepd_row(j),dpept_row(j),dpepa_row(j),dpepz_row(j)
       write(6,02) 'p  e+=',ppos_row(j)
       write(6,02) 'p cou=',pcou_row(j), &
                dpcoud_row(j),dpcout_row(j),dpcoua_row(j),dpcouz_row(j)


       write(6,*)  ' '
       write(6,02) 'e tot=',etot_row(j), &
                    ded_row(j),det_row(j),dea_row(j),dez_row(j)
       write(6,02) 'e gas=',egas_row(j), &
                 degasd_row(j),degast_row(j),degasa_row(j),degasz_row(j)
       write(6,02) 'e rad=',erad_row(j), &
                deradd_row(j),deradt_row(j),derada_row(j),deradz_row(j)
       write(6,02) 'e ion=',eion_row(j), &
                deiond_row(j),deiont_row(j),deiona_row(j),deionz_row(j)
       write(6,02) 'e  e-=',eele_row(j), &
                deepd_row(j),deept_row(j),deepa_row(j),deepz_row(j)
       write(6,02) 'e  e+=',epos_row(j)
       write(6,02) 'e cou=',ecou_row(j), &
                decoud_row(j),decout_row(j),decoua_row(j),decouz_row(j)

       write(6,*)  ' '
       write(6,02) 's tot=',stot_row(j), &
                    dsd_row(j),dst_row(j),dsa_row(j),dsz_row(j)
       write(6,02) 's/xka=',stot_row(j)/xka, &
             dsd_row(j)/xka,dst_row(j)/xka,dsa_row(j)/xka,dsz_row(j)/xka
       write(6,02) 's gas=',sgas_row(j), &
                 dsgasd_row(j),dsgast_row(j),dsgasa_row(j),dsgasz_row(j)
       write(6,02) 's rad=',srad_row(j), &
                dsradd_row(j),dsradt_row(j),dsrada_row(j),dsradz_row(j)
       write(6,02) 's ion=',sion_row(j), &
                dsiond_row(j),dsiont_row(j),dsiona_row(j),dsionz_row(j)
       write(6,02) 's  e-=',sele_row(j), &
                dsepd_row(j),dsept_row(j),dsepa_row(j),dsepz_row(j)
       write(6,02) 's  e+=',spos_row(j)
       write(6,02) 's cou=',scou_row(j), &
                dscoud_row(j),dscout_row(j),dscoua_row(j),dscouz_row(j)


! specific heats, and ratio of electostatic to thermal energy
! the 3 gammas and the sound speed for both the gas and the total
       write(6,*)  ' '
       write(6,02) 'cv  =',cv_row(j)/(kerg*avo)*abar_row(1), &
                    dcvdd_row(j),dcvdt_row(j), &
                    dcvda_row(j),dcvdz_row(j)
       write(6,02) 'cp  =',cp_row(j), &
                    dcpdd_row(j),dcpdt_row(j), &
                    dcpda_row(j),dcpdz_row(j)
       write(6,02) 'gam1=',gam1_row(j), &
                    dgam1dd_row(j),dgam1dt_row(j), &
                    dgam1da_row(j),dgam1dz_row(j)
       write(6,02) 'gam2=',gam2_row(j), &
                    dgam2dd_row(j),dgam2dt_row(j), &
                    dgam2da_row(j),dgam2dz_row(j)
       write(6,02) 'gam3=',gam3_row(j), &
                    dgam3dd_row(j),dgam3dt_row(j), &
                    dgam3da_row(j),dgam3dz_row(j)
       write(6,02) 'cs  =',cs_row(j), &
                    dcsdd_row(j),dcsdt_row(j), &
                    dcsda_row(j),dcsdz_row(j)

       write(6,*)  ' '
       write(6,02) 'cvgas=',cv_gas_row(j)/(kerg*avo)*abar_row(1), &
                    dcv_gasdd_row(j),dcv_gasdt_row(j), &
                    dcv_gasda_row(j),dcv_gasdz_row(j)
       write(6,02) 'cpgas=',cp_gas_row(j), &
                    dcp_gasdd_row(j),dcp_gasdt_row(j), &
                    dcp_gasda_row(j),dcp_gasdz_row(j)
       write(6,02) 'g1gas=',gam1_gas_row(j), &
                    dgam1_gasdd_row(j),dgam1_gasdt_row(j), &
                    dgam1_gasda_row(j),dgam1_gasdz_row(j)
       write(6,02) 'g2gas=',gam2_gas_row(j), &
                    dgam2_gasdd_row(j),dgam2_gasdt_row(j), &
                    dgam2_gasda_row(j),dgam2_gasdz_row(j)
       write(6,02) 'g3gas=',gam3_gas_row(j), &
                    dgam3_gasdd_row(j),dgam3_gasdt_row(j), &
                    dgam3_gasda_row(j),dgam3_gasdz_row(j)
       write(6,02) 'csgas=',cs_gas_row(j), &
                    dcs_gasdd_row(j),dcs_gasdt_row(j), &
                    dcs_gasda_row(j),dcs_gasdz_row(j)


! the thermodynamic consistency relations, these should all be
! at the floating point limit of zero
       write(6,*) ' '
       write(6,03) 'maxw1 =',dse_row(j),'maxw2 =',dpe_row(j), &
                   'maxw3 =',dsp_row(j)

! number density of ions and its derivatives
       write(6,03) 'xni   =',xni_row(j),  'xnim  =',xnim_row(j)
       write(6,03) 'dxnidd=',dxned_row(j),'dxnidt=',dxnet_row(j), &
                   'dxnida=',dxnea_row(j),'dxnidz=',dxnez_row(j)

! ion chemical potential and its derivatives
       write(6,03) 'etaion=',etaion_row(j)
       write(6,03) 'detaid=',detaid_row(j),'detait=',detait_row(j), &
                   'detaia=',detaia_row(j),'detaiz=',detaiz_row(j)


! number density of electrons+positrons and its derivatives
       write(6,03) 'xnele =',xne_row(j),'xnpos =',xnp_row(j), &
                   'xnem  =',xnem_row(j)
       write(6,03) 'dxnedd=',dxned_row(j),'dxnedt=',dxnet_row(j), &
                   'dxneda=',dxnea_row(j),'dxnedz=',dxnez_row(j)


! electron chemical potential, positron chemical potential and its derivatives
       write(6,03) 'etaele=',etaele_row(j),'etapos=',etapos_row(j)
       write(6,03) 'detadd=',detad_row(j),'detadt=',detat_row(j), &
                   'detada=',detaa_row(j),'detadz=',detaz_row(j)

       write(6,03) 'zeff  =',zeff_row(j), &
                   'ionzd =',zeff_row(j)/zbar_row(j), &
                   'plasg =',plasg_row(j)

! end of pipeline loop
      enddo

      return
      end

