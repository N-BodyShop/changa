  program hydro
    use eos_def
    use eos_lib
    use chem_def
    use chem_lib
    use const_lib
    use crlibm_lib


    implicit none

    double precision :: massCentral, rCentral
    double precision, dimension(:), allocatable :: mass, logR, logRho, logT, entropy, pressure, temp, logP, dlogTdlogP, HFraction
    double precision :: rho, T, P, energy, gamma, mstar, rstar, mtarget, rtarget, mprev, rprev, junk, h
    double precision :: T2, rho2, energy2, gamma2, S, T_guess, Xc, Xin
    double precision :: P_start, P_test, dP_test, P_prev, rhoc
    double precision :: Tc_start, Tc_test, dTc_test, Tc_prev

    double precision, parameter :: G=6.67d-8, msolar=2d33, rsolar=7d10, kB=1.38d-16, mproton=1.6d-24
    integer, parameter :: noeqs = 3
    double precision, parameter :: mtol = 1d-2, rtol= 1d-2
    integer :: bins, i, firstP, firstT
    integer :: handle
    read(*,*) massCentral, rCentral
    write(*,"(2(E16.10,1X))") massCentral, rCentral

    read(*,*) bins
    allocate( mass(bins), logR(bins), logrho(bins), logT(bins), logP(bins))
    allocate( entropy(bins), pressure(bins), temp(bins), dlogTdlogP(bins), HFraction(bins))
    do i = 1, bins
      read(*,*) mass(i), logR(i), logRho(i), logT(i), logP(i), dlogTdlogP(i), HFraction(i)
      HFraction(i) = 0.73D0
    enddo

    massCentral = (massCentral)*msolar
    mass = mass*msolar
    rCentral = rsolar*rCentral
    h = rCentral
    call mesaeos_init

    T_guess = 1e5
    do i = 1, bins
      rho = 1e1**logRho(i)
      T = 1e1**logT(i)
      temp(i) = T
      call find_HFraction( mass(i), Xin)
      call mesaeos_DTget( rho, T, Xin, pressure(i), energy, entropy(i), gamma)
      !call mesaeos_PSget( pressure(i), entropy(i), T_guess, rho2, energy2, T2, S, gamma2)
      !write(*,*) mass(i), T2/T-1D0, T_guess, rho2/rho - 1D0
    enddo

    rstar = 0.
    rtarget = 1d1**logR(bins)*rsolar
    mtarget = mass(bins)

    Tc_test = 0.5*G*massCentral/(rCentral)*mp/kB
    if( massCentral <= 0.) then
      Tc_test = 2e7
    endif

    Tc_prev = 0.
    rstar = 0.
    rprev = rstar
    P_test = pressure(1)
    P_prev = 0.
    mstar = 0.
    mprev = mstar

    rhoc = 1e1**logRho(1)
    firstT = 1

    do while( abs(rstar/rtarget - 1D0) > rtol )

      firstP = 1
      !P_test = pressure(1)
      P_prev = 0.
      mstar = 0.
      mprev = mstar
      call find_HFraction( mass(1), Xc)
      call mesaeos_DTget( rhoc, Tc_test, Xc, P_test, energy, junk, gamma)
      write(0,*) 'rtol=', rstar/rtarget - 1D0
      do while( abs(mstar/mtarget-1D0) > mtol)
      !

        call integrate_star( P_test, Tc_test, mstar, rstar, 0)
        call mesaeos_PTget( P_test, Tc_test, Xc, rhoc, energy, junk, gamma)

        write(0,*) 'mtol=', mstar/mtarget - 1.

        if( abs(mstar/mtarget-1D0) > mtol) then
          if (firstP .eq. 1) then
            firstP = 0
            dP_test = P_test*0.01
            P_prev = P_test
            mprev = mstar
            P_test = P_test + dP_test
          else
            dP_test = (P_test - P_prev)/(mstar-mprev)*(mtarget-mstar)
            if(dP_test > 0.) then
              dP_test = min(P_test*0.25, dP_test)
            endif
            if(dP_test < 0.) then
              dP_test = max(-P_test*0.1, dP_test)
            endif
            P_prev = P_test
            mprev = mstar
            P_test = P_test + dP_test
          endif
        endif
      enddo
      if( abs(rstar/rtarget - 1D0) > rtol ) then
        if (firstT .eq. 1) then
          firstT = 0
          dTc_test = Tc_test*0.01
          Tc_prev = Tc_test
          rprev = rstar
          Tc_test = Tc_test + dTc_test
        else
          dTc_test = (Tc_test - Tc_prev)/(rstar-rprev)*(rtarget-rstar)
          if(dTc_test > 0.) then
            dTc_test = min(Tc_test*0.2, dTc_test)
          endif
          if(dTc_test < 0.) then
            dTc_test = max(-Tc_test*0.1, dTc_test)
          endif
          Tc_prev = Tc_test
          rprev = rstar
          Tc_test = Tc_test + dTc_test
        endif
      endif
      !write(*,*) "Found star", rstar/rtarget, mstar/mtarget, Tc_test, P_test


    enddo
    !call mesaeos_DTget( 1e1**logRho(1), T_test, Xc, P, energy, junk, gamma)
    call integrate_star( P_test, Tc_test, mstar, rstar, 1)


    call mesaeos_shutdown
  contains
    subroutine mesaeos_init
      implicit none
      character (len=256) :: eos_file_prefix
      integer :: ierr
      logical, parameter :: use_cache = .true.
      character (len=255) :: my_mesa_dir

      integer :: i

      ierr = 0
      call get_environment_variable("MESA_DIR", my_mesa_dir)
      call const_init(my_mesa_dir,ierr)
      if (ierr /= 0) then
        write(*,*) 'const_init failed'
        stop 1
      end if

      call crlibm_init

      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
        write(*,*) 'failed in chem_init'
        stop 1
      end if

      eos_file_prefix = 'mesa'

      call eos_init(eos_file_prefix, '', '', '', use_cache, ierr)
      if (ierr /= 0) then
        write(*,*) 'eos_init failed in Setup_eos'
        stop 1
      end if

      write(*,*) 'loading eos tables'

      handle = alloc_eos_handle(ierr)
      if (ierr /= 0) then
        write(*,*) 'failed trying to allocate eos handle'
        stop 1
      end if
    end subroutine mesaeos_init

    subroutine mesaeos_shutdown()
      implicit none
      call free_eos_handle(handle)
      call eos_shutdown
    end subroutine mesaeos_shutdown

    subroutine mesaeos_chem( num_chem_isos, species, net_iso, chem_id, Xin, X, Y, Z, abar, zbar, xa)
      implicit none
      integer, intent(in) :: num_chem_isos, species
      double precision, intent(IN) :: Xin
      double precision, intent(OUT) :: X, Z, Y, abar, zbar
      double precision :: z2bar, ye
      integer, pointer, dimension(:) :: net_iso
      integer, pointer, dimension(:) :: chem_id
      double precision, dimension(species) :: xa
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
      double precision :: xz, frac, dabar_dx(species), dzbar_dx(species), sumx, mass_correction, dmc_dx(species)

      double precision, parameter :: Zfrac_C = 0.173312d0
      double precision, parameter :: Zfrac_N = 0.053177d0
      double precision, parameter :: Zfrac_O = 0.482398d0
      double precision, parameter :: Zfrac_Ne = 0.098675d0

      X = Xin
      Z = 0.02


      net_iso(:) = 0

      chem_id(h1) = ih1; net_iso(ih1) = h1
      chem_id(he4) = ihe4; net_iso(ihe4) = he4
      chem_id(c12) = ic12; net_iso(ic12) = c12
      chem_id(n14) = in14; net_iso(in14) = n14
      chem_id(o16) = io16; net_iso(io16) = o16
      chem_id(ne20) = ine20; net_iso(ine20) = ne20
      chem_id(mg24) = img24; net_iso(img24) = mg24

      Y = 1 - (X + Z)

      xa(h1) = X
      xa(he4) = Y
      xa(c12) = Z * Zfrac_C
      xa(n14) = Z * Zfrac_N
      xa(o16) = Z * Zfrac_O
      xa(ne20) = Z * Zfrac_Ne
      xa(species) = 1 - sum(xa(1:species-1))

      call composition_info( &
            species, chem_id, xa, X, Y, xz, abar, zbar, z2bar, ye, mass_correction, &
            sumx, dabar_dx, dzbar_dx, dmc_dx)

      return
    end subroutine mesaeos_chem

    subroutine mesaeos_DTget( Rho, T, Xin,press, energy, entropy, gamma)
       implicit none

       double precision :: X, Z, Y, abar, zbar
       integer, parameter :: species = 7
       integer, pointer, dimension(:) :: net_iso, chem_id
       double precision :: xa(species)

       double precision, intent(in) :: Rho, T, Xin
       double precision, intent(out) :: energy, press, entropy, gamma
       double precision :: Pgas,log10Rho
       double precision, dimension(num_eos_basic_results) :: res,d_dlnd, d_dlnT, d_dabar, d_dzbar
       integer :: ierr
       double precision :: log10T

       allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
       if (ierr /= 0) stop 'allocate failed'

       call mesaeos_chem( num_chem_isos, species, net_iso, chem_id, Xin, X, Y, Z, abar, zbar, xa)


       ! get a set of results for given energy and density
       call eosDT_get( &
             handle, Z, X, abar, zbar, &
             species, chem_id, net_iso, xa, &
             Rho, log10_cr(Rho), T, log10_cr(T), &
             res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)

       ! the indices for the results are defined in eos_def.f
       press = exp_cr(res(i_lnPgas)) !+ Radiation_Pressure(T)
       gamma = res(i_gamma1)
       energy = exp_cr(res(i_lnE))
       entropy = exp_cr(res(i_lnS))

       deallocate(net_iso, chem_id)

       if (ierr /= 0) then
          write(*,*) 'bad result from eos_get DT'
          stop 1
       end if

    end subroutine mesaeos_DTget

    subroutine mesaeos_PTget( press, T, Xin, rho, energy, entropy, gamma)
       implicit none

       double precision :: X, Z, Y, abar, zbar
       integer, parameter :: species = 7
       integer, pointer, dimension(:) :: net_iso, chem_id
       double precision :: xa(species)

       double precision, intent(in) :: press, T, Xin
       double precision, intent(out) :: energy, rho, entropy, gamma
       double precision :: Pgas,log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
       double precision, dimension(num_eos_basic_results) :: res, &
                    d_dlnRho_const_T, d_dlnT_const_Rho, d_dabar_const_TRho, d_dzbar_const_TRho

       integer :: ierr
       double precision :: log10T

       allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
       if (ierr /= 0) stop 'allocate failed'

       call mesaeos_chem( num_chem_isos, species, net_iso, chem_id, Xin, X, Y, Z, abar, zbar, xa)

       Pgas = press !- Radiation_Pressure(T)
       ! get a set of results for given energy and density
       call eosPT_get( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pgas, log10_cr(Pgas), T, log10_cr(T), &
               rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_const_TRho, d_dzbar_const_TRho, ierr)

       ! the indices for the results are defined in eos_def.f
       gamma = res(i_gamma1)
       energy = exp_cr(res(i_lnE))
       entropy = exp_cr(res(i_lnS))

       deallocate(net_iso, chem_id)

       if (ierr /= 0) then
          write(*,*) 'bad result from eos_get PT', press, T
          stop 1
       end if

    end subroutine mesaeos_PTget

    subroutine mesaeos_PSget( press, entropy, T_guess, Xin, rho, energy, T, newEntropy, gamma)
       implicit none

       double precision :: X, Z, Y, abar, zbar
       integer, parameter :: species = 7
       integer, pointer, dimension(:) :: net_iso, chem_id
       double precision :: xa(species)

       double precision, intent(in) :: press, entropy, T_guess, Xin
       double precision, intent(out) :: energy, rho, T, gamma, newEntropy
       double precision :: Pgas, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
       double precision, dimension(num_eos_basic_results) :: res, &
              d_dlnRho_const_T, d_dlnT_const_Rho, &
              d_dabar_const_TRho, d_dzbar_const_TRho
       integer :: ierr
       integer :: max_iter = 1000
       double precision :: log10T
       integer :: which_other, eos_calls
       double precision :: other_value, other_tol
       double precision, parameter :: logT_tol = 1e-3
       double precision :: s1, s2, lnS1, lnS2, logT_result
       double precision :: rho1, rho2, e1, e2, g1, g2
       double precision :: logT_guess, logT1, logT2
       logT_guess = log10_cr(T_guess)
       logT1 = 1D0
       logT2 = 7D0

       allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
       if (ierr /= 0) stop 'allocate failed'
       call mesaeos_chem( num_chem_isos, species, net_iso, chem_id, Xin, X, Y, Z, abar, zbar, xa)

       call mesaeos_PTget( press, 1e1**logT1, Xin, rho1, e1, s1, g1)
       call mesaeos_PTget( press, 1e1**logT2, Xin, rho2, e2, s2, g2)

       lnS1 = log_cr(s1)
       lnS2 = log_cr(s2)

       which_other = i_lnS
       other_value = log_cr(entropy)
       other_tol = 1e-3

       call eosPT_get_T( &
                handle, Z, X, abar, zbar, &
                species, chem_id, net_iso, xa, &
                log10_cr(press), which_other, other_value, &
                logT_tol, other_tol, max_iter, logT_guess, &
                logT1, logT2, s1, s2, &
                logT_result, rho, log10Rho, &
                dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
                res, d_dlnRho_const_T, d_dlnT_const_Rho, &
                d_dabar_const_TRho, d_dzbar_const_TRho, &
                eos_calls, ierr)

      if (ierr /= 0) then
        write(*,*) 'bad result from eos_get PS'
        stop 1
      end if

       deallocate(net_iso, chem_id)
       T = 1e1**logT_result

       call mesaeos_PTget( press, T, Xin, rho, energy, newEntropy, gamma)



    end subroutine mesaeos_PSget

    subroutine find_entropy( m, s, grad)
      implicit none
      double precision, intent(IN) :: m
      double precision, intent(OUT) :: s, grad
      integer :: i

      s = 0.
      if( m < mass(1)) then
        s = entropy(1)
        grad = dlogTdlogP(1)
        return
      endif

      if( m > mass(bins)) then
        s = entropy(bins)
        grad = dlogTdlogP(bins)
        return
      endif

      do i = 1, bins-1
        if( m .le. mass(i+1) .and. m .ge. mass(i) ) then
          s = entropy(i) + (entropy(i+1)-entropy(i))/(mass(i+1)-mass(i))*(m-mass(i))
          grad = dlogTdlogP(i) + (dlogTdlogP(i+1)-dlogTdlogP(i))/(mass(i+1)-mass(i))*(m-mass(i))
          return
        endif
      end do

      return
    end subroutine find_entropy

    subroutine find_Hfraction( m, X)
      implicit none
      double precision, intent(IN) :: m
      double precision, intent(OUT) :: X
      integer :: i

      s = 0.
      if( m < mass(1)) then
        X = HFraction(1)
        return
      endif

      if( m > mass(bins)) then
        X = HFraction(bins)
        return
      endif

      do i = 1, bins-1
        if( m .le. mass(i+1) .and. m .ge. mass(i) ) then
          X = HFraction(i) + (HFraction(i+1)-HFraction(i))/(mass(i+1)-mass(i))*(m-mass(i))
          return
        endif
      end do

      return
    end subroutine find_HFraction

    subroutine integrate_star( Pb, Tc, mstar, rstar, write_out)
      implicit none
      double precision, intent(IN) :: Pb, Tc
      integer, intent(IN) :: write_out
      double precision, intent(OUT) :: mstar, rstar

      double precision, dimension(noeqs) :: y, dydx, height
      double precision :: Sc, p, grad, Xc, Xin
      double precision :: m, dr, r, rho, e, T, s, gamma, m2

      call find_entropy( massCentral, Sc, grad)
      call find_Hfraction( massCentral, Xc)
      p = Pb
      m = 0D0
      dr = max(rCentral,1e-2*rstar)*1e-5
      r = dr
      T = Tc

      do while( p > 1e-8*Pb .and. T > 4e3)
        m2 = 0D0
        if(r >= rCentral) m2=massCentral
        call find_Hfraction( m+m2, Xin)
        call mesaeos_PTget( p, T, Xin, rho, e, s, gamma)

        if( write_out .eq. 1) then
          write(*,fmt="(6(E16.8,1X))") m/msolar, r/rsolar, rho, T, p, Xin
        endif
        y(1) = p
        y(2) = m
        y(3) = T
        call rk4step(r, r+dr, y, height, T)
        r = r + dr
        p = y(1)
        m = y(2)
        T = y(3)
        dr = 1D-2*min(min(abs(height(1)), abs(height(2))), abs(height(3)))
        mstar = m+m2
        rstar = r
      end do


      return

    end subroutine integrate_star

    double precision function spline( r, h) result( spl)
      implicit none
      double precision, intent(IN) :: r, h
      double precision :: u, u2, u3
      u = r/h
      u2 = u*u
      u3 = u2*u

      ! From Hernquist and Katz 1989
      if( u < 1D0) then
        spl = -4D0/3D0 + u2*(6D0/5D0 - 0.5D0*u)
      else if (u < 2D0) then
        spl = 1D0/15D0/u3 - 8D0/3D0 + 3d0*u - 6D0/5D0*u2 + 1D0/6D0*u3
      else
        spl = -1D0/u3
      endif

      spl = -u*spl

      return
    end function spline

    subroutine derivatives( x, y, dydx, T_guess)
      implicit none
      double precision, intent(IN) :: x, T_guess
      double precision, dimension(noeqs), intent(IN) :: y
      double precision, dimension(noeqs), intent(OUT) :: dydx
      double precision :: mass, p, r, rho, e, T, s, gamma, m2, grad, hfrac

      r = x
      p = y(1)
      mass = y(2)
      T = y(3)
      m2 = 0D0
      if(r >= rCentral) m2=massCentral
      call find_entropy( mass+m2, s, grad)
      call find_HFraction( mass+m2, hfrac)
      call mesaeos_PTget( p, T, hfrac, rho, e, s, gamma)
      dydx(1) = -G*mass/(r*r)*rho - G*massCentral/(h*h)*spline(r,h)*rho
      dydx(2) = 4.*3.141*rho*r*r
      dydx(3) = T/p*grad*dydx(1)

      return
    end subroutine derivatives

    SUBROUTINE rk4step(x1,x2,y,height, T_guess)
      IMPLICIT NONE
      double precision, DIMENSION(:), INTENT(INOUT) :: y
      double precision, INTENT(IN) :: x1, x2, T_guess
      double precision, dimension(:), intent(OUT) :: height

      INTEGER :: ndum
      double precision :: h6,hh,xh,h
      double precision, DIMENSION(size(y)) :: dym,dyt,yt
      double precision, dimension(size(y)) :: dydx

      h = x2-x1

      hh=h*0.5D0
      h6=h/6.0D0
      xh=x1+hh

      call derivatives(x1, y, dydx, T_guess)
      yt=y+hh*dydx

      call derivatives(xh,yt,dyt, T_guess)
      yt=y+hh*dyt

      call derivatives(xh,yt,dym, T_guess)
      yt=y+h*dym

      dym=dyt+dym

      call derivatives(x2,yt,dyt, T_guess)
      y=y+h6*(dydx+dyt+2D0*dym)

      height = 6D0*y/(dydx+dyt+2D0*dym)

    END SUBROUTINE rk4step

  end program hydro
