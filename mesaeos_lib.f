      module mesaeos_lib
         use eos_def
         use eos_lib
         use chem_def
         use chem_lib
         use const_lib
         use crlibm_lib
         implicit none
         integer :: handle

         integer, parameter :: species = 7
         integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7

         double precision, parameter :: Zfrac_C = 0.173312d0
         double precision, parameter :: Zfrac_N = 0.053177d0
         double precision, parameter :: Zfrac_O = 0.482398d0
         double precision, parameter :: Zfrac_Ne = 0.098675d0

      contains

      subroutine mesaeos_initialize_chem( net_iso, chem_id, xa, X, Y, Z, use_solar, fc12, fn14, fo16, fne20) 
         implicit none
         integer, pointer, dimension(:) :: net_iso, chem_id
         double precision, intent(in) :: X, Y, Z
         double precision, intent(in) :: fc12, fn14, fo16, fne20
         double precision, dimension(species), intent(out) :: xa
         integer, intent(in) :: use_solar
         net_iso(:) = 0

         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24

         xa(:) = 0.
         xa(h1) = X
         xa(he4) = Y
         if( use_solar .eq. 1) then 
            xa(c12) = Z * Zfrac_C
            xa(n14) = Z * Zfrac_N
            xa(o16) = Z * Zfrac_O
            xa(ne20) = Z * Zfrac_Ne
         else 
            xa(c12) = Z * fc12
            xa(n14) = Z * fn14
            xa(o16) = Z * fo16
            xa(ne20) = Z * fne20
         end if

         xa(species) = 1 - sum(xa(1:species-1))

      end subroutine mesaeos_initialize_chem

      subroutine mesaeos_DEget( Rho, energy, T_guess, Xin, Zin, use_solar, fc12, fn14, fo16, fne20, T, press, gamma) bind(c)
         implicit none

         double precision :: X, Y, Z, abar, zbar, z2bar, ye
         integer, pointer, dimension(:) :: net_iso, chem_id
         double precision :: xa(species)

         integer, intent(in) :: use_solar
         double precision, intent(in) :: fc12, fn14, fo16, fne20

         double precision, intent(in) :: Rho, T_guess, energy, Xin, Zin
         double precision, intent(out) :: T, press, gamma
         double precision :: Pgas,log10Rho, Prad
         double precision :: dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         double precision, dimension(num_eos_basic_results) :: res, d_dlnRho_const_T, 
     >                       d_dlnT_const_Rho, d_dabar_const_TRho, d_dzbar_const_TRho
         integer :: ierr
         double precision :: log10T

         double precision :: xz, frac, dabar_dx(species), dzbar_dx(species), sumx,
     >         mass_correction, dmc_dx(species)


         allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
         if (ierr /= 0) stop 'allocate failed'
         X = Xin
         Z = Zin
         Y = 1 - (X + Z)

         call mesaeos_initialize_chem( net_iso, chem_id, xa, X, Y, Z, use_solar, fc12, fn14, fo16, fne20)

         call composition_info(
     >         species, chem_id, xa, X, Y, xz, abar, zbar, z2bar, ye, mass_correction,
     >         sumx, dabar_dx, dzbar_dx, dmc_dx)

         T = T_guess 

         ! get a set of results for given energy and density
         call eosDE_get( 
     >         handle, Z, X, abar, zbar, 
     >         species, chem_id, net_iso, xa, 
     >         energy, log10_cr(energy), Rho, log10_cr(Rho), log10_cr(T_guess), 
     >         T, log10T, res, d_dlnRho_const_T, d_dlnT_const_Rho, 
     >         d_dabar_const_TRho, d_dzbar_const_TRho, 
     >         dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, 
     >         dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, 
     >         ierr)

         ! the indices for the results are defined in eos_def.f
         Prad = crad*T*T*T*T
         Pgas = exp_cr(res(i_lnPgas))
         press = Pgas + Prad
         gamma = res(i_gamma1)
         if( gamma < 1.) then
            write(*,*) "gamma is < 1)" , Prad, Pgas, Rho, T, gamma
         endif
         deallocate(net_iso, chem_id)

         if (ierr /= 0) then
            write(*,*) 'bad result from eos_get DE'
            write(*,*) rho, energy, T_guess, Xin, Zin
            stop 1
         end if

      end subroutine mesaeos_DEget

      subroutine mesaeos_DTget( Rho, T, Xin, Zin, use_solar, fc12, fn14, fo16, fne20, press, energy, gamma) bind(c)
         implicit none

         double precision :: X, Y, Z, abar, zbar, z2bar, ye
         integer, pointer, dimension(:) :: net_iso, chem_id
         double precision :: xa(species)

         integer, intent(in) :: use_solar
         double precision, intent(in) :: fc12, fn14, fo16, fne20


         double precision, intent(in) :: Rho, T, Xin, Zin 
         double precision, intent(out) :: energy, press, gamma
         double precision :: Pgas, log10Rho, Prad
         double precision, dimension(num_eos_basic_results) :: res,d_dlnd, d_dlnT, d_dabar, d_dzbar 
         integer :: ierr
         double precision :: log10T

         double precision :: xz, frac, dabar_dx(species), dzbar_dx(species), sumx,
     >         mass_correction, dmc_dx(species)


         allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
         if (ierr /= 0) stop 'allocate failed'
         X = Xin
         Z = Zin
         Y = 1 - (X + Z)

         call mesaeos_initialize_chem( net_iso, chem_id, xa, X, Y, Z, use_solar, fc12, fn14, fo16, fne20)

         call composition_info(
     >         species, chem_id, xa, X, Y, xz, abar, zbar, z2bar, ye, mass_correction,
     >         sumx, dabar_dx, dzbar_dx, dmc_dx)

         ! get a set of results for given energy and density
         call eosDT_get(
     >         handle, Z, X, abar, zbar, 
     >         species, chem_id, net_iso, xa,
     >         Rho, log10_cr(Rho), T, log10_cr(T), 
     >         res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         
         ! the indices for the results are defined in eos_def.f
         Prad = crad*T*T*T*T
         Pgas = exp_cr(res(i_lnPgas))
         press = Pgas + Prad
         gamma = res(i_gamma1)
         energy = exp_cr(res(i_lnE))

         deallocate(net_iso, chem_id)

         if (ierr /= 0) then
            write(*,*) 'bad result from eos_get DT'
            stop 1
         end if

      end subroutine mesaeos_DTget



      subroutine mesaeos_init( input_string) bind(c)
         use iso_c_binding, only: C_CHAR, c_null_char
         implicit none
         character (len=256) :: eos_file_prefix
         integer :: ierr
         logical, parameter :: use_cache = .true.
         character (len=256) :: my_mesa_dir

         character (kind=c_char, len=1), dimension (256), intent (in) ::
     >       input_string
         integer :: i
         type (EoS_General_Info), pointer :: rq

         my_mesa_dir = " "
         loop_string: do i=1, 256
            if ( input_string (i) == c_null_char ) then
               exit loop_string
            else
               my_mesa_dir (i:i) = input_string (i)
            end if
         end do loop_string

         ierr = 0

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

         write(*,*) 'turning off DT2'

         call eos_ptr( handle, rq, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed calling eos_ptr'
            stop
         end if
       
         rq% use_eosDT2 = .false.


      end subroutine mesaeos_init


      subroutine mesaeos_shutdown() bind(c)
         call free_eos_handle(handle)
         call eos_shutdown
      end subroutine mesaeos_shutdown
      
      end module mesaeos_lib
