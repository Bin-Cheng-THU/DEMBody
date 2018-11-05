!****************************************************************************
!
!  module global: 
!
!****************************************************************************

	module global

    use, intrinsic :: iso_fortran_env, dp=>real128
	implicit none

    character(20) :: filename="Bennu.bt"
	real(dp), parameter :: density=1.26d12, g=6.674184d-20, SPINPD=0d0 !m, s, kg  自转周期按s算 

	real(dp), parameter :: pi=3.1415926535897932384626433832795d0

	integer :: dispflag=1, ORBNUM=10000

	real(dp) :: relerr=1.0d-13, abserr=1.0d-13, tstart=0.0d0, tend=3877200.0d0
	real(dp) :: inicon(6)=(/0.0d0, 0.0d0, 0.0d0,&
                                    0.0d0, 0.0d0, 0.0d0/)

	real(dp) :: omgcon(3)=(/0.0d0, 0.0d0, 1.0d0/)
	real(dp) :: gcfactor
	real(dp) :: hill=2.89d4

	!polyhedron gravity field basic data

	integer :: N_v, N_e, N_f

	real(dp), allocatable :: v_data(:,:), ne_data(:,:), le_data(:), nf_data(:,:)

	real(dp), allocatable :: dv_data(:,:), dle_data(:)

	integer,allocatable :: e_data(:,:), f_data(:,:)

	!mass geometry data

	real(dp) :: T0, T1(3), T2(6), mass, rc0(3), Jmx0(3,3), Jmx(3,3)

	!normalization data
	real(dp) :: MXDIM, kappa !km, s, kg

	end module global
