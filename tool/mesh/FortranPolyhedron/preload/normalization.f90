!****************************************************************************
!
!  subroutine normalization: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine normalization()

	use global

	implicit none

	real(dp) :: tmp(6)

	tmp(1)=MAXVAL(v_data(:,1))

	tmp(2)=-MAXVAL(-v_data(:,1))

	tmp(3)=MAXVAL(v_data(:,2))

	tmp(4)=-MAXVAL(-v_data(:,2))

	tmp(5)=MAXVAL(v_data(:,3))

	tmp(6)=-MAXVAL(-v_data(:,3))

	MXDIM=MAX(tmp(1)-tmp(2),tmp(3)-tmp(4),tmp(5)-tmp(6))

	allocate(dv_data(N_v,3))
	allocate(dle_data(N_e))

	dv_data=v_data
	v_data=v_data/MXDIM

	dle_data=le_data
	le_data=le_data/MXDIM

	kappa=g*density*SPINPD**2

	return

	end subroutine normalization