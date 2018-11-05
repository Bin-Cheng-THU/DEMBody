!****************************************************************************
!
!  subroutine varsout: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine varsout()

	use global

	implicit none

	integer :: i

	open(unit=50, file='NUM_VEF.TXT')
	write(50,'(I10)') N_v
	write(50,'(I10)') N_e
	write(50,'(I10)') N_f
	close(50)

	open(unit=51, file='V_DATA.TXT')
	open(unit=54, file='DV_DATA.TXT')
	do i=1,N_v
	   write(51,198) v_data(i,:)
	   write(54,198) dv_data(i,:)
	end do
	close(51)
	close(54)
	
	open(unit=52, file='E_DATA.TXT')
	open(unit=55, file='LE_DATA.TXT')
	open(unit=56, file='DLE_DATA.TXT')
	open(unit=57, file='NE_DATA.TXT')
	do i=1,N_e
	   write(52,199) e_data(i,:)
	   write(55,'(e40.30e3)') le_data(i)
	   write(56,'(e40.30e3)') dle_data(i)
	   write(57,202) ne_data(i,:)
	end do
	close(52)
	close(55)
	close(56)
	close(57)

	open(unit=53, file='F_DATA.TXT')
	open(unit=58, file='NF_DATA.TXT')
	do i=1,N_f
	   write(53,201) f_data(i,:)
	   write(58,198) nf_data(i,:)	   
	end do
	close(53)
	close(58)

	open(unit=59, file='NORMDATA.TXT')
    write(59,'(e40.30e3)') MXDIM
    write(59,'(e40.30e3)') kappa	   
	close(59)
	
198 format(e40.30e3,' ',e40.30e3,' ',e40.30e3)
199 format(I10,' ',I10,' ',I10,' ',I10)
201 format(I10,' ',I10,' ',I10)
202 format(e40.30e3,' ',e40.30e3,' ',e40.30e3,' ',e40.30e3,' ',e40.30e3,' ',e40.30e3)

	return

	end subroutine varsout







