!****************************************************************************
!
!  subroutine fileload: load the vertices coordinates and faces vertices data
!
!  input: fname-the asteroid polyhedron file name
!
!  output: 
!  there must not be empty row in the file
!
!****************************************************************************

	subroutine fileload()

	use global	
	
	implicit none

	integer :: i,tmp

	open(unit=5, file=filename, status='old')

	read(5,*) N_v,N_f

	allocate(v_data(N_v,3))

	allocate(f_data(N_f,3))

	do i=1,N_v

	   read(5,*) v_data(i,1),v_data(i,2),v_data(i,3)

    end do
    
    v_data = v_data / 1000.0D0

	do i=1,N_f

       read(5,*) f_data(i,1),f_data(i,2),f_data(i,3)

    end do

    f_data = f_data + 1
    
	close(5)

	return

	end subroutine fileload