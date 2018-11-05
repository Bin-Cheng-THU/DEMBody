!****************************************************************************
!
!  subroutine precal: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine precal()

	use global

	implicit none

	logical flag

	integer :: i,j,k

	integer :: tmp_data(3*N_f,2),tmp_idx(3*N_f/2,2)

	real(dp) :: lam
	
	real(dp) :: tmp(3),ra(3),rb(3),rc(3),rp(3),ap(3)

	N_e=N_f*3/2

	allocate(e_data(N_e,4))

!generate edge data e_data, format: vertex 1, vertex 2, face 1, face 2

	do i=1,N_f

	  tmp_data(3*i-2,:)=f_data(i,1:2)
	  tmp_data(3*i-1,:)=f_data(i,2:3)
	  tmp_data(3*i,:)=f_data(i,3:1:-2)

	end do

	k=0

	do i=1,N_f*3-1

	  do j=i+1,N_f*3

	    flag=((tmp_data(j,1)==tmp_data(i,1)).and.(tmp_data(j,2)==tmp_data(i,2))).or.&
		     ((tmp_data(j,1)==tmp_data(i,2)).and.(tmp_data(j,2)==tmp_data(i,1)))
		  
		if (flag) then

		  k=k+1
	      e_data(k,1:2)=tmp_data(i,:)
		  e_data(k,3)=(i-1)/3+1
		  e_data(k,4)=(j-1)/3+1

		end if

	  end do
	   
	end do

!generate edge length data le_data

	allocate(le_data(N_e))

	do i=1,N_e
	
	  tmp=v_data(e_data(i,1),:)-v_data(e_data(i,2),:)
	  le_data(i)=sqrt(dot_product(tmp,tmp))

	end do
	
!generate edge normal vector data ne_data

	allocate(ne_data(N_e,6))

	do i=1,N_e

	  rb=v_data(e_data(i,1),:)
	  rc=v_data(e_data(i,2),:)
	  
	  do j=1,2

		k=sum(f_data(e_data(i,j+2),:))-sum(e_data(i,1:2))
		ra=v_data(k,:)
		lam=dot_product(ra-rc,rb-rc)/dot_product(rb-rc,rb-rc)
		rp=lam*rb+(1-lam)*rc
		ap=rp-ra
		ne_data(i,3*j-2:3*j)=ap/sqrt(dot_product(ap,ap))

	  end do

	end do

!generate face normal vector data nf_data

	allocate(nf_data(N_f,3))

	do i=1,N_f

	  ra=v_data(f_data(i,1),:)
	  rb=v_data(f_data(i,2),:)
	  rc=v_data(f_data(i,3),:)
	  call cross(rb-ra,rc-rb,tmp)
	  nf_data(i,:)=tmp/sqrt(dot_product(tmp,tmp))
	  
	end do

!	call prevarsout

	return

	end subroutine precal

!****************************************************************************
!
!  subroutine prevarsout: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine prevarsout()

	use global 

	implicit none

	integer :: i

	open(unit=177, file="V_DATA.TXT")

	do i=1,N_v

	  write(177,"(3(e40.30e3))") v_data(i,:)

	end do

	close(177)
	
	open(unit=178, file="NE_DATA.TXT")

	do i=1,N_e

	  write(178,"(6(e40.30e3))") ne_data(i,:)

	end do

	close(178)

	open(unit=179, file="LE_DATA.TXT")

	do i=1,N_e

	  write(179,"(e40.30e3)") le_data(i)

	end do

	close(179)

	open(unit=180, file="NF_DATA.TXT")

	do i=1,N_f

	  write(180,"(3(e40.30e3))") nf_data(i,:)

	end do

	close(180)
	
	open(unit=181, file="F_DATA.TXT")

	do i=1,N_f

	  write(181,"(3(I12))") f_data(i,:)

	end do

	close(181)

	open(unit=182, file="E_DATA.TXT")

	do i=1,N_e

	  write(182,"(4(I12))") e_data(i,:)

	end do

	close(182)

	return

	end subroutine prevarsout
