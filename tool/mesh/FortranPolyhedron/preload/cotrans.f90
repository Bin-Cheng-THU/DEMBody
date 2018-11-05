!****************************************************************************
!
!  subroutine cotrans: 
!
!  input: 
!
!  output: 
!
!****************************************************************************
INCLUDE 'link_fnl_static.h'
!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'

	subroutine cotrans()
    
	use global

	use imsl_libraries

	implicit none

	integer :: i

	real(dp) :: Jmx1(3,3), rcmx(3,1), rclen, tmp1(3,3), tmp2(3,3), tranmx(3,3), eigs(3), tmpv(3), tmpe
    real(dp) :: eyeR128(3,3) = (/1.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,1.0D0/)
    real(real64) :: Jmx1R64(3,3),eigsR64(3),tranmxR64(3,3),tmpeR64,JmxR64(3,3)
    
	do i=1,3
	  rcmx(i,1)=rc0(i)
	end do

	rclen=sqrt(dot_product(rc0,rc0))
    
	tmp1=(rclen**2)*eyeR128

	!tmp2=rcmx.xt.rcmx
    tmp2=MATMUL(rcmx,TRANSPOSE(rcmx))
    
	Jmx1=Jmx0-mass*(tmp1-tmp2)

    Jmx1R64 = Jmx1
    tranmxR64 = tranmx
	eigsR64=eig(Jmx1R64,v=tranmxR64)
    eigs = eigsR64
    Jmx1 = Jmx1R64
    tranmx = tranmxR64

	if ((eigs(1)>eigs(2)).and.(eigs(1)>eigs(3))) then

	  tmpv=tranmx(:,3)
	  tranmx(:,3)=tranmx(:,1)
	  tranmx(:,1)=tmpv

	  tmpe=eigs(3)
	  eigs(3)=eigs(1)
	  eigs(1)=tmpe

	  if (eigs(2)<eigs(1)) then

	  tmpv=tranmx(:,2)
	  tranmx(:,2)=tranmx(:,1)
	  tranmx(:,1)=tmpv

	  tmpe=eigs(2)
	  eigs(2)=eigs(1)
	  eigs(1)=tmpe

	  end if

	else if (eigs(2)>eigs(3)) then
	  
	  tmpv=tranmx(:,3)
	  tranmx(:,3)=tranmx(:,2)
	  tranmx(:,2)=tmpv

	  tmpe=eigs(3)
	  eigs(3)=eigs(2)
	  eigs(2)=tmpe

	  if (eigs(2)<eigs(1)) then

	  tmpv=tranmx(:,2)
	  tranmx(:,2)=tranmx(:,1)
	  tranmx(:,1)=tmpv

      tmpe=eigs(2)
	  eigs(2)=eigs(1)
	  eigs(1)=tmpe

	  end if

	else
	  
	  if(eigs(2)<eigs(1)) then

	  tmpv=tranmx(:,2)
	  tranmx(:,2)=tranmx(:,1)
	  tranmx(:,1)=tmpv

	  tmpe=eigs(2)
	  eigs(2)=eigs(1)
	  eigs(1)=tmpe

	  end if

    end if

    tranmxR64 = tranmx
	tmpeR64=det(tranmxR64)
    tmpe = tmpeR64

	if (tmpe<0.0) then

	  tranmx(:,1)=-tranmx(:,1)

    end if

    eigsR64 = eigs
	JmxR64=diag(eigsR64)
    Jmx = JmxR64

	do i=1,N_v

	  tmpv=v_data(i,:)
	  tmpv=tmpv-rc0
	  !v_data(i,:)=tranmx.tx.tmpv
      v_data(i,:)=MATMUL(TRANSPOSE(tranmx),tmpv)

	end do

!	call volintout

!	call vfdataout

	return

	end subroutine cotrans


!****************************************************************************
!
!  subroutine volintout: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine volintout()

	use global

	implicit none

	integer :: i

	open(unit=164, file="VOLINT_DATA.TXT")

	write(164,*) "Based on the initial fixed system, zeroth moment T0:"
	write(164,"(3x,e20.12e3)") T0
	
	write(164,*) "Based on the initial fixed system, first monent T1=(Tx,Ty,Tz):"
	write(164,"(3(3x,e20.12e3))") T1
	
	write(164,*) "Based on the initial fixed system, second monent T2=(Txx,Tyy,Tzz,Txy,Tyz,xz):"	
	write(164,"(6(3x,e20.12e3))") T2

	write(164,*) "The volume of asteroid:"
	write(164,"(3x,e20.12e3)") T0

	write(164,*) "The bulk density of asteroid:"
	write(164,"(3x,e20.12e3)") density

	write(164,*) "The mass of asteroid:"
	write(164,"(3x,e20.12e3)") mass

	write(164,*) "Based on the initial fixed system, mass center coordinates Rc0=(Xc0,yc0,Zc0):"
	write(164,"(3(3x,e20.12e3))") rc0

	write(164,*) "Based on the initial fixed system, inertial matrix J0:"

	do i=1,3

	  write(164,"(3(3x,e20.12e3))") Jmx0(i,:)

	end do

	write(164,*) "Based on the new fixed system, mass center coordinates Rc=(Xc,yc,Zc):"
	write(164,"(3(3x,e20.12e3))") 0.0d0, 0.0d0, 0.0d0

	write(164,*) "Based on the new fixed system, inertial matrix J:"

	do i=1,3

	  write(164,"(3(3x,e20.12e3))") Jmx(i,:)

	end do

	close(164)

	return

	end subroutine volintout

!****************************************************************************
!
!  subroutine vfdataout: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine vfdataout()

	use global

	implicit none

	integer :: i

	open(unit=366, file="V_DATA.TXT")

	do i=1,N_v

	  write(366,"(3(e40.30e3))") v_data(i,:)

	end do

	close(366)

	open(unit=377, file="F_DATA.TXT")

	do i=1,N_f

	  write(377,"(3(I12))") f_data(i,:)

	end do

	close(377)

	end subroutine vfdataout
 