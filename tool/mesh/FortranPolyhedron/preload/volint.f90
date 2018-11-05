!****************************************************************************
!
!  subroutine volint: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine volint()

	use global

	implicit none

	call compvolumeint

	mass=density*T0

	rc0=T1/T0

	Jmx0(1,1)=density*(T2(2)+T2(3))
	Jmx0(2,2)=density*(T2(3)+T2(1))
	Jmx0(3,3)=density*(T2(1)+T2(2))
	Jmx0(1,2)=-density*T2(4)
	Jmx0(2,1)=-density*T2(4)
	Jmx0(2,3)=-density*T2(5)
	Jmx0(3,2)=-density*T2(5)
	Jmx0(1,3)=-density*T2(6)
	Jmx0(3,1)=-density*T2(6)

	return

	end subroutine volint

!****************************************************************************
!
!  subroutine compvolumeint: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine compvolumeint()

	use global

	implicit none

	integer :: i,A,B,C

	real(dp) :: nx,ny,nz,tmp

	real(dp) :: wfs(N_f),nfs(N_f,3), F(12)!F=(Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca)

	T0=0.0d0
	T1=(/0.0d0, 0.0d0, 0.0d0/)
	T2=(/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

	call compwfnf(wfs,nfs)

	do i=1,N_f

	  nx=abs(nfs(i,1))
	  ny=abs(nfs(i,2))
	  nz=abs(nfs(i,3))
	  if ((nx>ny).and.(nx>nz)) then 
	    C=1
	  else if (ny>nz) then 
	    C=2
	  else 
	    C=3
	  end if
	  A=mod(C,3)+1
	  B=mod(A,3)+1

	  call compfaceint(wfs(i),nfs(i,:),i,A,B,C,F)

	  if (A==1) then
	    tmp=F(1)
	  else if (B==1) then
	    tmp=F(2)
	  else
	    tmp=F(3)
	  end if

	  T0=T0+nfs(i,1)*tmp

	  T1(A)=T1(A)+nfs(i,A)*F(4)
	  T1(B)=T1(B)+nfs(i,B)*F(5)
	  T1(C)=T1(C)+nfs(i,C)*F(6)
	  T2(A)=T2(A)+nfs(i,A)*F(7)
	  T2(B)=T2(B)+nfs(i,B)*F(8)
	  T2(C)=T2(C)+nfs(i,C)*F(9)
	  T2(3+A)=T2(3+A)+nfs(i,A)*F(10)
	  T2(3+B)=T2(3+B)+nfs(i,B)*F(11)
	  T2(3+C)=T2(3+C)+nfs(i,C)*F(12)

	end do

	T1=T1/2.0d0
	T2(1:3)=T2(1:3)/3.0d0
	T2(4:6)=T2(4:6)/2.0d0

	return

	end subroutine compvolumeint

!****************************************************************************
!
!  subroutine compwfnf: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine compwfnf(wfs,nfs)

	use global

	implicit none

	integer :: i

	real(dp) :: wfs(N_f), nfs(N_f,3),ra(3),rb(3),rc(3),tmp(3)

	do i=1,N_f
	  
	  ra=v_data(f_data(i,1),:)
	  rb=v_data(f_data(i,2),:)
	  rc=v_data(f_data(i,3),:)
	  call cross(rb-ra,rc-rb,tmp)
	  tmp=tmp/sqrt(dot_product(tmp,tmp))
	  wfs(i)=-dot_product(ra,tmp)
	  nfs(i,:)=tmp
	  
	end do

	return

	end subroutine compwfnf

!****************************************************************************
!
!  subroutine compfaceint: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine compfaceint(w,nf,id,A,B,C,F)
    use, intrinsic :: iso_fortran_env, dp=>real128
	implicit none

	integer :: A,B,C,id

	real(dp) :: w,nf(3),F(12),P(10) !p=(P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb)

	real(dp) :: k(4)

	call compprojectionint(A,B,id,P)

	k(1)=1.0d0/nf(C)
	k(2)=k(1)*k(1)
	k(3)=k(2)*k(1)
	k(4)=k(3)*k(1)

	F(1)=k(1)*P(2)
	F(2)=k(1)*P(3)
	F(3)=-k(2)*(nf(A)*P(2)+nf(B)*P(3)+w*P(1))

	F(4)=k(1)*P(4)
	F(5)=k(1)*P(6)
	F(6)=k(3)*((nf(A)**2)*P(4)+2.0d0*nf(A)*nf(B)*P(5)&
		 +(nf(B)**2)*P(6)+w*(2.0d0*(nf(A)*P(2)+nf(B)*P(3))+w*P(1)))

	F(7)=k(1)*P(7)
	F(8)=k(1)*P(10)
	F(9)=-k(4)*((nf(A)**3)*P(7)+3.0d0*(nf(A)**2)*nf(B)*P(8)& 
	     +3.0d0*nf(A)*(nf(B)**2)*P(9)+(nf(B)**3)*P(10)&
	     +3.0d0*w*((nf(A)**2)*P(4)+2.0d0*nf(A)*nf(B)*P(5)+(nf(B)**2)*P(6))&
	     +w*w*(3.0d0*(nf(A)*P(2)+nf(B)*P(3))+w*P(1)))

	F(10)=k(1)*P(8)
	F(11)=-k(2)*(nf(A)*P(9)+nf(B)*P(10)+w*P(6))
	F(12)=k(3)*((nf(A)**2)*P(7)+2.0d0*nf(A)*nf(B)*P(8)+(nf(B)**2)*P(9)&
	      +w*(2.0d0*(nf(A)*P(4)+nf(B)*P(5))+w*P(2)))

	return

	end subroutine compfaceint

!****************************************************************************
!
!  subroutine compprojectionint: 
!
!  input: 
!
!  output: 
!
!****************************************************************************

	subroutine compprojectionint(A,B,id,P)

	use global

	implicit none

	integer :: A,B,id,i
	real(dp) :: P(10)
	real(dp) :: a0,a1,da
	real(dp) :: b0,b1,db
	real(dp) :: a0_2,a0_3,a0_4,b0_2,b0_3,b0_4
	real(dp) :: a1_2,a1_3,b1_2,b1_3
	real(dp) :: C1,Ca,Caa,Caaa,Cb,Cbb,Cbbb
	real(dp) :: Cab,Kab,Caab,Kaab,Cabb,Kabb

	P=(/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

	do i=1,3

	  a0=v_data(f_data(id,i),A)
      b0=v_data(f_data(id,i),B)
      a1=v_data(f_data(id,mod(i,3)+1),A)
      b1=v_data(f_data(id,mod(i,3)+1),B)
      da=a1-a0
      db=b1-b0

      a0_2=a0*a0
	  a0_3=a0_2*a0
	  a0_4=a0_3*a0

      b0_2=b0*b0
	  b0_3=b0_2*b0
	  b0_4=b0_3*b0

      a1_2=a1*a1
	  a1_3=a1_2*a1 
      b1_2=b1*b1
	  b1_3=b1_2*b1

      C1=a1+a0
      Ca=a1*C1+a0_2
	  Caa=a1*Ca+a0_3
	  Caaa=a1*Caa+a0_4
      Cb=b1*(b1+b0)+b0_2
	  Cbb=b1*Cb+b0_3
	  Cbbb=b1*Cbb+b0_4
      Cab=3.0d0*a1_2+2.0d0*a1*a0+a0_2
	  Kab=a1_2+2.0d0*a1*a0+3.0d0*a0_2
      Caab=a0*Cab+4.0d0*a1_3
	  Kaab=a1*Kab+4.0d0*a0_3
      Cabb=4.0d0*b1_3+3.0d0*b1_2*b0+2.0d0*b1*b0_2+b0_3
      Kabb=b1_3+2.0d0*b1_2*b0+3.0d0*b1*b0_2+4.0d0*b0_3

!p=(P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb)

      P(1)=P(1)+db*C1
      P(2)=P(2)+db*Ca
      P(4)=P(4)+db*Caa
      P(7)=P(7)+db*Caaa
      P(3)=P(3)+da*Cb
      P(6)=P(6)+da*Cbb
      P(10)=P(10)+da*Cbbb
      P(5)=P(5)+db*(b1*Cab+b0*Kab)
      P(8)=P(8)+db*(b1*Caab+b0*Kaab)
      P(9)=P(9)+da*(a1*Cabb+a0*Kabb)

	end do

	P(1)=P(1)/2.0d0
    P(2)=P(2)/6.0d0
    P(4)=P(4)/12.0d0
    P(7)=P(7)/20.0d0
    P(3)=-P(3)/6.0d0
    P(6)=-P(6)/12.0d0
    P(10)=-P(10)/20.0d0
    P(5)=P(5)/24.0d0
    P(8)=P(8)/60.0d0
    P(9)=-P(9)/60.0d0

	return

	end subroutine compprojectionint