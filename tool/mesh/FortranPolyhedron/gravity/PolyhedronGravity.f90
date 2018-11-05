    !****************************************************************************
    !
    !  module PolyhedronGravity: 
    !
    !  variables:  
    !    C_density, C_GravConst
    !    gcfactor
    !    N_v, N_e, N_f
    !    v_data, ne_data, le_data, nf_data, dv_data, dle_data
    !    e_data, f_data
    !    MXDIM, kappa
    !
    !  contains:
    !    varsload():   load polyhedron gravity field basic data
    !    calgf():      calculate the gravity
    !    cross():      calculate the cross product
    !****************************************************************************

    module PolyhedronGravity
    use, intrinsic :: iso_fortran_env, dp=>real128
    implicit none

    !  constant parameters
    real(dp), parameter :: C_density=1.26d12, C_GravConst=6.674184d-20 !km, s, kg
    real(dp) :: gcfactor

    !  polyhedron gravity field basic data
    integer :: N_v, N_e, N_f
    real(dp), allocatable :: v_data(:,:), ne_data(:,:), le_data(:), nf_data(:,:)
    real(dp), allocatable :: dv_data(:,:), dle_data(:)
    integer,allocatable :: e_data(:,:), f_data(:,:)

    !  normalization data
    real(dp) :: MXDIM, kappa !km, s, kg

    contains

    !****************************************************************************
    !
    !  subroutine varsload:   load data
    !
    !  input:  none
    !
    !  output:  polyhedron gravity field basic data
    !
    !****************************************************************************
    subroutine varsload()
    implicit none

    integer :: i

    open(unit=50, file='NUM_VEF.TXT', status='old')
    read(50,*) N_v
    read(50,*) N_e
    read(50,*) N_f
    close(50)

    allocate(v_data(N_v,3))
    allocate(f_data(N_f,3))
    allocate(dv_data(N_v,3))
    allocate(dle_data(N_e))
    allocate(e_data(N_e,4))
    allocate(le_data(N_e))
    allocate(ne_data(N_e,6))
    allocate(nf_data(N_f,3))

    open(unit=51, file='V_DATA.TXT', status='old')
    open(unit=54, file='DV_DATA.TXT', status='old')
    do i=1,N_v
        read(51,*) v_data(i,:)
        read(54,*) dv_data(i,:)
    end do
    close(51)
    close(54)

    open(unit=52, file='E_DATA.TXT', status='old')
    open(unit=55, file='LE_DATA.TXT', status='old')
    open(unit=56, file='DLE_DATA.TXT', status='old')
    open(unit=57, file='NE_DATA.TXT', status='old')
    do i=1,N_e
        read(52,*) e_data(i,:)
        read(55,*) le_data(i)
        read(56,*) dle_data(i)
        read(57,*) ne_data(i,:)
    end do
    close(52)
    close(55)
    close(56)
    close(57)

    open(unit=53, file='F_DATA.TXT', status='old')
    open(unit=58, file='NF_DATA.TXT', status='old')
    do i=1,N_f
        read(53,*) f_data(i,:)
        read(58,*) nf_data(i,:)	   
    end do
    close(53)
    close(58)

    open(unit=59, file='NORMDATA.TXT', status='old')
    read(59,*) MXDIM
    read(59,*) kappa	   
    close(59)

    return
    end subroutine varsload

    !****************************************************************************
    !
    !  subroutine calgf:  calculate the gravity
    !
    !  input: rpt is field point vector, unit: m
    !
    !  output: gf is the gravitation, unit: kg, m, s
    !
    !****************************************************************************
    subroutine calgf(rpt,gf)
    implicit none

    integer :: i, j
    real(dp) :: rpt(3), fpt(3), gf(3)
    real(dp) :: lpt(N_v,4), erle(N_e,3), frof(N_f,3)
    real(dp) :: tmp, le, re(3), tmpf(3,4), fr(3,3), omegaf, tmpv(3), F1(3), F2(3)

    fpt = rpt
    fpt = fpt/1000.0D0/MXDIM
    
    ! $OMP PARALLEL DO Firstprivate(fpt) Private(i) Shared(lpt,v_data)
    do i=1,N_v
        lpt(i,1:3)=v_data(i,:)-fpt
        lpt(i,4)=sqrt(dot_product(lpt(i,1:3),lpt(i,1:3)))
    end do
    ! $OMP END PARALLEL DO
    
    ! $OMP PARALLEL DO Private(i,le,re) Shared(lpt,e_data,erle,nf_data,ne_data)
    do i=1,N_e
        tmp=sum(lpt(e_data(i,1:2),4))
        le=log((tmp+le_data(i))/(tmp-le_data(i)))
        re=lpt(e_data(i,1),1:3)
        erle(i,:)=nf_data(e_data(i,3),:)*dot_product(re,ne_data(i,1:3))&
        +nf_data(e_data(i,4),:)*dot_product(re,ne_data(i,4:6))
        erle(i,:)=erle(i,:)*le
    end do
    ! $OMP END PARALLEL DO

    do i=1,3
        F1(i)=sum(erle(:,i))
    end do    
    !do i=1,3
    !    call kahan_sum(erle(:,i),F1(i))
    !end do

    gcfactor=0.0d0

    ! $OMP PARALLEL DO Private(i,tmpf,fr,tmpv,tmp,omegaf) Shared(lpt,f_data,frof,nf_data,gcfactor)
    do i=1,N_f
        tmpf=lpt(f_data(i,:),:)
        do j=1,3
            fr(j,:)=tmpf(j,1:3)/tmpf(j,4)
        end do
        call cross(fr(2,:),fr(3,:),tmpv)
        tmp=1.0d0+dot_product(fr(1,:),fr(2,:))+dot_product(fr(2,:),fr(3,:))+dot_product(fr(3,:),fr(1,:))
        omegaf=2.0d0*atan2(dot_product(tmpv,fr(1,:)),tmp)
        frof(i,:)=dot_product(tmpf(1,1:3),nf_data(i,:))*nf_data(i,:)
        frof(i,:)=frof(i,:)*omegaf
        gcfactor=gcfactor+omegaf
    end do
    ! $OMP END PARALLEL DO

    do i=1,3
        F2(i)=sum(frof(:,i))
    end do
    !do i=1,3
    !    call kahan_sum(frof(:,i),F2(i))
    !end do

    gf=-F1+F2
    gf = gf*1000.0D0*MXDIM*C_density*C_GravConst
    
    return
    end subroutine calgf
    
    !****************************************************************************
    !
    !  subroutine cross: calculate the cross of vectors
    !
    !  input: a, b
    !
    !  output: c
    !
    !****************************************************************************
    subroutine cross(a,b,c)
    implicit none

    real(dp) :: a(3),b(3),c(3)

    c=(/a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)/)

    return
    end subroutine cross
    
    !****************************************************************************
    !
    !  subroutine kahan_sum: calculate the sum of vectors
    !
    !  input: a
    !
    !  output: b
    !
    !****************************************************************************
    subroutine kahan_sum(a,b)
    implicit none

    real(dp),intent(in),dimension(:) :: a
    real(dp) :: b
    integer :: I
    real(dp) :: sum,c,y,t
    
    sum = 0.0D0
    c = 0.0D0
    do I = 1,size(a,1)
        y = a(I) - c
        t = sum + y
        c = (t - sum) - y
        sum = t
    end do
    
    b = sum + c
    
    return
    end subroutine kahan_sum
    
    end module PolyhedronGravity