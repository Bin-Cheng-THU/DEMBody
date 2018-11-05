    !  DPASolver_beta.f90 
    !
    !  FUNCTIONS:
    !	DPASolver_beta      - Entry point of console application.
    !
    !	This version is to solve relatively long term integration.
    !
    !****************************************************************************
    !
    !  PROGRAM: OALISolver_beta Ìí¼Ó¾«È·¼ì²âÅö×²
    !
    !  PURPOSE: 
    !
    !****************************************************************************

    program OALISolver_beta

    use global

    implicit none

    call fileload
    write(*,*) "Model data loading finished!"

    call volint
    write(*,*) "Volume integration computation finished!"

    call cotrans
    write(*,*) "Coordinates transfer finished!"

    call precal
    write(*,*) "Calculating beforehand finished!"

    call normalization
    write(*,*) "Normalization finished!"

    call varsout
    write(*,*) "Varies out!"

    call volintout
    write(*,*) "Completed!"

    end program OALISolver_beta

    !****************************************************************************
    !
    !  subroutine cross: 
    !
    !  input: 
    !
    !  output: 
    !
    !****************************************************************************

    subroutine cross(a,b,c)
    use, intrinsic :: iso_fortran_env, dp=>real128
    implicit none

    real(dp) :: a(3),b(3),c(3)

    c=(/a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)/)

    return

    end subroutine cross


    !****************************************************************************
    !
    !  subroutine linspace: 
    !
    !  input: 
    !
    !  output: 
    !
    !****************************************************************************

    subroutine linspace(x1,xn,N,xlin)
    use, intrinsic :: iso_fortran_env, dp=>real128
    implicit none

    integer :: N

    real(dp) :: x1, xn, xlin(N)

    integer :: i

    real(dp) :: dx

    dx=(xn-x1)/(N-1)

    do i=1,N

    xlin(i)=x1+(i-1)*dx

    end do

    return

    end subroutine linspace