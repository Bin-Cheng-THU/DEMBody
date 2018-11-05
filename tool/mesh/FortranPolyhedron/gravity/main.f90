    program main
    use PolyhedronGravity
    use omp_lib
    implicit none

    real(dp) :: fpt1(3) =(/2000.0D0,  2000.0D0,	0.0D0/)  !E4
    real(dp) :: pv1(3),torque(3)
    real(8) :: ostart,oend
    
    call varsload
    write(*,*) "Varies load finished!数据装载完毕,good luck!!"


    ostart = omp_get_wtime()
    call calgf(fpt1,pv1)
    oend = omp_get_wtime()
    write(*,*) oend - ostart
    call cross(fpt1,pv1,torque)
    write(*,*) torque
    write(*,*) pv1

    end program main