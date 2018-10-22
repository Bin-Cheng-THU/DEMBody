    !********************************************************************
    !     DEMBody 4.4
    !     ***********
    !
    !     DEM & NBODY CODE
    !     ----------------
    !
    !     Developed by Bin Cheng, Tsinghua.
    !     .................................
    !
    !********************************************************************
    program main

    use global
    use omp_lib
    implicit none   
    
    !$OMP PARALLEL
    write(*,*) "Parallel Succeed"
    !$OMP END PARALLEL

    write(*,*) "Remember to change the Makefile !!!"
    write(*,*) "Remember to change input files  !!!"
    write(*,*) "Remember to change system files !!!"
    
    !  Read input parameters and perform initial setup.
    call start
    
    write(*,*) "< Begin to calculate, good luck!"
    
#ifdef openmp
    T1 = omp_get_wtime()
#endif

    !  Calculate total energy and produce output.
3   call output

    !  Advance solutions until next output.
4   call intgrt
    go to 3
    end
