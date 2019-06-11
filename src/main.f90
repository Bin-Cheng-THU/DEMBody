    !********************************************************************
    !     DEMBody 6.1
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
    integer :: thread_num
    
    write(*,*) '**************************************************'
    write(*,*) 'Copyright (c) 2016-2019 by The School of Aerospace'
    write(*,*) 'Enginnering, Tsinghua University'
    write(*,*) 
    write(*,*) 'DEMBody'
    write(*,*)
    write(*,*) 'Bin Cheng, chengb16@mails.tsinghua.edu.cn'
    write(*,*) '**************************************************'

    !$OMP PARALLEL
    thread_num = OMP_get_num_threads()
    !$OMP END PARALLEL
    
    write(*,*)
    write(*,'(A14,I2,A6)') " < Execute on ",thread_num," CPUs..."
    
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
