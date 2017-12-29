    !********************************************************************
    !     DEMBody 2.0
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
    implicit none
    
    !$OMP PARALLEL
    write(*,*) "Parallel Succeed"
    !$OMP END PARALLEL

    write(*,*) "Remember to change ParticleNumber !!!"
    write(*,*) "Remember to change the input.txt  !!!"
    write(*,*) "Remember to change the output.txt !!!"
    write(*,*) "Remember to change the CPU - time !!!"
    write(*,*) "Remember to change the Link-List  !!!"
    write(*,*) "Remember to change the Node-link  !!!"
    write(*,*) "Remember to change the Wall-Ground!!!"
    
    !  Read input parameters and perform initial setup.
    call start
    write(*,*) "Begin to calculate, good luck!"
    call CPU_TIME(T1)

    !  Calculate total energy and produce output.
3   call output

    !  Advance solutions until next output.
4   call intgrt
    go to 3
    end
