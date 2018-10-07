    !********************************************************************
    !     DEMBody 4.2
    !     ***********
    !
    !     Polynomial initialization.
    !     --------------------------
    !
    !********************************************************************
    SUBROUTINE start()

    use global
    implicit none
    
    !  Initialize global scalars, counters & useful constants.
    call initialSystem

    !  Load Contact History
    call loadRestartData

    !  Read input parameters and set initial conditions: Body(I),X(K,I),Xdot(K,I); I=1,N & K=1,3.
    call loadData

    RETURN
    END