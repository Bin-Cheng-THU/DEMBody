    !********************************************************************
    !     DEMBody 5.2
    !     ***********
    !
    !     Lattice generate.
    !     -------------------------------
    !     Lattice coordinates to global ID
    !
    !     
    !
    !********************************************************************
    subroutine latticeXY2ID(t_idx,t_idy,t_idz,t_ID,t_LatNx,t_LatNy)
    
    implicit none
    integer :: t_idx,t_idy,t_idz
    integer :: t_LatNx,t_LatNy
    integer :: t_ID

    t_ID = t_idx + (t_idy-1)*t_LatNx + (t_idz-1)*t_LatNx*t_LatNy

    end subroutine