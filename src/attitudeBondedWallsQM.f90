    !********************************************************************
    !     DEMBody 4.3
    !     ***********
    !
    !     N-body integrator of Quaternion.
    !     -------------------------------
    !     @Refresh Attitude Matrix of Bonded Walls
    !
    !********************************************************************
    subroutine attitudeBondedWallsQM()
    
    use global
    implicit none

    !  Attitude Matrix in inertial frame
    bondedWallMatI(1,1) = bondedWallQ(1)**2-bondedWallQ(2)**2-bondedWallQ(3)**2+bondedWallQ(4)**2
    bondedWallMatI(2,1) = 2.0D0*(bondedWallQ(1)*bondedWallQ(2)-bondedWallQ(3)*bondedWallQ(4))
    bondedWallMatI(3,1) = 2.0D0*(bondedWallQ(1)*bondedWallQ(3)+bondedWallQ(2)*bondedWallQ(4))
    bondedWallMatI(1,2) = 2.0D0*(bondedWallQ(1)*bondedWallQ(2)+bondedWallQ(3)*bondedWallQ(4))
    bondedWallMatI(2,2) = -bondedWallQ(1)**2+bondedWallQ(2)**2-bondedWallQ(3)**2+bondedWallQ(4)**2
    bondedWallMatI(3,2) = 2.0D0*(bondedWallQ(2)*bondedWallQ(3)-bondedWallQ(1)*bondedWallQ(4))
    bondedWallMatI(1,3) = 2.0D0*(bondedWallQ(1)*bondedWallQ(3)-bondedWallQ(2)*bondedWallQ(4))
    bondedWallMatI(2,3) = 2.0D0*(bondedWallQ(2)*bondedWallQ(3)+bondedWallQ(1)*bondedWallQ(4))
    bondedWallMatI(3,3) = -bondedWallQ(1)**2-bondedWallQ(2)**2+bondedWallQ(3)**2+bondedWallQ(4)**2

    !  Attitude Matrix in body frame
    bondedWallMatB(1,1) = bondedWallQ(1)**2-bondedWallQ(2)**2-bondedWallQ(3)**2+bondedWallQ(4)**2
    bondedWallMatB(2,1) = 2.0D0*(bondedWallQ(1)*bondedWallQ(2)+bondedWallQ(3)*bondedWallQ(4))
    bondedWallMatB(3,1) = 2.0D0*(bondedWallQ(1)*bondedWallQ(3)-bondedWallQ(2)*bondedWallQ(4))
    bondedWallMatB(1,2) = 2.0D0*(bondedWallQ(1)*bondedWallQ(2)-bondedWallQ(3)*bondedWallQ(4))
    bondedWallMatB(2,2) = -bondedWallQ(1)**2+bondedWallQ(2)**2-bondedWallQ(3)**2+bondedWallQ(4)**2
    bondedWallMatB(3,2) = 2.0D0*(bondedWallQ(2)*bondedWallQ(3)+bondedWallQ(1)*bondedWallQ(4))
    bondedWallMatB(1,3) = 2.0D0*(bondedWallQ(1)*bondedWallQ(3)+bondedWallQ(2)*bondedWallQ(4))
    bondedWallMatB(2,3) = 2.0D0*(bondedWallQ(2)*bondedWallQ(3)-bondedWallQ(1)*bondedWallQ(4))
    bondedWallMatB(3,3) = -bondedWallQ(1)**2-bondedWallQ(2)**2+bondedWallQ(3)**2+bondedWallQ(4)**2
    
    end subroutine