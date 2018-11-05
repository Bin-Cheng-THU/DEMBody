    !********************************************************************
    !     DEMBody 4.6
    !     ***********
    !
    !     Initialization of global scalars.
    !     ---------------------------------
    !
    !********************************************************************  
    SUBROUTINE initialSystem()

    use global
    use loadFile
    use omp_lib
    implicit none
    
    integer I,J,K,L
    integer :: wallFlag
    integer :: nRow
    integer :: idx,idy,idz
    integer :: idgx,idgy,idgz
    real(8) :: lx,ly,lz
    real(8) :: ERR
    integer :: numNode
    integer :: numGrid
    integer :: idNode
    real(8) :: pointTriMesh(3,3)
    real(8) :: centerLattice(3),RV(3),RVL
    integer :: TagTriMesh(3)
    integer :: idTriMesh(3,3)
    integer :: idRange(2,3)
    integer :: index
    type(trimeshLattice),pointer :: Temp
    type(trimeshLattice),pointer :: TempH
    real(8) :: ostart,oend

    !  Initialize parameters and set useful constants.
    Time = 0.0D0
    Tnext = 0.0D0
    Step = 0
    CheckPointTnext = 0.0D0
    refreshLattice = .true.
    refreshNum = 0
    
    !  Initialize control parameters
    write(*,*) "< System file loading..."    
    open (1000,FILE="../Input/systemControl.dembody",STATUS='OLD',BLANK='NULL',POSITION='REWIND')

    read (1000,*) vsDEMBody,vsDEMBody  !  DEM version
    if (vsDEMBody .NE. VERSION) then
        write(*,*) "DEMBody Version Mismatching!"
        stop
    end if
    read (1000,*) 
    read (1000,*) LatDx,LatDy,LatDz    !  Parallel Lattice grid interval
    read (1000,*) LatNx,LatNy,LatNz    !  Parallel Lattice grid number
    read (1000,*) LatMx,LatMy,LatMz    !  Parallel Lattice grid origin
    read (1000,*) verlet               !  verlet distance
    read (1000,*)
    read (1000,*) Deltat,Tcrit         !  output time interval and termination time
    read (1000,*) CheckPointDt         !  check point time interval
    read (1000,*) Dt                   !  time step of second-step integral
    read (1000,*)                                            
    read (1000,*) isPlanet             !  whether use Planet Gravity function
    read (1000,*) isRotSystem          !  whether use Rotary System function 
    read (1000,*) isQuaternion         !  whether intergrate Quaternion
    read (1000,*) isContactWall        !  whether use Contactable Walls
    read (1000,*) isMovingWall         !  whether use Moving Walls
    read (1000,*) isBondedWall         !  whether use Bonded Walls
    read (1000,*) isTriMeshWall        !  whether use TriMesh Walls
    read (1000,*) isBondedTriMeshWall  !  whether use Bonded TriMesh Walls
    read (1000,*) isFunnelWall         !  whether use Funnel Walls
    read (1000,*) isPeriodic           !  whether use Periodic function
    read (1000,*) isGravBody           !  whether use Gravity Body
    read (1000,*) isGravTriMesh        !  whether use Gravity TriMesh
    read (1000,*) MAX_ACC              !  maximum contact acceleration
    read (1000,*)                         
    read (1000,*) G(1),G(2),G(3)       !  the totle gravity    
    
    wallFlag = 1
    
    !  initial the contactable walls
    if (isContactWall) then
        write(*,*) '< is Contactable walls, loading...'
        read (1000,*) 
        read (1000,*) contactWallNum
        allocate (contactWallTag(contactWallNum))
        allocate (contactWallPoint(3,contactWallNum))
        allocate (contactWallVector(3,contactWallNum))
        do I = 1,contactWallNum
            contactWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (contactWallPoint(K,I),K=1,3),(contactWallVector(K,I),K=1,3)
        end do
    else
        read (1000,*) 
        read (1000,*)
    end if
    
    !  initial the bonded walls
    if (isBondedWall) then
        write(*,*) '< is Bonded walls, loading...'
        read (1000,*)
        read (1000,*) (bondedWallX(K),K=1,3),(bondedWallXdot(K),K=1,3),(bondedWallW(K),K=1,3)
        read (1000,*) (bondedWallQ(K),K=1,4)
        read (1000,*) bondedWallBody,(bondedWallInertia(K),K=1,3)
        read (1000,*) bondedWallNum
        allocate (bondedWallTag(bondedWallNum))
        allocate (bondedWallPoint(3,bondedWallNum))
        allocate (bondedWallVectorN(3,bondedWallNum))
        allocate (bondedWallVectorTx(3,bondedWallNum))
        allocate (bondedWallVectorTy(3,bondedWallNum))
        allocate (bondedWallLx(bondedWallNum))
        allocate (bondedWallLy(bondedWallNum))
        do I = 1,bondedWallNum
            bondedWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (bondedWallPoint(K,I),K=1,3),(bondedWallVectorN(K,I),K=1,3),(bondedWallVectorTx(K,I),K=1,3),(bondedWallVectorTy(K,I),K=1,3),bondedWallLx(I),bondedWallLy(I)
        end do
        call attitudeQ2M(bondedWallQ,bondedWallMatI,bondedWallMatB)
        !  unified input parameters
        do K = 1,3
            bondedWallWB(K) = bondedWallMatI(K,1)*bondedWallW(1) + bondedWallMatI(K,2)*bondedWallW(2) + bondedWallMatI(K,3)*bondedWallW(3)
        end do
        !  load bonded wall mesh
        open (1500,File="../Input/bondedWallPoint.vtk")
        nRow = GetFileN(1500)
        allocate (bondedWallMeshPoint(3,nRow))
        do I = 1,nRow
            read (1500,*) (bondedWallMeshPoint(K,I),K=1,3)
        end do
        close(1500)
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the trimesh walls
    if (isTriMeshWall) then
        write(*,*) '< is TriMesh walls, loading...'
        read (1000,*)
        read (1000,*) trimeshWallNum
        read (1000,*)
        allocate (trimeshWallTag(trimeshWallNum))
        allocate (trimeshWallPoint(3,trimeshWallNum))
        allocate (trimeshWallVectorN(3,trimeshWallNum))
        allocate (trimeshWallVectorTx(3,trimeshWallNum))
        allocate (trimeshWallVectorTy(3,trimeshWallNum))
        allocate (trimeshWallLength(4,trimeshWallNum))
        !  load trimesh wall
        open (1500,File="../Input/trimeshWall.mesh")
        do I = 1,trimeshWallNum
            trimeshWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1500,*) (trimeshWallPoint(K,I),K=1,3),(trimeshWallVectorTx(K,I),K=1,3),(trimeshWallVectorTy(K,I),K=1,3)
            !  initial normal vector
            trimeshWallVectorN(1,I) = trimeshWallVectorTx(2,I)*trimeshWallVectorTy(3,I)-trimeshWallVectorTx(3,I)*trimeshWallVectorTy(2,I)
            trimeshWallVectorN(2,I) = trimeshWallVectorTx(3,I)*trimeshWallVectorTy(1,I)-trimeshWallVectorTx(1,I)*trimeshWallVectorTy(3,I)
            trimeshWallVectorN(3,I) = trimeshWallVectorTx(1,I)*trimeshWallVectorTy(2,I)-trimeshWallVectorTx(2,I)*trimeshWallVectorTy(1,I)
            ERR = trimeshWallVectorN(1,I)*trimeshWallVectorN(1,I) + trimeshWallVectorN(2,I)*trimeshWallVectorN(2,I) + trimeshWallVectorN(3,I)*trimeshWallVectorN(3,I)
            do K = 1,3
                trimeshWallVectorN(K,I) = trimeshWallVectorN(K,I)/sqrt(ERR)
            end do
            !  initial length vector
            trimeshWallLength(1,I) = trimeshWallVectorTx(1,I)*trimeshWallVectorTx(1,I) + trimeshWallVectorTx(2,I)*trimeshWallVectorTx(2,I) + trimeshWallVectorTx(3,I)*trimeshWallVectorTx(3,I)
            trimeshWallLength(2,I) = trimeshWallVectorTx(1,I)*trimeshWallVectorTy(1,I) + trimeshWallVectorTx(2,I)*trimeshWallVectorTy(2,I) + trimeshWallVectorTx(3,I)*trimeshWallVectorTy(3,I)
            trimeshWallLength(3,I) = trimeshWallVectorTy(1,I)*trimeshWallVectorTy(1,I) + trimeshWallVectorTy(2,I)*trimeshWallVectorTy(2,I) + trimeshWallVectorTy(3,I)*trimeshWallVectorTy(3,I)
            trimeshWallLength(4,I) = trimeshWallLength(2,I)*trimeshWallLength(2,I) - trimeshWallLength(1,I)*trimeshWallLength(3,I)
        end do
        close(1500)
        !!********************************Test for trimeshLattice******************************************
        !ostart = omp_get_wtime()
        !open(3000,File="trimeshLattice_test.csv")
        !write(3000,*) 'X',',','Y',',','Z',',','N'
        !  allocate trimesh to Lattice
        LatNum = LatNx*LatNy*LatNz 
        allocate(trimeshDEM(LatNum))
        do I = 1,LatNum
            trimeshDEM(I)%No = 0
            nullify(trimeshDEM(I)%next)
        end do
        do I = 1,trimeshWallNum
            !  triangle vertexes
            pointTriMesh(1,1) = trimeshWallPoint(1,I)
            pointTriMesh(2,1) = trimeshWallPoint(2,I)
            pointTriMesh(3,1) = trimeshWallPoint(3,I)
            pointTriMesh(1,2) = trimeshWallPoint(1,I) + trimeshWallVectorTx(1,I)
            pointTriMesh(2,2) = trimeshWallPoint(2,I) + trimeshWallVectorTx(2,I)
            pointTriMesh(3,2) = trimeshWallPoint(3,I) + trimeshWallVectorTx(3,I)
            pointTriMesh(1,3) = trimeshWallPoint(1,I) + trimeshWallVectorTy(1,I)
            pointTriMesh(2,3) = trimeshWallPoint(2,I) + trimeshWallVectorTy(2,I)
            pointTriMesh(3,3) = trimeshWallPoint(3,I) + trimeshWallVectorTy(3,I)
            !  transform points to Inner Mesh
            do J = 1,3
                if (pointTriMesh(1,J) .LT. -LatMx) then
                    pointTriMesh(1,J) = -LatMx
                else if (pointTriMesh(1,J) .GT. -LatMx+(LatDx*LatNx)) then
                    pointTriMesh(1,J) = -LatMx+(LatDx*LatNx)
                end if
                if (pointTriMesh(2,J) .LT. -LatMy) then
                    pointTriMesh(2,J) = -LatMy
                else if (pointTriMesh(2,J) .GT. -LatMy+(LatDy*LatNy)) then
                    pointTriMesh(2,J) = -LatMy+(LatDy*LatNy)
                end if
                if (pointTriMesh(3,J) .LT. -LatMz) then
                    pointTriMesh(3,J) = -LatMz
                else if (pointTriMesh(3,J) .GT. -LatMz+(LatDz*LatNz)) then
                    pointTriMesh(3,J) = -LatMz+(LatDz*LatNz)
                end if
            end do
            !  find id range
            do J = 1,3
                TagTriMesh(J) = floor((pointTriMesh(1,J)+LatMx)/LatDx) + 1 + floor((pointTriMesh(2,J)+LatMy)/LatDy)*LatNx + floor((pointTriMesh(3,J)+LatMz)/LatDz)*(LatNx*LatNy)
            end do
            do J = 1,3
                idTriMesh(3,J) = int((TagTriMesh(J)-1)/(LatNx*LatNy))+1
                idTriMesh(2,J) = int((TagTriMesh(J)-(idTriMesh(3,J)-1)*(LatNx*LatNy)-1)/LatNx)+1
                idTriMesh(1,J) = TagTriMesh(J)-(idTriMesh(3,J)-1)*(LatNx*LatNy)-(idTriMesh(2,J)-1)*LatNx
            end do
            idRange(1,1) = min(idTriMesh(1,1),idTriMesh(1,2),idTriMesh(1,3))
            idRange(1,1) = max(idRange(1,1)-1,1)
            idRange(2,1) = max(idTriMesh(1,1),idTriMesh(1,2),idTriMesh(1,3))
            idRange(2,1) = min(idRange(2,1)+1,LatNx)
            idRange(1,2) = min(idTriMesh(2,1),idTriMesh(2,2),idTriMesh(2,3))
            idRange(1,2) = max(idRange(1,2)-1,1)
            idRange(2,2) = max(idTriMesh(2,1),idTriMesh(2,2),idTriMesh(2,3))
            idRange(2,2) = min(idRange(2,2)+1,LatNy)
            idRange(1,3) = min(idTriMesh(3,1),idTriMesh(3,2),idTriMesh(3,3))
            idRange(1,3) = max(idRange(1,3)-1,1)
            idRange(2,3) = max(idTriMesh(3,1),idTriMesh(3,2),idTriMesh(3,3))
            idRange(2,3) = min(idRange(2,3)+1,LatNx)
            do J = 1,(idRange(2,1)-idRange(1,1)+1)
                do K = 1,(idRange(2,2)-idRange(1,2)+1)
                    do L = 1,(idRange(2,3)-idRange(1,3)+1)
                        centerLattice(1) = dble(idRange(1,1)+J-1-1+0.5)*LatDx - LatMx
                        centerLattice(2) = dble(idRange(1,2)+K-1-1+0.5)*LatDy - LatMy
                        centerLattice(3) = dble(idRange(1,3)+L-1-1+0.5)*LatDz - LatMz
                        RV(1) = centerLattice(1) - trimeshWallPoint(1,I)
                        RV(2) = centerLattice(2) - trimeshWallPoint(2,I)
                        RV(3) = centerLattice(3) - trimeshWallPoint(3,I)
                        RVL = abs(RV(1)*trimeshWallVectorN(1,I) + RV(2)*trimeshWallVectorN(2,I) + RV(3)*trimeshWallVectorN(3,I))
                        if (RVL .LE. 2.0*LatDx) then
                            index = (idRange(1,1)+J-1) + (idRange(1,2)+K-1-1)*LatNx + (idRange(1,3)+L-1-1)*LatNx*LatNy
                            trimeshDEM(index)%No = trimeshDEM(index)%No + 1
                            Temp => trimeshDEM(index)
                            do while (associated(Temp%next))
                                Temp => Temp%next
                            end do
                            allocate(TempH)
                            TempH = trimeshLattice(I,NULL())
                            Temp%next => TempH
                            !write(3000,"(F10.5,A2,F10.5,A2,F10.5,A2,I6)") (idRange(1,1)+J-1-1+0.5)*LatDx-LatMx,',',(idRange(1,2)+K-1-1+0.5)*LatDy-LatMy,',',(idRange(1,3)+L-1-1+0.5)*LatDz-LatMz,',',I
                        end if
                    end do
                end do     
            end do
        end do
        !open(3000,File='trimeshLattice_num.txt')
        !do I = 1,LatNum
        !    write(3000,*) trimeshDEM(I)%No
        !end do
        !close(3000)
        !oend = omp_get_wtime()
        !write(*,*) oend-ostart
    else
        read (1000,*)
        read (1000,*)
    end if

    !  initial the bonded trimesh walls
    if (isBondedTriMeshWall) then
        write(*,*) '< is Bonded TriMesh walls, loading...'
        read (1000,*)
        read (1000,*) (bondedTriMeshWallX(K),K=1,3),(bondedTriMeshWallXdot(K),K=1,3),(bondedTriMeshWallW(K),K=1,3)
        read (1000,*) (bondedTriMeshWallQ(K),K=1,4)
        read (1000,*) bondedTriMeshWallBody,(bondedTriMeshWallInertia(K),K=1,3)
        read (1000,*) bondedTriMeshWallNum
        read (1000,*)
        allocate (bondedTriMeshWallTag(bondedTriMeshWallNum))
        allocate (bondedTriMeshWallPoint(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallVectorN(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallVectorTx(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallVectorTy(3,bondedTriMeshWallNum))
        allocate (bondedTriMeshWallLength(4,bondedTriMeshWallNum))
        !  load trimesh wall
        open (1500,File="../Input/bondedTriMeshWall.mesh")
        do I = 1,bondedTriMeshWallNum
            bondedTriMeshWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1500,*) (bondedTriMeshWallPoint(K,I),K=1,3),(bondedTriMeshWallVectorTx(K,I),K=1,3),(bondedTriMeshWallVectorTy(K,I),K=1,3)
            !  initial normal vector
            bondedTriMeshWallVectorN(1,I) = bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTy(3,I)-bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTy(2,I)
            bondedTriMeshWallVectorN(2,I) = bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTy(1,I)-bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTy(3,I)
            bondedTriMeshWallVectorN(3,I) = bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTy(2,I)-bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTy(1,I)
            ERR = bondedTriMeshWallVectorN(1,I)*bondedTriMeshWallVectorN(1,I) + bondedTriMeshWallVectorN(2,I)*bondedTriMeshWallVectorN(2,I) + bondedTriMeshWallVectorN(3,I)*bondedTriMeshWallVectorN(3,I)
            do K = 1,3
                bondedTriMeshWallVectorN(K,I) = bondedTriMeshWallVectorN(K,I)/sqrt(ERR)
            end do
            !  initial length vector
            bondedTriMeshWallLength(1,I) = bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTx(1,I) + bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTx(2,I) + bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTx(3,I)
            bondedTriMeshWallLength(2,I) = bondedTriMeshWallVectorTx(1,I)*bondedTriMeshWallVectorTy(1,I) + bondedTriMeshWallVectorTx(2,I)*bondedTriMeshWallVectorTy(2,I) + bondedTriMeshWallVectorTx(3,I)*bondedTriMeshWallVectorTy(3,I)
            bondedTriMeshWallLength(3,I) = bondedTriMeshWallVectorTy(1,I)*bondedTriMeshWallVectorTy(1,I) + bondedTriMeshWallVectorTy(2,I)*bondedTriMeshWallVectorTy(2,I) + bondedTriMeshWallVectorTy(3,I)*bondedTriMeshWallVectorTy(3,I)
            bondedTriMeshWallLength(4,I) = bondedTriMeshWallLength(2,I)*bondedTriMeshWallLength(2,I) - bondedTriMeshWallLength(1,I)*bondedTriMeshWallLength(3,I)
        end do
        close(1500)
        !!********************************Test for bondedTriMeshLattice******************************************
        !  allocate trimesh to Lattice
        LatNum = LatNx*LatNy*LatNz 
        allocate(bondedTriMeshDEM(LatNum))
        do I = 1,LatNum
            bondedTriMeshDEM(I)%No = 0
            nullify(bondedTriMeshDEM(I)%next)
        end do
        do I = 1,bondedTriMeshWallNum
            !  triangle vertexes
            pointTriMesh(1,1) = bondedTriMeshWallPoint(1,I)
            pointTriMesh(2,1) = bondedTriMeshWallPoint(2,I)
            pointTriMesh(3,1) = bondedTriMeshWallPoint(3,I)
            pointTriMesh(1,2) = bondedTriMeshWallPoint(1,I) + bondedTriMeshWallVectorTx(1,I)
            pointTriMesh(2,2) = bondedTriMeshWallPoint(2,I) + bondedTriMeshWallVectorTx(2,I)
            pointTriMesh(3,2) = bondedTriMeshWallPoint(3,I) + bondedTriMeshWallVectorTx(3,I)
            pointTriMesh(1,3) = bondedTriMeshWallPoint(1,I) + bondedTriMeshWallVectorTy(1,I)
            pointTriMesh(2,3) = bondedTriMeshWallPoint(2,I) + bondedTriMeshWallVectorTy(2,I)
            pointTriMesh(3,3) = bondedTriMeshWallPoint(3,I) + bondedTriMeshWallVectorTy(3,I)
            !  transform points to Inner Mesh
            do J = 1,3
                if (pointTriMesh(1,J) .LT. -LatMx) then
                    pointTriMesh(1,J) = -LatMx
                else if (pointTriMesh(1,J) .GT. -LatMx+(LatDx*LatNx)) then
                    pointTriMesh(1,J) = -LatMx+(LatDx*LatNx)
                end if
                if (pointTriMesh(2,J) .LT. -LatMy) then
                    pointTriMesh(2,J) = -LatMy
                else if (pointTriMesh(2,J) .GT. -LatMy+(LatDy*LatNy)) then
                    pointTriMesh(2,J) = -LatMy+(LatDy*LatNy)
                end if
                if (pointTriMesh(3,J) .LT. -LatMz) then
                    pointTriMesh(3,J) = -LatMz
                else if (pointTriMesh(3,J) .GT. -LatMz+(LatDz*LatNz)) then
                    pointTriMesh(3,J) = -LatMz+(LatDz*LatNz)
                end if
            end do
            !  find id range
            do J = 1,3
                TagTriMesh(J) = floor((pointTriMesh(1,J)+LatMx)/LatDx) + 1 + floor((pointTriMesh(2,J)+LatMy)/LatDy)*LatNx + floor((pointTriMesh(3,J)+LatMz)/LatDz)*(LatNx*LatNy)
            end do
            do J = 1,3
                idTriMesh(3,J) = int((TagTriMesh(J)-1)/(LatNx*LatNy))+1
                idTriMesh(2,J) = int((TagTriMesh(J)-(idTriMesh(3,J)-1)*(LatNx*LatNy)-1)/LatNx)+1
                idTriMesh(1,J) = TagTriMesh(J)-(idTriMesh(3,J)-1)*(LatNx*LatNy)-(idTriMesh(2,J)-1)*LatNx
            end do
            idRange(1,1) = min(idTriMesh(1,1),idTriMesh(1,2),idTriMesh(1,3))
            idRange(1,1) = max(idRange(1,1)-1,1)
            idRange(2,1) = max(idTriMesh(1,1),idTriMesh(1,2),idTriMesh(1,3))
            idRange(2,1) = min(idRange(2,1)+1,LatNx)
            idRange(1,2) = min(idTriMesh(2,1),idTriMesh(2,2),idTriMesh(2,3))
            idRange(1,2) = max(idRange(1,2)-1,1)
            idRange(2,2) = max(idTriMesh(2,1),idTriMesh(2,2),idTriMesh(2,3))
            idRange(2,2) = min(idRange(2,2)+1,LatNy)
            idRange(1,3) = min(idTriMesh(3,1),idTriMesh(3,2),idTriMesh(3,3))
            idRange(1,3) = max(idRange(1,3)-1,1)
            idRange(2,3) = max(idTriMesh(3,1),idTriMesh(3,2),idTriMesh(3,3))
            idRange(2,3) = min(idRange(2,3)+1,LatNx)
            do J = 1,(idRange(2,1)-idRange(1,1)+1)
                do K = 1,(idRange(2,2)-idRange(1,2)+1)
                    do L = 1,(idRange(2,3)-idRange(1,3)+1)
                        centerLattice(1) = dble(idRange(1,1)+J-1-1+0.5)*LatDx - LatMx
                        centerLattice(2) = dble(idRange(1,2)+K-1-1+0.5)*LatDy - LatMy
                        centerLattice(3) = dble(idRange(1,3)+L-1-1+0.5)*LatDz - LatMz
                        RV(1) = centerLattice(1) - bondedTriMeshWallPoint(1,I)
                        RV(2) = centerLattice(2) - bondedTriMeshWallPoint(2,I)
                        RV(3) = centerLattice(3) - bondedTriMeshWallPoint(3,I)
                        RVL = abs(RV(1)*bondedTriMeshWallVectorN(1,I) + RV(2)*bondedTriMeshWallVectorN(2,I) + RV(3)*bondedTriMeshWallVectorN(3,I))
                        if (RVL .LE. 2.0*LatDx) then
                            index = (idRange(1,1)+J-1) + (idRange(1,2)+K-1-1)*LatNx + (idRange(1,3)+L-1-1)*LatNx*LatNy
                            bondedTriMeshDEM(index)%No = bondedTriMeshDEM(index)%No + 1
                            Temp => bondedTriMeshDEM(index)
                            do while (associated(Temp%next))
                                Temp => Temp%next
                            end do
                            allocate(TempH)
                            TempH = trimeshLattice(I,NULL())
                            Temp%next => TempH
                            !write(3000,"(F10.5,A2,F10.5,A2,F10.5,A2,I6)") (idRange(1,1)+J-1-1+0.5)*LatDx-LatMx,',',(idRange(1,2)+K-1-1+0.5)*LatDy-LatMy,',',(idRange(1,3)+L-1-1+0.5)*LatDz-LatMz,',',I
                        end if
                    end do
                end do     
            end do
        end do
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the funnel walls
    if (isFunnelWall) then
        write(*,*) '< is Funnel walls, loading...'
        read (1000,*) 
        read (1000,*) funnelWallNum
        allocate (funnelWallTag(funnelWallNum))
        allocate (funnelWallPoint(3,funnelWallNum))
        allocate (funnelWallVector(3,funnelWallNum))
        allocate (funnelWallRadius(2,funnelWallNum))
        allocate (funnelWallLength(funnelWallNum))
        do I = 1,funnelWallNum
            funnelWallTag(I) = NMAX + wallFlag
            wallFlag = wallFlag + 1
            read (1000,*) (funnelWallPoint(K,I),K=1,3),(funnelWallVector(K,I),K=1,3),(funnelWallRadius(K,I),K=1,2),funnelWallLength(I)
        end do
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the geometry boundary of the Particle Box
    if (isPeriodic) then
        write(*,*) '< is Periodic boundary, loading...'
        read (1000,*)
        read (1000,*) (PlaSx1p(K),K=1,3),(PlaSx1v(K),K=1,3)
        read (1000,*) (PlaSx2p(K),K=1,3),(PlaSx2v(K),K=1,3)
        read (1000,*) (PlaSy1p(K),K=1,3),(PlaSy1v(K),K=1,3)
        read (1000,*) (PlaSy2p(K),K=1,3),(PlaSy2v(K),K=1,3)
        read (1000,*) LenBoxX
        read (1000,*) LenBoxY
        read (1000,*) gamma
        LenBoxX = PlaSx1p(1) - PlaSx2p(1)
        LenBoxY = PlaSy1p(2) - PlaSy2p(2)
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the gravity body
    if (isGravBody) then
        write(*,*) '< is Gravity Body, loading...'
        read (1000,*)
        gravBodyTag = NMAX + wallFlag
        read (1000,*) (gravBodyX(K),K=1,3),(gravBodyXdot(K),K=1,3),(gravBodyW(K),K=1,3)
        read (1000,*) gravBodyBody,gravBodyR,gravBodyInertia
        gravBodyQ(1) = 0.0D0
        gravBodyQ(2) = 0.0D0
        gravBodyQ(3) = 0.0D0
        gravBodyQ(4) = 1.0D0
    else
        read (1000,*)
        read (1000,*)
    end if
    
    !  initial the properties of Saturn and Pan
    if (isPlanet) then
        write(*,*) '< is Planet system, loading'
        read (1000,*)
        read (1000,*) muS, muP            !  Gravity constant of Saturn and Pan
        read (1000,*) omega               !  Angular velocity of Pan's initial revolution (suppose tidal locking and anticlockwise, usually negative)
        read (1000,*) (rOrig(K),K=1,3)    !  Saturn to Pan's initial position, i.e., the origin point
        omega = int(omega/abs(omega))*sqrt(abs(muS/rOrig(2)/rOrig(2)/rOrig(2)))
    else
        read (1000,*)
        read (1000,*)
    end if
        
    !  initial the rotary system
    if (isRotSystem) then
        write(*,*) '< is Rotary System, loading...'
        read (1000,*)
        read (1000,*) sysOmega, sysGrav
    else
        read (1000,*)
        read (1000,*)
    end if

    !  initial the grav trimesh
    if (isGravTriMesh) then
        write(*,*) '< is Grav TriMesh, loading...'
        read (1000,*)
        read (1000,*) TriLatDx,TriLatDy,TriLatDz
        read (1000,*) TriLatMx,TriLatMy,TriLatMz
        read (1000,*) TriLatNx,TriLatNy,TriLatNz
        read (1000,*)
        open (2000,FILE="../Input/gravTriMesh.force",STATUS='OLD',BLANK='NULL',POSITION='REWIND')
        numNode = (TriLatNx+1)*(TriLatNy+1)*(TriLatNz+1)
        numGrid = TriLatNx*TriLatNy*TriLatNz
        allocate(GravTriMeshNode(3,numNode))
        do I = 1,numNode
            read (2000,*) (GravTriMeshNode(K,I),K=1,3)
        end do
        close(2000)
        
        allocate(GravTriMeshGrid(numGrid))
        do I = 1,numGrid
            idz = int((I-1)/(TriLatNx*TriLatNy))+1
            idy = int((I-(idz-1)*(TriLatNx*TriLatNy)-1)/TriLatNx)+1
            idx = I-(idz-1)*(TriLatNx*TriLatNy)-(idy-1)*TriLatNx
            lx = dble(idx-1)*TriLatDx - TriLatMx
            ly = dble(idy-1)*TriLatDy - TriLatMy
            lz = dble(idz-1)*TriLatDz - TriLatMz            
            GravTriMeshGrid(I)%PositionCenter(1) = lx + 0.5D0*TriLatDx
            GravTriMeshGrid(I)%PositionCenter(2) = ly + 0.5D0*TriLatDy
            GravTriMeshGrid(I)%PositionCenter(3) = lz + 0.5D0*TriLatDz
            idNode = idx + (idy-1)*(TriLatNx+1) + (idz-1)*(TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(1) = idNode
            GravTriMeshGrid(I)%NodeID(2) = idNode + 1
            GravTriMeshGrid(I)%NodeID(3) = idNode + 2 + TriLatNx
            GravTriMeshGrid(I)%NodeID(4) = idNode + 1 + TriLatNx
            GravTriMeshGrid(I)%NodeID(5) = idNode + (TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(6) = idNode + 1 + (TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(7) = idNode + 2 + TriLatNx + (TriLatNx+1)*(TriLatNy+1)
            GravTriMeshGrid(I)%NodeID(8) = idNode + 1 + TriLatNx + (TriLatNx+1)*(TriLatNy+1)
        end do
    else
        read (1000,*)
        read (1000,*)
    end if
    close(1000)

#ifdef LatticeSearch
    !  Initialize Parallel Lattice.
    write(*,*) "< Parallel Lattice initializing..."
    
    LatNum = LatNx*LatNy*LatNz             !  Paraller Number
    allocate(DEM(LatNum))
    do I = 1,LatNum
        idz = int((I-1)/(LatNx*LatNy))+1
        idy = int((I-(idz-1)*(LatNx*LatNy)-1)/LatNx)+1
        idx = I-(idz-1)*(LatNx*LatNy)-(idy-1)*LatNx
        lx = dble(idx-1)*LatDx - LatMx
        ly = dble(idy-1)*LatDy - LatMy
        lz = dble(idz-1)*LatDz - LatMz
        !DEM(I)%ID(1) = idx                 !  IDx of Lattice
        !DEM(I)%ID(2) = idy                 !  IDy of Lattice
        !DEM(I)%ID(3) = idz                 !  IDz of Lattice
#ifdef nConfined        
        DEM(I)%PositionD(1) = lx           !  PositionDx of Lattice
        DEM(I)%PositionD(2) = ly           !  PositionDy of Lattice
        DEM(I)%PositionD(3) = lz           !  PositionDz of Lattice
        DEM(I)%PositionU(1) = lx + LatDx   !  PositionUx of Lattice
        DEM(I)%PositionU(2) = ly + LatDy   !  PositionUy of Lattice
        DEM(I)%PositionU(3) = lz + LatDz   !  PositionUz of Lattice
#endif      
        DEM(I)%NeighborID = 0
        !  Left lattice
        if (idx .NE. 1) then
            DEM(I)%NeighborID(1) = I - 1
        end if
        !  Left Upper lattice
        if (idx.NE.1 .AND. idy.NE.LatNy) then
            DEM(I)%NeighborID(2) = I - 1 + LatNx
        end if
        !  Upper lattice
        if (idy .NE. LatNy) then
            DEM(I)%NeighborID(3) = I + LatNx
        end if
        !  Right Upper lattice
        if (idx.NE.LatNx .AND. idy.NE.LatNy) then
            DEM(I)%NeighborID(4) = I + 1 + LatNx
        end if
        !  Center Covered lattice
        if (idz .NE. LatNz) then
            DEM(I)%NeighborID(5) = I + LatNx*LatNy
        end if
        !  Left Covered lattice
        if (idx.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(6) = I - 1 + LatNx*LatNy
        end if
        !  Left Upper Covered lattice
        if (idx.NE.1 .AND. idy.NE.LatNy .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(7) = I - 1 + LatNx + LatNx*LatNy
        end if
        !  Upper Covered lattice
        if (idy.NE.LatNy .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(8) = I + LatNx + LatNx*LatNy
        end if
        !  Right Upper Covered lattice
        if (idx.NE.LatNx .AND. idy.NE.LatNy .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(9) = I + 1 + LatNx + LatNx*LatNy
        end if
        !  Right Covered lattice
        if (idx.NE.LatNx .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(10) = I + 1 + LatNx*LatNy
        end if
        !  Right Below Covered lattice
        if (idx.NE.LatNx .AND. idy.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(11) = I + 1 - LatNx + LatNx*LatNy
        end if
        !  Below Covered lattice
        if (idy.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(12) = I - LatNx + LatNx*LatNy
        end if
        !  Left Below Covered lattice
        if (idx.NE.1 .AND. idy.NE.1 .AND. idz.NE.LatNz) then
            DEM(I)%NeighborID(13) = I - 1 - LatNx + LatNx*LatNy
        end if
        !  Right lattice
        if (idx .NE. LatNx) then
            DEM(I)%NeighborID(14) = I + 1
        end if
        !  Right Below lattice
        if (idx.NE.LatNx .AND. idy.NE.1) then
            DEM(I)%NeighborID(15) = I + 1 - LatNx
        end if
        !  Below lattice
        if (idy .NE. 1) then
            DEM(I)%NeighborID(16) = I - LatNx
        end if
        !  Left Below lattice
        if (idx.NE.1 .AND. idy.NE.1) then
            DEM(I)%NeighborID(17) = I - 1 - LatNx
        end if
        !  Center Under lattice
        if (idz .NE. 1) then
            DEM(I)%NeighborID(18) = I - LatNx*LatNy
        end if
        !  Right Under lattice
        if (idx.NE.LatNx .AND. idz.NE.1) then
            DEM(I)%NeighborID(19) = I + 1 - LatNx*LatNy
        end if
        !  Right Below Under lattice
        if (idx.NE.LatNx .AND. idy.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(20) = I + 1 - LatNx - LatNx*LatNy
        end if
        !  Below Under lattice
        if (idy.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(21) = I - LatNx - LatNx*LatNy
        end if
        !  Left Below Under lattice
        if (idx.NE.1 .AND. idy.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(22) = I - 1 - LatNx - LatNx*LatNy
        end if
        !  Left Under lattice
        if (idx.NE.1 .AND. idz.NE.1) then
            DEM(I)%NeighborID(23) = I - 1 - LatNx*LatNy
        end if
        !  Left Upper Under lattice
        if (idx.NE.1 .AND. idy.NE.LatNy .AND. idz.NE.1) then
            DEM(I)%NeighborID(24) = I - 1 + LatNx - LatNx*LatNy
        end if
        !  Upper Under lattice
        if (idy.NE.LatNy .AND. idz.NE.1) then
            DEM(I)%NeighborID(25) = I + LatNx - LatNx*LatNy
        end if
        !  Right Upper Under lattice
        if (idx.NE.LatNx .AND. idy.NE.LatNy .AND. idz.NE.1) then
            DEM(I)%NeighborID(26) = I + 1 + LatNx - LatNx*LatNy
        end if
    end do

    !  Initialize Neighbor Nodelink.
    allocate(IDInner(LatNum))
    allocate(tailInner(LatNum))
    do I = 1,LatNum
        IDInner(I)%No = 0                 !  Array of Lattice's inner particles
        nullify(IDInner(I)%next)
    end do 
#endif

    !  Initialize Hertz list. Linklist-Compressed.
    write(*,*) "< Hertz list initializing..."
    
    allocate(Head(NMAX))
    do I = 1,NMAX
        Head(I)%No = 0                    !  Length of Linklist
        !Head(I)%recordTime = 0.0D0        !  record time
        Head(I)%Hertz(1) = 0.0D0          !  Hertz(1)
        Head(I)%Hertz(2) = 0.0D0          !  Hertz(2)
        Head(I)%Hertz(3) = 0.0D0          !  Hertz(3)
        Head(I)%Mrot(1) = 0.0D0           !  Mrot(1)
        Head(I)%Mrot(2) = 0.0D0           !  Mrot(2)
        Head(I)%Mrot(3) = 0.0D0           !  Mrot(3)
        Head(I)%Mtwist(1) = 0.0D0         !  Mtwist(1)
        Head(I)%Mtwist(2) = 0.0D0         !  Mtwist(2)
        Head(I)%Mtwist(3) = 0.0D0         !  Mtwist(3)
        Head(I)%is_touching = .false.     !  Whether touch or not
        Head(I)%is_slipping = .false.     !  Whether slip or not
        Head(I)%is_rolling = .false.      !  Whether roll or not
        Head(I)%is_twisting = .false.     !  Whether twist or not
        nullify(Head(I)%prev)             !  Point to Prev 
        nullify(Head(I)%next)             !  Point to Next
    end do

    return
    end