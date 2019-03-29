    !********************************************************************
    !     DEMBody 5.1
    !     ***********
    !
    !     Initialization of bondedTriMesh.
    !     ---------------------------------
    !
    !********************************************************************
    SUBROUTINE initialBondedTriMesh()

    use global
    implicit none

    integer I,J,K,L
    real(8) :: pointTriMesh(3,3)
    integer :: TagTriMesh(3)
    integer :: idTriMesh(3,3)
    integer :: idRange(2,3)
    real(8) :: centerLattice(3),RV(3),RVL
    integer :: index
    type(trimeshLattice),pointer :: Temp
    type(trimeshLattice),pointer :: TempH

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
            !  in case of TriMeshes out of mesh
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
        idRange(2,3) = min(idRange(2,3)+1,LatNz)
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

    end