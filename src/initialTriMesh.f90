    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     Initialization of TriMesh.
    !     Allocate trimesh into DEM lattice
    !     ---------------------------------
    !
    !********************************************************************
    SUBROUTINE initialTriMesh()
    
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
                    !  check whether this lattice is near the trimesh
                    if (RVL .LE. 2.0*LatDx) then
                        index = (idRange(1,1)+J-1) + (idRange(1,2)+K-1-1)*LatNx + (idRange(1,3)+L-1-1)*LatNx*LatNy
                        trimeshDEM(index)%No = trimeshDEM(index)%No + 1
                        !  search from the beginning，should be slow and replaced by tail method
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
    
    end