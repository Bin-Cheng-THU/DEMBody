    !********************************************************************
    !     DEMBody 4.4
    !     ***********
    !
    !     Grav TriMesh force.
    !     --------------------------
    !
    !     @Gravity of Grav TriMesh
    !
    !     Note that considering the movement of Pan nucleus, 
    !     the gravity of Pan must be transferred from the origin point.
    !     So we use Gravity Body instead.
    !
    !********************************************************************
    subroutine forceGravTriMesh()

    use global
    implicit none

    integer I,J,K
    integer :: index                 !  Index of the GravTriMeshGrid containing the Particle
    real(8) :: gravNode(3,8)         !  gravity of Nodes
    real(8) :: coord(3)              !  coordinates of the Particle
    real(8) :: shapeFunction(8)      !  shape function in Finite Element Method

    !$OMP PARALLEL DO PRIVATE(I,J,K,index,gravNode,coord,shapeFunction)
    do I = 1,N
        
        index = floor((X(1,I)+TriLatMx)/TriLatDx) + 1 + floor((X(2,I)+TriLatMy)/TriLatDy)*TriLatNx + floor((X(3,I)+TriLatMz)/TriLatDz)*(TriLatNx*TriLatNy)
        
        do J = 1,8
            do K = 1,3
                gravNode(K,J) = GravTriMeshNode(K,GravTriMeshGrid(index)%NodeID(J))
            end do
        end do
        
        coord(1) = (X(1,I)-GravTriMeshGrid(index)%PositionCenter(1))/TriLatDx*2.0D0
        coord(2) = (X(2,I)-GravTriMeshGrid(index)%PositionCenter(2))/TriLatDy*2.0D0
        coord(3) = (X(3,I)-GravTriMeshGrid(index)%PositionCenter(3))/TriLatDz*2.0D0
        
        shapeFunction(1) = -(coord(1)-1.0D0)*(coord(2)-1.0D0)*(coord(3)-1.0D0)/8.0D0
        shapeFunction(2) =  (coord(1)+1.0D0)*(coord(2)-1.0D0)*(coord(3)-1.0D0)/8.0D0
        shapeFunction(3) = -(coord(1)+1.0D0)*(coord(2)+1.0D0)*(coord(3)-1.0D0)/8.0D0
        shapeFunction(4) =  (coord(1)-1.0D0)*(coord(2)+1.0D0)*(coord(3)-1.0D0)/8.0D0
        shapeFunction(5) =  (coord(1)-1.0D0)*(coord(2)-1.0D0)*(coord(3)+1.0D0)/8.0D0
        shapeFunction(6) = -(coord(1)+1.0D0)*(coord(2)-1.0D0)*(coord(3)+1.0D0)/8.0D0
        shapeFunction(7) =  (coord(1)+1.0D0)*(coord(2)+1.0D0)*(coord(3)+1.0D0)/8.0D0
        shapeFunction(8) = -(coord(1)-1.0D0)*(coord(2)+1.0D0)*(coord(3)+1.0D0)/8.0D0
        
        do K = 1,3
            F(K,I) = F(K,I) + shapeFunction(1)*gravNode(K,1) + shapeFunction(2)*gravNode(K,2) + shapeFunction(3)*gravNode(K,3) &
            + shapeFunction(4)*gravNode(K,4) + shapeFunction(5)*gravNode(K,5) + shapeFunction(6)*gravNode(K,6) &
            + shapeFunction(7)*gravNode(K,7) + shapeFunction(8)*gravNode(K,8)
        end do     
    end do
    !$OMP END PARALLEL DO        

    return
    end