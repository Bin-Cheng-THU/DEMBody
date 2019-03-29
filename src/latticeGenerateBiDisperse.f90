    !********************************************************************
    !     DEMBody 5.1
    !     ***********
    !
    !     Generate lattice based on partition of BiDisperse particles.
    !     ---------------------------------
    !     @Allocate biDisperse particles into Paralle Lattice.
    !
    !     @Strategy: allocate biDisperse particles into Parallel Lattice, so for particles in Lattice
    !     @    A, you just need to search contact of biDisperse particles in Lattice A.
    !
    !********************************************************************
    SUBROUTINE latticeGenerateBiDisperse()

    use global
    use omp_lib
    implicit none

    integer I,J,K,L
    integer :: TagCenter
    integer :: idCenter(3)
    integer :: idRange(2,3)
    integer :: index
    
    real(8) :: o1,o2 
    type(biDisperseLattice),pointer :: Temp
    type(biDisperseLattice),pointer :: TempH

    !o1 = omp_get_wtime()
    !  initialized Parallel Lattice
    !$OMP PARALLEL DO PRIVATE(I,Temp) SCHEDULE(GUIDED)
    do I = 1,LatNum
        !  reset Inner quantities    
        do while(associated(biDisperseDEM(I)%next))
            Temp => biDisperseDEM(I)%next
            biDisperseDEM(I)%next => biDisperseDEM(I)%next%next
            deallocate(Temp)
        end do
        biDisperseDEM(I)%No = 0
        !  assign tail of biDisperseDEM to tail
        biDisperseDEMtail(I)%next => biDisperseDEM(I)
    end do
    !$OMP END PARALLEL DO
    !o2 = omp_get_wtime()
    !write(*,*) "Lattice biDisperseDEM Empty",(o2-o1)
    
    !o1 = omp_get_wtime()
    do I = 1,biDisperseNum
        !  find id range
        TagCenter = floor((biDisperseX(1,I)+LatMx)/LatDx) + 1 + floor((biDisperseX(2,I)+LatMy)/LatDy)*LatNx + floor((biDisperseX(3,I)+LatMz)/LatDz)*(LatNx*LatNy)
        !  Tag to id
        idCenter(3) = int((TagCenter-1)/(LatNx*LatNy))+1
        idCenter(2) = int((TagCenter-(idCenter(3)-1)*(LatNx*LatNy)-1)/LatNx)+1
        idCenter(1) = TagCenter-(idCenter(3)-1)*(LatNx*LatNy)-(idCenter(2)-1)*LatNx
        !  id range
        idRange(1,1) = max(idCenter(1)-biDisperseScale-1,1)
        idRange(2,1) = min(idCenter(1)+biDisperseScale+1,LatNx)
        idRange(1,2) = max(idCenter(2)-biDisperseScale-1,1)
        idRange(2,2) = min(idCenter(2)+biDisperseScale+1,LatNy)
        idRange(1,3) = max(idCenter(3)-biDisperseScale-1,1)
        idRange(2,3) = min(idCenter(3)+biDisperseScale+1,LatNz)
        do J = 1,(idRange(2,1)-idRange(1,1)+1)
            do K = 1,(idRange(2,2)-idRange(1,2)+1)
                do L = 1,(idRange(2,3)-idRange(1,3)+1)
                    index = (idRange(1,1)+J-1) + (idRange(1,2)+K-1-1)*LatNx + (idRange(1,3)+L-1-1)*LatNx*LatNy
                    biDisperseDEM(index)%No = biDisperseDEM(index)%No + 1
                    !TempH => biDisperseDEM(index)
                    !do while (associated(TempH%next))
                    !    TempH => TempH%next
                    !end do
                    allocate(Temp)
                    Temp = biDisperseLattice(I,NULL())
                    !  Using the order of biDisperse particles to accelerate the Lattice generate (same strategy as in latticeGenerate)
                    biDisperseDEMtail(index)%next%next => Temp
                    biDisperseDEMtail(index)%next => Temp
                end do
            end do     
        end do
    end do
    !o2 = omp_get_wtime()
    !write(*,*) "Lattice biDisperseDEM generate",(o2-o1)    
    
    !  reset Position array
    !$OMP PARALLEL DO PRIVATE(I,K)
    do I = 1,biDisperseNum
        do K = 1,3
            biDisperseXT(K,I) = 0.0D0
        end do
    end do
    !$OMP END PARALLEL DO    

    end