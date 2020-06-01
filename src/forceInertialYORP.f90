    !********************************************************************
    !     DEMBody 8.4
    !     ***********
    !
    !     YORP system force in inertial frame.
    !     --------------------------
    !
    !     @Gravity
    !     @No non-inertial forces
    !
    !     Note that here we use a simply sphere to conduct the gravity of
    !     the asteroid interior. It should be replaced with more accurate
    !     method instead (maybe polyhedron method?).
    !
    !     Self-gravity between large boulders are calculated by Traverse Method;
    !     We ignore the self-gravity between fine particles, and the gravity between
    !     these particles and large boulders in interior is calculated by sphere-assumption.
    !     The position of this sphere is refreshed according to the COM of large boulders. 
    !
    !********************************************************************
    subroutine forceInertialYORP()

    use global
    implicit none

    integer I,J,K
    
    real(8) :: dist(3),distS,distL,center

    !$OMP PARALLEL DO PRIVATE(I,K,dist,distS,distL,center)
    do I = 1,N
        !  Distance vector
        do K = 1,3
            dist(K) = X(K,I) - inertialYORPcenter(K)
        end do
        distL = dist(1)**2 + dist(2)**2 + dist(3)**2
        distS = sqrt(distL)
        if (distS .GE. inertialYORPradius) then
            center = GravConst*inertialYORPmass/distL/distS
        else
            center = 4.0D0/3.0D0*PI*GravConst*inertialYORPdensity
        end if
        do K = 1,3
            F(K,I) = F(K,I) - center*dist(K)
        end do
    end do
    !$OMP END PARALLEL DO

    return
    end