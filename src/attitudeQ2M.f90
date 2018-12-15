    !********************************************************************
    !     DEMBody 5.0
    !     ***********
    !
    !     N-body integrator of Quaternion.
    !     -------------------------------
    !     Quaternion to Attitute Matrix
    !
    !     @Refresh Attitude Matrix
    !     
    !
    !********************************************************************
    subroutine attitudeQ2M(t_quaternion,t_MatI,t_MatB)
    
    implicit none
    real(8) :: t_quaternion(4)
    real(8) :: t_MatI(3,3)
    real(8) :: t_MatB(3,3)

    !  Attitude Matrix in inertial frame
    t_MatI(1,1) = t_quaternion(1)**2-t_quaternion(2)**2-t_quaternion(3)**2+t_quaternion(4)**2
    t_MatI(2,1) = 2.0D0*(t_quaternion(1)*t_quaternion(2)-t_quaternion(3)*t_quaternion(4))
    t_MatI(3,1) = 2.0D0*(t_quaternion(1)*t_quaternion(3)+t_quaternion(2)*t_quaternion(4))
    t_MatI(1,2) = 2.0D0*(t_quaternion(1)*t_quaternion(2)+t_quaternion(3)*t_quaternion(4))
    t_MatI(2,2) = -t_quaternion(1)**2+t_quaternion(2)**2-t_quaternion(3)**2+t_quaternion(4)**2
    t_MatI(3,2) = 2.0D0*(t_quaternion(2)*t_quaternion(3)-t_quaternion(1)*t_quaternion(4))
    t_MatI(1,3) = 2.0D0*(t_quaternion(1)*t_quaternion(3)-t_quaternion(2)*t_quaternion(4))
    t_MatI(2,3) = 2.0D0*(t_quaternion(2)*t_quaternion(3)+t_quaternion(1)*t_quaternion(4))
    t_MatI(3,3) = -t_quaternion(1)**2-t_quaternion(2)**2+t_quaternion(3)**2+t_quaternion(4)**2

    !  Attitude Matrix in body frame
    t_MatB(1,1) = t_quaternion(1)**2-t_quaternion(2)**2-t_quaternion(3)**2+t_quaternion(4)**2
    t_MatB(2,1) = 2.0D0*(t_quaternion(1)*t_quaternion(2)+t_quaternion(3)*t_quaternion(4))
    t_MatB(3,1) = 2.0D0*(t_quaternion(1)*t_quaternion(3)-t_quaternion(2)*t_quaternion(4))
    t_MatB(1,2) = 2.0D0*(t_quaternion(1)*t_quaternion(2)-t_quaternion(3)*t_quaternion(4))
    t_MatB(2,2) = -t_quaternion(1)**2+t_quaternion(2)**2-t_quaternion(3)**2+t_quaternion(4)**2
    t_MatB(3,2) = 2.0D0*(t_quaternion(2)*t_quaternion(3)+t_quaternion(1)*t_quaternion(4))
    t_MatB(1,3) = 2.0D0*(t_quaternion(1)*t_quaternion(3)+t_quaternion(2)*t_quaternion(4))
    t_MatB(2,3) = 2.0D0*(t_quaternion(2)*t_quaternion(3)-t_quaternion(1)*t_quaternion(4))
    t_MatB(3,3) = -t_quaternion(1)**2-t_quaternion(2)**2+t_quaternion(3)**2+t_quaternion(4)**2

    end subroutine