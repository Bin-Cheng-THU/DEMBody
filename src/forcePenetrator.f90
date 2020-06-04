    !********************************************************************
    !     DEMBody 6.4
    !     ***********
    !
    !     Force for penetrators.
    !     --------------------------
    !      
    !     @Using rolling friction similar to LIGGGHTS
    !     @Using Hertz-Mindlin contact model similar to Wada
    !     @No Cohesion model in walls
    !     @No Rolling model in bonded walls
    !     @Using damping model applicable for ice ball
    !     
    !     @1 for cylindical penetrator
    !     @2 for conical penetrator
    !     @3 for hemispherical penetrator
    !
    !********************************************************************
    subroutine forcePenetrator()

    use global
    implicit none

    select case (peneType)
    case (1)
    call forcePeneCylinder

    case (2)
    call forcePeneCylinder

    case (3)
    call forcePeneCylinder

    case default
    write (*,*) 'Unknow Penetrator Type!'
    stop

    end select

    return
    end
