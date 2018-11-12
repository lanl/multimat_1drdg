!!------------------------------------------------------------------------------
!!----- Pressure non-equilibrium multi-material (Reconstruction-limiting module)
!!----- by
!!----- Aditya K Pandare
!!------------------------------------------------------------------------------

MODULE reconstruction_mm6eq

USE glob_var
USE eos

implicit none

CONTAINS

!-------------------------------------------------------------------------------
!----- P0P1 reconstruction:
!-------------------------------------------------------------------------------

subroutine reconstruction_p0p1(ucons)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  !--- central difference reconstruction (least-squares for uniform meshes)
  do ie = 1,imax
    ucons(2,:,ie) = 0.25 * (ucons(1,:,ie+1) - ucons(1,:,ie-1))
  end do !ie

  if (g_nlim .eq. 1) then
    call min_superbee(ucons)

  elseif (g_nlim .eq. 2) then
    call sconsistent_superbee(ucons)

  elseif (g_nlim .eq. 3) then
    call sconsistent_oversuperbee(ucons)

  elseif (g_nlim .ne. 0) then
    write(*,*) "Error: incorrect limiter index in control file: ", g_nlim
    stop

  end if

end subroutine reconstruction_p0p1

!-------------------------------------------------------------------------------
!----- superbee limiter:
!-------------------------------------------------------------------------------

subroutine min_superbee(ucons)

integer :: ie, ieqn, ifc
real*8  :: theta(g_neqns)
real*8  :: uneigh(g_tdof,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1)

  do ie = 1,imax

    uneigh(:,:,-1) = ucons(:,:,ie-1)
    uneigh(:,:,0)  = ucons(:,:,ie)
    uneigh(:,:,1)  = ucons(:,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, uneigh, theta)

     do ieqn = 1,g_neqns
       theta(ieqn) = minval(theta)
     end do !ieqn

    ! 2. limit 2nd dofs
    do ieqn = 1,g_neqns
      ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
    end do !ieqn

  end do !ie

end subroutine min_superbee

!-------------------------------------------------------------------------------
!----- system-consistent superbee limiter:
!-------------------------------------------------------------------------------

subroutine sconsistent_superbee(ucons)

integer :: ie, ieqn, ifc
real*8  :: theta(g_neqns)
real*8  :: al1, rho1, rho2, vel, rhoe1, rhoe2
real*8  :: uneigh(g_tdof,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1)

  do ie = 1,imax

    uneigh(:,:,-1) = ucons(:,:,ie-1)
    uneigh(:,:,0)  = ucons(:,:,ie)
    uneigh(:,:,1)  = ucons(:,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, uneigh, theta)

    al1 = ucons(1,1,ie)

    ! 2. Obtain consistent limiter functions for the equation system
    !    Interface detection
    if ( (al1 .gt. 2.0*alphamin) .and. (al1 .lt. 1.0-2.0*alphamin) ) then

      !        get primitive variables
      rho1 = ucons(1,2,ie)/al1
      rho2 = ucons(1,3,ie)/(1.0-al1)
      vel  = ucons(1,4,ie)/( ucons(1,2,ie) + ucons(1,3,ie) )
      rhoe1 = ucons(1,5,ie)/al1
      rhoe2 = ucons(1,6,ie)/(1.0-al1)

      !   i.   Volume fraction: Keep limiter function the same
      ucons(2,1,ie) = theta(1) * ucons(2,1,ie)
      !   ii.  Continuity:
      ucons(2,2,ie) = rho1 * ucons(2,1,ie)
      ucons(2,3,ie) = - rho2 * ucons(2,1,ie)
      !   iii. Momentum:
      ucons(2,4,ie) = vel * (ucons(2,2,ie) + ucons(2,3,ie))
      !   iv.  Energy:
      ucons(2,5,ie) = rhoe1 * ucons(2,1,ie)
      ucons(2,6,ie) = - rhoe2 * ucons(2,1,ie)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end subroutine sconsistent_superbee

!-------------------------------------------------------------------------------
!----- system-consistent overbee+superbee limiter:
!-------------------------------------------------------------------------------

subroutine sconsistent_oversuperbee(ucons)

integer :: ie, ieqn, ifc
real*8  :: theta(g_neqns), theta_al
real*8  :: al1, rho1, rho2, vel, rhoe1, rhoe2
real*8  :: uneigh(g_tdof,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1)

  do ie = 1,imax

    uneigh(:,:,-1) = ucons(:,:,ie-1)
    uneigh(:,:,0)  = ucons(:,:,ie)
    uneigh(:,:,1)  = ucons(:,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, uneigh, theta)

    al1 = ucons(1,1,ie)

    ! 2. Obtain consistent limiter functions for the equation system
    !    Interface detection
    if ( (al1 .gt. 10.0*alphamin) .and. (al1 .lt. 1.0-10.0*alphamin) ) then
!    if ( al1*(1.0-al1) .gt. 0.1 ) then

      !   0.   compressive limiting for volume fraction
      call overbee_fn(uneigh, theta_al)

      !        Get primitive variables
      rho1 = ucons(1,2,ie)/al1
      rho2 = ucons(1,3,ie)/(1.0-al1)
      vel  = ucons(1,4,ie)/( ucons(1,2,ie) + ucons(1,3,ie) )
      rhoe1 = ucons(1,5,ie)/al1
      rhoe2 = ucons(1,6,ie)/(1.0-al1)

      !   i.   Volume fraction: Keep limiter function the same
      ucons(2,1,ie) = theta_al * ucons(2,1,ie)
      !   ii.  Continuity:
      ucons(2,2,ie) = rho1 * ucons(2,1,ie)
      ucons(2,3,ie) = - rho2 * ucons(2,1,ie)
      !   iii. Momentum:
      ucons(2,4,ie) = vel * (ucons(2,2,ie) + ucons(2,3,ie))
      !   iv.  Energy:
      ucons(2,5,ie) = rhoe1 * ucons(2,1,ie)
      ucons(2,6,ie) = - rhoe2 * ucons(2,1,ie)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end subroutine sconsistent_oversuperbee

!-------------------------------------------------------------------------------
!----- Superbee limiter for n equations individually
!----- this sub actually calculates the limiter function according to superbee
!-------------------------------------------------------------------------------

subroutine superbee_fn(neq,ucons,theta)

integer, intent(in) :: neq
real*8,  intent(in) :: ucons(g_tdof,g_neqns,-1:1)

integer :: ieqn, ifc
real*8  :: ui, ug, umin, umax, diff, phi, theta(neq), thetal, beta_lim

  theta(:) = 1.0

  do ieqn = 1,neq

    !--- limiter
    ! beta = 2 : Superbee
    !      = 1 : Minmod
    beta_lim = 1.0

    ui = ucons(1,ieqn,0)

    ! find min and max in neighborhood
    umax = max( max(ucons(1,ieqn,-1), ui), ucons(1,ieqn,1) )
    umin = min( min(ucons(1,ieqn,-1), ui), ucons(1,ieqn,1) )

    do ifc = 1,2
      ! unlimited 2nd order solution
      ug = ucons(1,ieqn,0) + ((-1.0)**ifc) * ucons(2,ieqn,0)

      ! bounds
      diff = ug-ui
      if (diff > 1.0d-16) then
        phi = (umax-ui)/(2.0*diff)

      else if (diff < -1.0d-16) then
        phi = (umin-ui)/(2.0*diff)

      else
        phi = 1.0

      end if

      ! limiter function
      thetal = max( 0.0, max( min(beta_lim*phi, 1.0), min(phi, beta_lim) ) )
      theta(ieqn) = min(thetal, theta(ieqn))

    end do !ifc

  end do !ieqn

end subroutine superbee_fn

!-------------------------------------------------------------------------------
!----- Overbee limiter for n equations individually
!----- this sub actually calculates the limiter function according to overbee
!-------------------------------------------------------------------------------

subroutine overbee_fn(ucons,theta)

real*8,  intent(in) :: ucons(g_tdof,g_neqns,-1:1)

integer :: ifc
real*8  :: ui, ug, umin, umax, diff, phi, theta, thetal

  theta = 2.0

  ui = ucons(1,1,0)

  ! find min and max in neighborhood
  umax = max( max(ucons(1,1,-1), ui), ucons(1,1,1) )
  umin = min( min(ucons(1,1,-1), ui), ucons(1,1,1) )

  do ifc = 1,2
    ! unlimited 2nd order solution
    ug = ucons(1,1,0) + ((-1.0)**ifc) * ucons(2,1,0)

    ! bounds
    diff = ug-ui
    if (diff > 1.0d-16) then
      phi = (umax-ui)/(2.0*diff)

    else if (diff < -1.0d-16) then
      phi = (umin-ui)/(2.0*diff)

    else
      phi = 1.0

    end if

    ! limiter function
    thetal = max( 0.0, min(2.0*phi, 2.0) )
    theta = min(thetal, theta)

  end do !ifc

end subroutine overbee_fn

!-------------------------------------------------------------------------------

END MODULE reconstruction_mm6eq
