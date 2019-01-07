!!------------------------------------------------------------------------------
!!----- Pressure non-equilibrium multi-material (Reconstruction-limiting module)
!!----- by
!!----- Aditya K Pandare
!!------------------------------------------------------------------------------

MODULE reconstruction_mm6eq

USE glob_var
USE gquadrature
USE eos

implicit none

CONTAINS

!-------------------------------------------------------------------------------
!----- P0 reconstruction (does nothing):
!-------------------------------------------------------------------------------

subroutine reconstruction_p0(ucons)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  ie = 0

end subroutine reconstruction_p0

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

  !--- limit reconstructed solution
  call limiting_p1(ucons)

end subroutine reconstruction_p0p1

!-------------------------------------------------------------------------------
!----- P1 reconstruction (only limiting):
!-------------------------------------------------------------------------------

subroutine reconstruction_p1(ucons)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  !--- limit reconstructed solution
  call limiting_p1(ucons)

end subroutine reconstruction_p1

!-------------------------------------------------------------------------------
!----- P1P2 reconstruction:
!-------------------------------------------------------------------------------

subroutine reconstruction_p1p2(ucons)

integer :: ig, j, ie, je, ieqn, ngauss
data       ngauss/2/

real*8  :: dxi, dxj, xci, xcj, xg, wi, &
           carea(2), weight(2), &
           b2, b3, b2t, b3t, &
           r1, r2, rhs, lhs
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  call rutope(1, ngauss, carea, weight)

  do ieqn = 1,g_neqns

  !--- 2-exact least-squares reconstruction
  do ie = 1,imax

    dxi = coord(ie+1)-coord(ie)
    xci = 0.5*(coord(ie+1)+coord(ie))
    lhs = 0.0
    rhs = 0.0

    ! neighbours
    do j = 1,2
      je = ie + ((-1)**j)

      ! internal cells only
      if ( (je.ge.1) .and. (je.le.imax) ) then

        dxj = coord(je+1)-coord(je)
        xcj = 0.5*(coord(je+1)+coord(je))
        b2t = 0.0
        b3t = 0.0

        ! quadrature
        do ig = 1,ngauss
          xg = carea(ig) * 0.5*dxj + xcj
          wi = 0.5 * weight(ig)
          b2 = p1basis(xg, xci, dxi)
          b3 = p2basis(xg, xci, dxi)
          b2t = b2t + wi * b2
          b3t = b3t + wi * b3
        end do !ig

        ! assemble lhs and rhs
        b2 = p1basis(xcj, xci, dxi)
        lhs = lhs + b3t*b3t + b2*b2
        r1  = ucons(1,ieqn,je) - ( ucons(1,ieqn,ie) + ucons(2,ieqn,ie)*b2t )
        r2  = ucons(2,ieqn,je)*dxi/dxj - ucons(2,ieqn,ie)
        rhs = rhs + b3t*(r1) + b2*(r2)

      end if
    end do !j

    ! get p2 dof
    ucons(3,ieqn,ie) = rhs/lhs

  end do !ie

  !ie = 1
  !dxi = 0.5 * (coord(ie+1)-coord(ie))
  !dxj = 0.5 * (coord(ie+2)-coord(ie+1))
  !ucons(3,ieqn,ie) = (dxi*dxi)*(ucons(2,ieqn,ie+1)/dxj - ucons(2,ieqn,ie)/dxi) / (dxj+dxi)

  !ie = imax
  !dxi = 0.5 * (coord(ie+1)-coord(ie))
  !dxj = 0.5 * (coord(ie)-coord(ie-1))
  !ucons(3,ieqn,ie) = (dxi*dxi)*(ucons(2,ieqn,ie)/dxi - ucons(2,ieqn,ie-1)/dxj) / (dxi+dxj)

  end do !ieqn

  !--- limit reconstructed solution
  call limiting_p2(ucons)
  call limiting_p1(ucons)

end subroutine reconstruction_p1p2

!-------------------------------------------------------------------------------
!----- P1 limiting:
!-------------------------------------------------------------------------------

subroutine limiting_p1(ucons)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  select case (g_nlim)

  case(0)

  case(1)
    call min_superbee(ucons)

  case(2)
    call sconsistent_superbee_p1(ucons)

  case(3)
    call sconsistent_oversuperbee(ucons)

  case(4)
    call weno_p1(ucons)

  case default
    write(*,*) "Error: incorrect p1-limiter index in control file: ", g_nlim
    call exit

  end select

end subroutine limiting_p1

!-------------------------------------------------------------------------------
!----- P2 limiting:
!-------------------------------------------------------------------------------

subroutine limiting_p2(ucons)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  select case (g_nlim)

  case(0)

  case(2)
    call sconsistent_superbee_p2(ucons)

  case default
    write(*,*) "Error: incorrect p2-limiter index in control file: ", g_nlim
    call exit

  end select

end subroutine limiting_p2

!-------------------------------------------------------------------------------
!----- superbee limiter:
!-------------------------------------------------------------------------------

subroutine min_superbee(ucons)

integer :: ie, ieqn, ifc
real*8  :: theta(g_neqns)
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1)

  do ie = 1,imax

    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, 2.0, uneigh, theta)

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
!----- system-consistent superbee limiter for P1 dofs:
!-------------------------------------------------------------------------------

subroutine sconsistent_superbee_p1(ucons)

integer :: ie, ieqn, ifc
real*8  :: al1, theta(g_neqns)
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1)

  do ie = 1,imax

    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, 2.0, uneigh, theta)

    al1 = ucons(1,1,ie)

    ! 2. Obtain consistent limiter functions for the equation system
    !    Interface detection
    if ( (g_nmatint .eq. 1) .and. &
         (al1 .gt. 10.0*g_alphamin) .and. (al1 .lt. 1.0-10.0*g_alphamin) ) then

      call intfac_limiting(ucons(:,:,ie), 0.0, 0.0, theta, theta(1))

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end subroutine sconsistent_superbee_p1

!-------------------------------------------------------------------------------
!----- system-consistent superbee limiter for P2 dofs:
!-------------------------------------------------------------------------------

subroutine sconsistent_superbee_p2(ucons)

integer :: ie, ieqn, ifc
real*8  :: dx2, al1, theta(g_neqns)
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1)

  do ie = 1,imax

    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / dx2

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,:,0)  = ucons(2:3,:,ie) / dx2

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / dx2

    ! 1. compute limiter function
    call superbee_fn(g_neqns, 2.0, uneigh, theta)

    al1 = ucons(1,1,ie)

    ! 2. Obtain consistent limiter functions for the equation system
    !    Interface detection
    if ( (g_nmatint .eq. 1) .and. &
         (al1 .gt. 10.0*g_alphamin) .and. (al1 .lt. 1.0-10.0*g_alphamin) ) then

      call intfac_limiting_p2(ucons(:,:,ie), theta(1))

    else

      do ieqn = 1,g_neqns
        ucons(3,ieqn,ie) = theta(ieqn) * ucons(3,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end subroutine sconsistent_superbee_p2

!-------------------------------------------------------------------------------
!----- system-consistent overbee+superbee limiter:
!-------------------------------------------------------------------------------

subroutine sconsistent_oversuperbee(ucons)

integer :: ie, ieqn, ifc
real*8  :: theta(g_neqns), theta_al, thrho(2)
real*8  :: al1, rho1, rho2, vel, rhoe1, rhoe2
real*8  :: rhoneigh(2,2,-1:1), &
           uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1)

  do ie = 1,imax

    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, 2.0, uneigh, theta)

    al1 = ucons(1,1,ie)

    ! 2. Obtain consistent limiter functions for the equation system
    !    Interface detection
    if ( (g_nmatint .eq. 1) .and. &
         (al1 .gt. 10.0*g_alphamin) .and. (al1 .lt. 1.0-10.0*g_alphamin) ) then

!      rhoneigh(1,1,-1) = ucons(1,2,ie-1)/ucons(1,1,ie-1)
!      rhoneigh(1,1,0)  = ucons(1,2,ie)/ucons(1,1,ie)
!      rhoneigh(1,1,1)  = ucons(1,2,ie+1)/ucons(1,1,ie+1)
!
!      rhoneigh(1,2,-1) = ucons(1,3,ie-1)/(1.0-ucons(1,1,ie-1))
!      rhoneigh(1,2,0)  = ucons(1,3,ie)/(1.0-ucons(1,1,ie))
!      rhoneigh(1,2,1)  = ucons(1,3,ie+1)/(1.0-ucons(1,1,ie+1))
!
!      rhoneigh(2,1,0) = ( theta(2)*ucons(2,2,ie) &
!                        - rhoneigh(1,1,0)*theta(1)*ucons(2,1,ie) ) / al1
!      rhoneigh(2,2,0) = ( theta(3)*ucons(2,3,ie) &
!                        + rhoneigh(1,2,0)*theta(1)*ucons(2,1,ie) ) / (1.0-al1)
!
!      call superbee_fn(2, 1.0, rhoneigh, thrho)
!
!      rhoneigh(2,1,0) = thrho(1) * rhoneigh(2,1,0)
!      rhoneigh(2,2,0) = thrho(2) * rhoneigh(2,2,0)
!
      rhoneigh(2,1,0) = 0.0
      rhoneigh(2,2,0) = 0.0

      !   compressive limiting for volume fraction
      call overbee_fn(uneigh, theta_al)
      call intfac_limiting(ucons(:,:,ie), rhoneigh(2,1,0), rhoneigh(2,2,0), theta, theta_al)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end subroutine sconsistent_oversuperbee

!-------------------------------------------------------------------------------
!----- Superbee limiter for n equations individually
!----- this sub actually calculates the limiter function according to superbee.
!----- the input to this function can be modified to limit P1 or P2 dofs, refer
!----- to sconsistent_superbee_p1 and sconsistent_superbee_p2 respectively
!-------------------------------------------------------------------------------

subroutine superbee_fn(neq,beta_lim,ucons,theta)

integer, intent(in) :: neq
real*8,  intent(in) :: beta_lim, ucons(2,neq,-1:1)

integer :: ieqn, ifc
real*8  :: ui, ug, umin, umax, diff, phi, theta(neq), thetal

  theta(:) = 1.0

  do ieqn = 1,neq

    !--- limiter
    ! beta = 2 : Superbee
    !      = 1 : Minmod

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
!----- Function that consistently applies limiter for near-interface cell
!-------------------------------------------------------------------------------

subroutine intfac_limiting(ucons, drho1dx, drho2dx, theta, theta_al)

real*8,  intent(in) :: theta(g_neqns), theta_al, drho1dx, drho2dx
real*8  :: al1, rho1, rho2, vel, rhoe1, rhoe2
real*8  :: ucons(g_tdof,g_neqns,1)

  !        get primitive variables
  al1 = ucons(1,1,1)
  rho1 = ucons(1,2,1)/al1
  rho2 = ucons(1,3,1)/(1.0-al1)
  vel  = ucons(1,4,1)/( ucons(1,2,1) + ucons(1,3,1) )
  rhoe1 = ucons(1,5,1)/al1
  rhoe2 = ucons(1,6,1)/(1.0-al1)

  !   i.   Volume fraction: Keep limiter function the same
  ucons(2,1,1) = theta_al * ucons(2,1,1)
  !   ii.  Continuity:
  ucons(2,2,1) = rho1 * ucons(2,1,1) &
                + al1  * drho1dx
  ucons(2,3,1) = - rho2 * ucons(2,1,1) &
                + (1.0-al1) * drho2dx
  !   iii. Momentum:
  ucons(2,4,1) = vel * (ucons(2,2,1) + ucons(2,3,1))
  !   iv.  Energy:
  ucons(2,5,1) = rhoe1 * ucons(2,1,1) &
                + al1 * 0.5*vel*vel * drho1dx
  ucons(2,6,1) = - rhoe2 * ucons(2,1,1) &
                + (1.0-al1) * 0.5*vel*vel * drho2dx

end subroutine intfac_limiting

!-------------------------------------------------------------------------------
!----- Function that consistently applies limiter to P2 dofs for
!----- near-interface cell
!-------------------------------------------------------------------------------

subroutine intfac_limiting_p2(ucons, theta_al)

real*8,  intent(in) :: theta_al
real*8  :: al1, rho1, rho2, vel, rhoe1, rhoe2
real*8  :: ucons(g_tdof,g_neqns,1)

  !        get primitive variables
  al1 = ucons(1,1,1)
  rho1 = ucons(1,2,1)/al1
  rho2 = ucons(1,3,1)/(1.0-al1)
  vel  = ucons(1,4,1)/( ucons(1,2,1) + ucons(1,3,1) )
  rhoe1 = ucons(1,5,1)/al1
  rhoe2 = ucons(1,6,1)/(1.0-al1)

  !   i.   Volume fraction: Keep limiter function the same
  ucons(3,1,1) = theta_al * ucons(3,1,1)
  !   ii.  Continuity:
  ucons(3,2,1) = rho1 * ucons(3,1,1)
  ucons(3,3,1) = - rho2 * ucons(3,1,1)
  !   iii. Momentum:
  ucons(3,4,1) = vel * (ucons(3,2,1) + ucons(3,3,1))
  !   iv.  Energy:
  ucons(3,5,1) = rhoe1 * ucons(3,1,1)
  ucons(3,6,1) = - rhoe2 * ucons(3,1,1)

end subroutine intfac_limiting_p2

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
!----- WENO limiter for P1
!-------------------------------------------------------------------------------

subroutine weno_p1(ucons)

integer :: ie,iel,ier,nsten,is,ieqn
real*8  :: theta(g_neqns), theta_al, al1, dxalp
real*8  :: wi,epsweno,wenocp1,wt, &
           dx,dxl,dxr,weight(3), &
           gradv(3),osc(3),gradu(g_neqns,imax), &
           ucons(g_tdof,g_neqns,0:imax+1)

  epsweno = 1.d-8
  wenocp1 = 10.0

  do ieqn = 1,g_neqns

    !--- determine WENO-limited first derivatives
    do ie = 1,imax

      iel = ie - 1
      ier = ie + 1

      dxl = 0.5 * (coord(ie)-coord(ie-1))
      dx  = 0.5 * (coord(ie+1)-coord(ie))
      dxr = 0.5 * (coord(ie+2)-coord(ie+1))

      !--- the nsten stencils
      gradv(1) = ucons(2,ieqn,ie)/dx

      if (ie == imax) then
        nsten = 2
        gradv(2) = ucons(2,ieqn,iel)/dxl

      elseif (ie == 1) then
        nsten = 2
        gradv(2) = ucons(2,ieqn,ier)/dxr

      else
        nsten = 3
        gradv(2) = ucons(2,ieqn,iel)/dxl
        gradv(3) = ucons(2,ieqn,ier)/dxr

      end if

      !--- oscillation indicators

      do is = 1,nsten
        osc(is) = dsqrt(gradv(is)*gradv(is))
      end do !is

      wt = 0.d0

      do is = 1,nsten
        if (is.eq.1) then
          wi = wenocp1
        else
          wi = 1.d0
        end if !is
        weight(is) = wi* (epsweno + osc(is))**(-2.d0)
        wt = wt + weight(is)
      end do !is

      !--- normalize the weight function
      do is = 1,nsten
        weight(is) = weight(is)/wt
      end do !is

      !--- reconstruct the limited gradient

      gradu(ieqn,ie) = 0.d0

      do is = 1,nsten
        gradu(ieqn,ie) = gradu(ieqn,ie) + weight(is)*gradv(is)
      end do !is

    end do !ie

  end do !ieqn

  !--- modify the first derivatives
  do ie = 1,imax

    dx = 0.5 * (coord(ie+1)-coord(ie))
    gradu(:,ie) = dx * gradu(:,ie)
    al1 = ucons(1,1,ie)

    ! Obtain consistent limiter functions for the equation system
    ! Interface detection
    if ( (g_nmatint .eq. 1) .and. &
         (al1 .gt. 10.0*g_alphamin) .and. (al1 .lt. 1.0-10.0*g_alphamin) ) then

      dxalp = ucons(2,1,ie)
      theta_al = gradu(1,ie)/( dxalp + dsign(1.0d-12,dxalp) )
      theta = 0.0
      call intfac_limiting(ucons(:,:,ie), 0.0, 0.0, theta, theta_al)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = gradu(ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end subroutine weno_p1

!-------------------------------------------------------------------------------

END MODULE reconstruction_mm6eq
