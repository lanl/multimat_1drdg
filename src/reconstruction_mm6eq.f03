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

subroutine reconstruction_p0(ucons, uprim)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  !--- reconstruct primitives from first-order solution
  call recons_primitives(ucons, uprim)

end subroutine reconstruction_p0

!-------------------------------------------------------------------------------
!----- P0P1 reconstruction:
!-------------------------------------------------------------------------------

subroutine reconstruction_p0p1(ucons, uprim)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  !--- central difference reconstruction (least-squares for uniform meshes)
  do ie = 1,imax
    ucons(2,:,ie) = 0.25 * (ucons(1,:,ie+1) - ucons(1,:,ie-1))
  end do !ie

  !--- reconstruct primitives from second-order solution
  call recons_primitives(ucons, uprim)

  call limiting_p1(ucons, uprim)

end subroutine reconstruction_p0p1

!-------------------------------------------------------------------------------
!----- P1 reconstruction (only limiting):
!-------------------------------------------------------------------------------

subroutine reconstruction_p1(ucons, uprim)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  !--- reconstruct primitives from second-order solution
  call recons_primitives(ucons, uprim)

  call limiting_p1(ucons, uprim)

end subroutine reconstruction_p1

!-------------------------------------------------------------------------------
!----- P1P2 reconstruction:
!-------------------------------------------------------------------------------

subroutine reconstruction_p1p2(ucons, uprim)

integer :: ig, j, ie, je, ieqn, ngauss
data       ngauss/2/

real*8  :: dxi, dxj, xci, xcj, xg, wi, &
           carea(2), weight(2), &
           b2, b3, b2t, b3t, &
           r1, r2, rhs, lhs
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

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

  !--- reconstruct primitives from third-order solution
  call recons_primitives(ucons, uprim)

  call limiting_p2(ucons, uprim)

end subroutine reconstruction_p1p2

!-------------------------------------------------------------------------------
!----- Reconstructing primitives from DG solution:
!-------------------------------------------------------------------------------

subroutine recons_primitives(ucons, uprim)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: ie, ifc, imat, ieqn
real*8  :: rhoavg, drhodx, uface(g_neqns), up_face(g_neqns), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    !--- 1. cell-average bulk velocity, material pressures and material momenta
    uface(:) = ucons(1,:,ie)
    call get_uprim_mm6eq(uface, up_face)
    rhoavg = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ie))
    uprim(1,vel_idx(nummat, 0),ie) = ucons(1,g_mmi%imome,ie) / rhoavg
    do imat = 1,nummat
      uprim(1,apr_idx(nummat, imat),ie) = ucons(1,imat,ie) &
        * up_face(g_mmi%irmin+imat-1)
      uprim(1,mmom_idx(nummat, imat),ie) = ucons(1,g_mmi%irmin+imat-1,ie) &
        * uprim(1,vel_idx(nummat, 0),ie)
    end do !imat

    !--- 2. "gradients" of material pressures, bulk velocity and material momenta
    if (g_nsdiscr .gt. 0) then
    uprim(2,:,ie) = 0.0
    do ifc = 1,2

      !--- rdgp0p1 or dgp1
      if ((g_nsdiscr .eq. 11) .or. (g_nsdiscr .eq. 1)) then
        do ieqn = 1,g_neqns
          uface(ieqn) = ucons(1,ieqn,ie) + ((-1.0)**ifc) * ucons(2,ieqn,ie)
        end do !ieqn

      !--- rdgp1p2
      elseif (g_nsdiscr .eq. 12) then
        do ieqn = 1,g_neqns
          uface(ieqn) = ucons(1,ieqn,ie) + ((-1.0)**ifc) * ucons(2,ieqn,ie) &
                                         + 1.0/3.0*ucons(3,ieqn,ie)
        end do !ieqn

      end if

      call get_uprim_mm6eq(uface, up_face)

      do imat = 1,nummat
        uprim(2,apr_idx(nummat, imat),ie) = uprim(2,apr_idx(nummat, imat),ie) &
          + 0.5 * ((-1.0)**ifc) * uface(imat)*up_face(g_mmi%irmin+imat-1)
      end do !imat

    end do !ifc

    drhodx = sum(ucons(2,g_mmi%irmin:g_mmi%irmax,ie))
    uprim(2,vel_idx(nummat, 0),ie) = ( ucons(2,g_mmi%imome,ie) &
                             - uprim(1,vel_idx(nummat, 0),ie) &
                             * sum(ucons(2,g_mmi%irmin:g_mmi%irmax,ie)) ) &
                           / rhoavg

    do imat = 1,nummat
      uprim(2,mmom_idx(nummat, imat),ie) = &
        uprim(1,vel_idx(nummat, 0),ie) * ucons(2,g_mmi%irmin+imat-1,ie) &
        + ucons(1,g_mmi%irmin+imat-1,ie) * uprim(2,vel_idx(nummat, 0),ie)
    end do !imat

    if (g_nsdiscr .eq. 12) then
      uprim(3,vel_idx(nummat, 0),ie) = ( ucons(3,g_mmi%imome,ie) &
        - (uprim(1,vel_idx(nummat, 0),ie)*sum(ucons(3,g_mmi%irmin:g_mmi%irmax,ie)) &
        + 2.0*drhodx*uprim(2,vel_idx(nummat, 0),ie)) ) / rhoavg
      write(*,*) " FATAL: 2nd order dofs not reconstructed for material momenta!"
      stop 1
    end if
    end if

  end do !ie

end associate

end subroutine recons_primitives

!-------------------------------------------------------------------------------
!----- P0 limiting:
!-------------------------------------------------------------------------------

subroutine limiting_p0(ucons, uprim)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

end subroutine limiting_p0

!-------------------------------------------------------------------------------
!----- P1 limiting:
!-------------------------------------------------------------------------------

subroutine limiting_p1(ucons, uprim)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  select case (g_nlim)

  case(0)

  case(1)
    call min_superbee(ucons, uprim)

  case(2)
    call superbee_p1(ucons, uprim)

  case(3)
    call oversuperbee(ucons, uprim)

  case default
    write(*,*) "Error: incorrect p1-limiter index in control file: ", g_nlim
    call exit

  end select

  call boundpreserve_alpha_p1(ucons)

end subroutine limiting_p1

!-------------------------------------------------------------------------------
!----- P2 limiting:
!-------------------------------------------------------------------------------

subroutine limiting_p2(ucons, uprim)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  select case (g_nlim)

  case(0)

  case(2)
    call superbee_p2(ucons, uprim)

  case(4)
    call weno_p2(ucons, uprim)

  case(5)
    call superbeeweno_p2(ucons, uprim)

  case default
    write(*,*) "Error: incorrect p2-limiter index in control file: ", g_nlim
    call exit

  end select

  call boundpreserve_alpha_p2(ucons)

end subroutine limiting_p2

!-------------------------------------------------------------------------------
!----- Bound preserving limiter for p1:
!-------------------------------------------------------------------------------

subroutine boundpreserve_alpha_p1(ucons)

integer :: ie, ig, ngauss, imat, iamax
real*8  :: careap(2), &
           alm, thal(g_mmi%nummat), thal1, thal2, &
           alm_min, alm_max, eps, al_eps, diff, &
           !arhom, arhom_min, &
           ucons(g_tdof,g_neqns,0:imax+1)

associate (nummat=>g_mmi%nummat)

  eps = 1.0d-16
  al_eps = 0.01*g_alphamin

  ngauss = 2
  careap(1) = -1.0
  careap(2) = 1.0

  do ie = 1,imax
  do imat = 1,nummat

    alm_min   = 10.0
    alm_max   = 0.0001 * eps
    !arhom_min = 1.0d20

    do ig = 1,ngauss

      ! reconstructed volume fraction
      alm = ucons(1,imat,ie) + careap(ig) * ucons(2,imat,ie)
      !arhom = ucons(1,g_mmi%irmin+imat-1,ie) + careap(ig) * ucons(2,g_mmi%irmin+imat-1,ie)

      alm_min = min(alm, alm_min)
      alm_max = max(alm, alm_max)
      !arhom_min = min(arhom, arhom_min)

    end do !ig

    if (alm_min .lt. al_eps) then
    diff = alm_min-ucons(1,imat,ie)
    if (dabs(diff) .le. eps) then
      thal1 = 1.0
    else
      thal1 = dabs((ucons(1,imat,ie)-al_eps)/(diff))
    end if
    else
    thal1 = 1.0
    end if

    if (alm_max .gt. 1.0-al_eps) then
    diff = alm_max-ucons(1,imat,ie)
    if (dabs(diff) .le. eps) then
      thal2 = 1.0
    else
      thal2 = dabs((ucons(1,imat,ie)-(1.0-al_eps))/(diff))
    end if
    else
    thal2 = 1.0
    end if

    thal(imat) = min(1.0, thal1, thal2)

    !diff = arhom_min-ucons(1,g_mmi%irmin+imat-1,ie)
    !if (dabs(diff) .le. eps) then
    !  thal1 = 1.0
    !else
    !  thal1 = dabs((ucons(1,g_mmi%irmin+imat-1,ie))/(diff))
    !end if
    !thal1 = min(1.0, thal1)

    !ucons(2,g_mmi%irmin+imat-1,ie) = thal1 * ucons(2,g_mmi%irmin+imat-1,ie)

  end do !imat

  iamax = maxloc(ucons(1,1:nummat,ie), 1)
  thal(1:nummat) = minval(thal)

  ucons(2,1:nummat,ie) = thal(1:nummat) * ucons(2,1:nummat,ie)

  !if ( (g_nmatint .eq. 1) .and. &
  !     (minval(thal) .lt. 1.0) .and. &
  !     interface_cell(ucons(1,iamax,ie)) ) then

    do imat = 1,nummat
      ucons(2,g_mmi%irmin+imat-1,ie) = thal(imat) * ucons(2,g_mmi%irmin+imat-1,ie)
      ucons(2,g_mmi%iemin+imat-1,ie) = thal(imat) * ucons(2,g_mmi%iemin+imat-1,ie)
    end do !imat

  !end if

  end do !ie

end associate

end subroutine boundpreserve_alpha_p1

!-------------------------------------------------------------------------------
!----- Bound preserving limiter for p2:
!-------------------------------------------------------------------------------

subroutine boundpreserve_alpha_p2(ucons)

integer :: ie, ig, ngauss, imat, iamax
real*8  :: b3, &
           carea(2), weight(2), careap(4), &
           alm, thal(g_mmi%nummat), thal1, thal2, &
           alm_min, alm_max, eps, al_eps, diff, &
           !arhom, arhom_min, &
           ucons(g_tdof,g_neqns,0:imax+1)

associate (nummat=>g_mmi%nummat)

  eps = 1.0d-16
  al_eps = 0.01*g_alphamin

  ngauss = 2
  call rutope(1, ngauss, carea, weight)

  ngauss = 4
  careap(1) = -1.0
  careap(2:3) = carea(1:2)
  careap(4) = 1.0

  do ie = 1,imax
  thal = 1.0

  do imat = 1,nummat

    alm_min   = 10.0
    alm_max   = 0.0001 * eps
    !arhom_min = 1.0d20

    do ig = 1,ngauss

      b3 = 0.5*careap(ig)*careap(ig) - 1.0/6.0

      ! reconstructed volume fraction
      alm = ucons(1,imat,ie) &
            + careap(ig) * ucons(2,imat,ie) &
            + b3 * ucons(3,imat,ie)
      !arhom = ucons(1,g_mmi%irmin+imat-1,ie) &
      !        + careap(ig) * ucons(2,g_mmi%irmin+imat-1,ie) &
      !        + b3 * ucons(3,g_mmi%irmin+imat-1,ie)

      alm_min = min(alm, alm_min)
      alm_max = max(alm, alm_max)
      !arhom_min = min(arhom, arhom_min)

    end do !ig

    if (alm_min .lt. al_eps) then
    diff = alm_min-ucons(1,imat,ie)
    if (dabs(diff) .le. eps) then
      thal1 = 1.0
    else
      thal1 = dabs((ucons(1,imat,ie)-al_eps)/(diff))
    end if
    else
    thal1 = 1.0
    end if

    if (alm_max .gt. 1.0-al_eps) then
    diff = alm_max-ucons(1,imat,ie)
    if (dabs(diff) .le. eps) then
      thal2 = 1.0
    else
      thal2 = dabs((ucons(1,imat,ie)-(1.0-al_eps))/(diff))
    end if
    else
    thal2 = 1.0
    end if

    thal(imat) = min(1.0, thal1, thal2)

    !diff = arhom_min-ucons(1,g_mmi%irmin+imat-1,ie)
    !if (dabs(diff) .le. eps) then
    !  thal1 = 1.0
    !else
    !  thal1 = dabs((ucons(1,g_mmi%irmin+imat-1,ie))/(diff))
    !end if
    !thal1 = min(1.0, thal1)

    !ucons(2,g_mmi%irmin+imat-1,ie) = thal1 * ucons(2,g_mmi%irmin+imat-1,ie)

  end do !imat

  iamax = maxloc(ucons(1,1:nummat,ie), 1)
  thal(1:nummat) = minval(thal)

  do imat = 1,nummat
    ucons(2:3,imat,ie) = thal(imat) * ucons(2:3,imat,ie)
  end do !imat

  !if ( (g_nmatint .eq. 1) .and. &
  !     (minval(thal) .lt. 1.0) .and. &
  !     interface_cell(ucons(1,iamax,ie)) ) then

    do imat = 1,nummat
      ucons(2:3,g_mmi%irmin+imat-1,ie) = thal(imat) * ucons(2:3,g_mmi%irmin+imat-1,ie)
      ucons(2:3,g_mmi%iemin+imat-1,ie) = thal(imat) * ucons(2:3,g_mmi%iemin+imat-1,ie)
    end do !imat

  !end if

  end do !ie

end associate

end subroutine boundpreserve_alpha_p2

!-------------------------------------------------------------------------------
!----- superbee limiter:
!-------------------------------------------------------------------------------

subroutine min_superbee(ucons, uprim)

integer :: ie, ieqn, ifc
real*8  :: theta(g_neqns)
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

  do ie = 1,imax

    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, 2.0, 1.0, uneigh, theta)

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

subroutine superbee_p1(ucons, uprim)

integer :: ie, ieqn, iamax, imat
real*8  :: almax, dalmax, theta(g_neqns), thetap(g_nprim), thetac
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)
logical :: tr_cell(g_neqns)

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    ! 1. detect interface/single-material cell
    iamax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,iamax,ie)
    dalmax = maxval(dabs(ucons(2,1:nummat,ie))/(0.5 * (coord(ie+1)-coord(ie))))

    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then
    ! 2a. Obtain consistent limiter functions for equation system at interface

      uneigh(1:2,1:nummat,-1) = ucons(1:2,1:nummat,ie-1)
      uneigh(1:2,1:nummat,0)  = ucons(1:2,1:nummat,ie)
      uneigh(1:2,1:nummat,1)  = ucons(1:2,1:nummat,ie+1)

      ! i. compute limiter function for volume fractions only
      call superbee_fn(nummat, 2.0, 1.0, uneigh(:,1:nummat,:), theta(1:nummat))

      ! ii. use common limiter function for all volume-fractions
      theta(1:nummat) = theta(iamax) !minval(theta(1:nummat))

      ! iii. consistent limiting
      call intfac_limiting(ucons(:,:,ie), theta, theta(iamax))

      !--- consistent limiting of primitives
      uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = ucons(2,1:nummat,ie) &
        * uprim(1,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)/ucons(1,1:nummat,ie)
      uprim(2,mmom_idx(nummat,1):mmom_idx(nummat,nummat),ie) = &
        ucons(1,g_mmi%irmin:g_mmi%irmax,ie)*uprim(1,vel_idx(nummat, 0),ie) &
        * ucons(2,1:nummat,ie) /ucons(1,1:nummat,ie)
      uprim(2,vel_idx(nummat, 0),ie) = 0.0
      !---

      !uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
      !uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
      !uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

      !! iv. monotonicity of primitives
      !call superbee_fn(g_nprim, 2.0, 1.0, uneigh(:,1:g_nprim,:), thetap)

      !uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = &
      !  minval(thetap(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
      !  * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      !uprim(2,mmom_idx(nummat,1):mmom_idx(nummat,nummat),ie) = &
      !  minval(thetap(mmom_idx(nummat,1):mmom_idx(nummat,nummat))) &
      !  * uprim(2,mmom_idx(nummat,1):mmom_idx(nummat,nummat),ie)
      !uprim(2,vel_idx(nummat, 0),ie) = thetap(vel_idx(nummat, 0)) &
      !  * uprim(2,vel_idx(nummat, 0),ie)

      !--- common for all equations
      !ucons(2,:,ie) = min(theta(iamax), minval(theta(nummat+1:))) * ucons(2,:,ie)

      !uprim(2,1:nummat,ie) = minval(thetap(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
      !  * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      !uprim(2,vel_idx(nummat, 0),ie) = min(theta(iamax), minval(theta(nummat+1:))) &
      !  * uprim(2,vel_idx(nummat, 0),ie)
      !---

      !--- separate-per-equation
      !ucons(2,1:nummat,ie) = theta(iamax) * ucons(2,1:nummat,ie)
      !thetac = minval(theta(g_mmi%irmin:))
      !ucons(2,g_mmi%irmin:,ie) = thetac * ucons(2,g_mmi%irmin:,ie)

      !uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = &
      !  minval(thetap(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
      !  * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      !uprim(2,vel_idx(nummat, 0),ie) = thetac * uprim(2,vel_idx(nummat, 0),ie)
      !---

      !--- additional limiting step to check TVD constraints
      !uneigh(1:2,g_mmi%irmin:g_mmi%iemax,-1) = ucons(1:2,g_mmi%irmin:g_mmi%iemax,ie-1)
      !uneigh(1:2,g_mmi%irmin:g_mmi%iemax,0)  = ucons(1:2,g_mmi%irmin:g_mmi%iemax,ie)
      !uneigh(1:2,g_mmi%irmin:g_mmi%iemax,1)  = ucons(1:2,g_mmi%irmin:g_mmi%iemax,ie+1)

      !! iv. check if limiting satisfies TVD conditions for conserved quantities
      !tr_cell = .false.
      !do ieqn = g_mmi%irmin,g_mmi%iemax
      !  tr_cell(ieqn) = troubled_cell(2.0, uneigh(:,ieqn,:), &
      !    (coord(ie+1)-coord(ie)))
      !end do !ieqn

      !if (any(tr_cell)) then
      !  call superbee_fn(g_neqns, 2.0, 1.0, uneigh, theta)
      !  theta(g_mmi%irmin:g_mmi%irmax) = minval(theta(g_mmi%irmin:g_mmi%irmax))
      !  theta(g_mmi%iemin:g_mmi%iemax) = minval(theta(g_mmi%iemin:g_mmi%iemax))

      !  do ieqn = g_mmi%irmin,g_mmi%iemax
      !    ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      !  end do !ieqn
      !end if
      !---

    else
    ! 2b. Obtain limiter functions for equation system in single-material cell

      uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
      uneigh(1:2,:,0)  = ucons(1:2,:,ie)
      uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

      ! i. compute limiter function
      call superbee_fn(g_neqns, 2.0, 1.0, uneigh, theta)

      uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
      uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
      uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

      ! ii. monotonicity of primitives
      call superbee_fn(g_nprim, 2.0, 1.0, uneigh(:,1:g_nprim,:), thetap)

      ! iii. use common limiter function for all volume-fractions
      theta(1:nummat) = theta(iamax) !minval(theta(1:nummat))

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

      uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = &
        minval(thetap(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
        * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      uprim(2,mmom_idx(nummat,1):mmom_idx(nummat,nummat),ie) = &
        minval(thetap(mmom_idx(nummat,1):mmom_idx(nummat,nummat))) &
        * uprim(2,mmom_idx(nummat,1):mmom_idx(nummat,nummat),ie)
      uprim(2,vel_idx(nummat, 0),ie) = thetap(vel_idx(nummat, 0)) &
        * uprim(2,vel_idx(nummat, 0),ie)

    end if

  end do !ie

end associate

end subroutine superbee_p1

!-------------------------------------------------------------------------------
!----- system-consistent superbee limiter for P2 dofs:
!-------------------------------------------------------------------------------

subroutine superbee_p2(ucons, uprim)

integer :: ie, ieqn, iamax, imat
real*8  :: dx2, almax, dalmax, theta2(g_neqns), theta1(g_neqns), &
           theta1c, theta2c, &
           thetap2(g_nprim), thetap1(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), alneigh(2,1,-1:1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    !--- 1. P2 derivative limiting

    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / dx2

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,:,0)  = ucons(2:3,:,ie) / dx2

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / dx2

    call superbee_fn(g_neqns, 2.0, 1.0, uneigh, theta2)

    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,1:g_nprim,-1) = uprim(2:3,:,ie-1) / dx2

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,1:g_nprim,0)  = uprim(2:3,:,ie) / dx2

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,1:g_nprim,1)  = uprim(2:3,:,ie+1) / dx2

    ! monotonicity of primitives
    call superbee_fn(g_nprim, 2.0, 1.0, uneigh(:,1:g_nprim,:), thetap2)

    !--- 2. P1 derivative limiting
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call superbee_fn(g_neqns, 2.0, 1.0, uneigh, theta1)

    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    ! monotonicity of primitives
    call superbee_fn(g_nprim, 2.0, 1.0, uneigh(:,1:g_nprim,:), thetap1)
    thetap1(apr_idx(nummat,1):apr_idx(nummat,nummat)) = &
      minval(thetap1(apr_idx(nummat,1):apr_idx(nummat,nummat)))
    thetap2(apr_idx(nummat,1):apr_idx(nummat,nummat)) = &
      minval(thetap2(apr_idx(nummat,1):apr_idx(nummat,nummat)))

    iamax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,iamax,ie)
    dalmax = max( maxval(dabs(ucons(2,1:nummat,ie) + ucons(3,1:nummat,ie))), &
                  maxval(dabs(ucons(2,1:nummat,ie) - ucons(3,1:nummat,ie))) ) /dx2

    ! use common limiter function for all volume-fractions
    theta1(1:nummat) = theta1(iamax) !minval(theta1(1:nummat))
    theta2(1:nummat) = minval(theta2(1:nummat))

    !--- 3. Obtain consistent limiter functions for the equation system
    !       with interface detection
    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      !--- separate-per-equation
      !ucons(2,1:nummat,ie) = max(theta1(iamax), theta2(iamax)) &
      !  * ucons(2,1:nummat,ie)
      !ucons(3,1:nummat,ie) = theta2(iamax) * ucons(3,1:nummat,ie)

      !theta1c = minval(theta1(g_mmi%irmin:))
      !theta2c = minval(theta2(g_mmi%irmin:))
      !ucons(2,g_mmi%irmin:,ie) = max(theta1c, theta2c) &
      !  * ucons(2,g_mmi%irmin:,ie)
      !ucons(3,g_mmi%irmin:,ie) = theta2c * ucons(3,g_mmi%irmin:,ie)

      !uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = &
      !  max(thetap1(apr_idx(nummat,1):apr_idx(nummat,nummat)), &
      !  thetap2(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
      !  * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      !uprim(3,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = &
      !  thetap2(apr_idx(nummat,1):apr_idx(nummat,nummat)) * &
      !  uprim(3,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      !uprim(2,vel_idx(nummat, 0),ie) = max(theta1c, theta2c) &
      !  * uprim(2,vel_idx(nummat, 0),ie)
      !uprim(3,vel_idx(nummat, 0),ie) = theta2c * uprim(3,vel_idx(nummat, 0),ie)
      !---

      !--- special interface treatment
      call intfac_limiting_p2(ucons(:,:,ie), theta1(iamax), theta2(iamax))
      do imat = 1,nummat
        uprim(2:3,apr_idx(nummat, imat),ie) = ucons(2:3,apr_idx(nummat, imat),ie) &
          * uprim(1,apr_idx(nummat, imat),ie)/ucons(1,imat,ie)
      end do !imat
      uprim(2:3,vel_idx(nummat, 0),ie) = 0.0
      !---

    else

      do ieqn = 1,g_neqns
        ucons(3,ieqn,ie) = theta2(ieqn) * ucons(3,ieqn,ie)
        ucons(2,ieqn,ie) = max(theta1(ieqn), theta2(ieqn)) * ucons(2,ieqn,ie)
      end do !ieqn

      uprim(3,:,ie) = thetap2(:) * uprim(3,:,ie)
      uprim(2,:,ie) = max(thetap1(:), thetap2(:)) * uprim(2,:,ie)

    end if

  end do !ie

end associate

end subroutine superbee_p2

!-------------------------------------------------------------------------------
!----- system-consistent weno limiter for P2:
!-------------------------------------------------------------------------------

subroutine weno_p2(ucons, uprim)

integer :: ie, ieqn, iamax, imat
real*8  :: dx2, almax, dalmax, theta2(g_neqns), theta1(g_neqns), &
           thetap(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), alneigh(2,1,-1:1), &
           uxxlim(g_neqns,0:imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  uxxlim = 0.0
  theta2 = 1.0
  theta1 = 1.0

  !--- 1. P2 derivative limiting
  do ie = 1,imax

    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / (dx2*dx2)

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,:,0)  = ucons(2:3,:,ie) / (dx2*dx2)

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / (dx2*dx2)

    call weno_fn(uneigh)
    uxxlim(:,ie) = uneigh(2,:,0) * dx2 * dx2

  end do !ie

  ucons(3,:,:) = uxxlim(:,:)

  uxxlim = 0.0

  !--- 2. P1 derivative limiting
  do ie = 1,imax

    uneigh(1,:,-1) = ucons(1,:,ie-1)
    uneigh(1,:,0)  = ucons(1,:,ie)
    uneigh(1,:,1)  = ucons(1,:,ie+1)

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(2,:,-1) = (ucons(2,:,ie-1) + 2.0*ucons(3,:,ie-1)) / dx2
    uneigh(2,:,0)  = ucons(2,:,ie) / dx2
    uneigh(2,:,1)  = (ucons(2,:,ie+1) - 2.0*ucons(3,:,ie+1)) / dx2

    call weno_fn(uneigh)
    uxxlim(:,ie) = uneigh(2,:,0) * dx2

  end do !ie

  ucons(2,:,:) = uxxlim(:,:)

  do ie = 1,imax

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    iamax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,iamax,ie)
    dalmax = maxval(ucons(2,1:nummat,ie)/dx2)

    !--- 3. Obtain consistent limiter functions for the equation system
    !       with interface detection
    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      call intfac_limiting_p2(ucons(:,:,ie), theta1(iamax), theta2(iamax))

    end if

  end do !ie

end associate

end subroutine weno_p2

!-------------------------------------------------------------------------------
!----- system-consistent superbee+weno limiter for P1P2:
!-------------------------------------------------------------------------------

subroutine superbeeweno_p2(ucons, uprim)

integer :: ie, ieqn, iamax, imat
real*8  :: dx2, almax, dalmax, theta1(g_neqns)
real*8  :: uneigh(2,g_neqns,-1:1), alneigh(2,1,-1:1), &
           uxxlim(g_neqns,0:imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  uxxlim = 0.0

  !--- 1. P2 derivative limiting
  do ie = 1,imax

    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / (dx2*dx2)

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,:,0)  = ucons(2:3,:,ie) / (dx2*dx2)

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / (dx2*dx2)

    call weno_fn(uneigh)
    uxxlim(:,ie) = uneigh(2,:,0) * dx2 * dx2

  end do !ie

  ucons(3,:,:) = uxxlim(:,:)

  do ie = 1,imax

    !--- 2. P1 derivative limiting
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call superbee_fn(g_neqns, 2.0, 1.0, uneigh, theta1)

    iamax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,iamax,ie)
    dalmax = maxval(ucons(2,1:nummat,ie)/dx2)

    ! use common limiter function for all volume-fractions
    theta1(1:nummat) = theta1(iamax)

    !--- 3. Obtain consistent limiter functions for the equation system
    !       with interface detection
    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      call intfac_limiting_p2(ucons(:,:,ie), theta1(iamax), 1.0)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta1(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

    !if ( dabs(sum(ucons(2,1:nummat,ie))) .gt. 1.0d-12 ) &
    !  print*, ie, " : ", sum(ucons(2,1:nummat,ie)), minval(ucons(1,1:nummat,ie))

  end do !ie

end associate

end subroutine superbeeweno_p2

!-------------------------------------------------------------------------------
!----- system-consistent overbee+superbee limiter:
!-------------------------------------------------------------------------------

subroutine oversuperbee(ucons, uprim)

integer :: ie, ieqn, iamax
real*8  :: theta(g_neqns), theta_al, thrho(2)
real*8  :: almax, dalmax, rho1, rho2, vel, rhoe1, rhoe2
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    ! 1. compute limiter function
    call superbee_fn(g_neqns, 2.0, 1.0, uneigh, theta)

    iamax = maxloc(ucons(1,1:nummat,ie), 1)
    theta(1:nummat) = theta(iamax)
    almax = ucons(1,iamax,ie)
    dalmax = maxval(ucons(2,1:nummat,ie)/(0.5 * (coord(ie+1)-coord(ie))))

    ! 2. Obtain consistent limiter functions for the equation system
    !    Interface detection
    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      !   compressive limiting for volume fraction
      call overbee_fn(uneigh(:,iamax,:), theta_al)
      call intfac_limiting(ucons(:,:,ie), theta, theta_al)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end associate

end subroutine oversuperbee

!-------------------------------------------------------------------------------
!----- Superbee limiter for n equations individually
!----- this sub actually calculates the limiter function according to superbee.
!----- the input to this function can be modified to limit P1 or P2 dofs, refer
!----- to superbee_p1 and superbee_p2 respectively
!-------------------------------------------------------------------------------

subroutine superbee_fn(neq,beta_lim,ascale,ucons,theta)

integer, intent(in) :: neq
real*8,  intent(in) :: beta_lim, ascale, ucons(2,neq,-1:1)

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
        phi = min(1.0, (umax-ui)/(2.0*diff))

      else if (diff < -1.0d-16) then
        phi = min(1.0, (umin-ui)/(2.0*diff))

      else
        phi = 1.0

      end if

      ! limiter function
      phi = phi/ascale
      thetal = max( 0.0, max( min(beta_lim*phi, 1.0), min(phi, beta_lim) ) )
      theta(ieqn) = min(thetal, theta(ieqn))

    end do !ifc

  end do !ieqn

end subroutine superbee_fn

!-------------------------------------------------------------------------------
!----- WENO limiter function
!-------------------------------------------------------------------------------

subroutine weno_fn(ucons)

integer :: ieqn, is, nsten, imat, min_cw
real*8  :: wi,epsweno,wenocp1,wt, &
           weight(3), al_weight(g_mmi%nummat,3), &
           gradv(3),osc(3),gradu, &
           ulim(2,g_neqns), &
           ucons(2,g_neqns,-1:1)

associate (nummat=>g_mmi%nummat)

  epsweno = 1.d-10
  wenocp1 = 200.0
  nsten = 3

  do ieqn = 1,g_neqns

    !--- the nsten stencils
    gradv(1) = ucons(2,ieqn,0)
    gradv(2) = ucons(2,ieqn,-1)
    gradv(3) = ucons(2,ieqn,1)

    !--- oscillation indicators
    do is = 1,nsten
      osc(is) = dsqrt(gradv(is)*gradv(is))
    end do !is

    !--- weight functions
    wt = 0.d0

    do is = 1,nsten
      if (is.eq.1) then
        wi = wenocp1
      else
        wi = 1.d0
      end if !is
      weight(is) = wi* (epsweno + osc(is))**(-4.d0)
      wt = wt + weight(is)
    end do !is

    !--- normalize the weight function
    do is = 1,nsten
      weight(is) = weight(is)/wt
    end do !is

    !--- reconstruct the limited gradient
    gradu = 0.d0

    do is = 1,nsten
      gradu = gradu + weight(is)*gradv(is)
    end do !is

    if (ieqn .lt. g_mmi%irmin) al_weight(ieqn,:) = weight(:)

    ulim(2,ieqn) = gradu

  end do !ieqn

  min_cw = minloc(al_weight(:,1), 1)
  al_weight(:,1) = al_weight(min_cw,1)
  al_weight(:,2) = al_weight(min_cw,2)
  al_weight(:,3) = al_weight(min_cw,3)

  do imat = 1,nummat
    ucons(2,imat,0) =   al_weight(imat,1) * ucons(2,imat,0) &
                      + al_weight(imat,2) * ucons(2,imat,-1) &
                      + al_weight(imat,3) * ucons(2,imat,1)
  end do !imat

  ucons(2,g_mmi%irmin:g_neqns,0) = ulim(2,g_mmi%irmin:g_neqns)

end associate

end subroutine weno_fn

!-------------------------------------------------------------------------------
!----- Overbee limiter for n equations individually
!----- this sub actually calculates the limiter function according to overbee
!-------------------------------------------------------------------------------

subroutine overbee_fn(ucons,theta)

real*8,  intent(in) :: ucons(g_tdof,1,-1:1)

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
!----- Function that consistently applies limiter for near-interface cell
!-------------------------------------------------------------------------------

subroutine intfac_limiting(ucons, theta, theta_al)

real*8,  intent(in) :: theta(g_neqns), theta_al

integer :: imat
real*8, dimension(g_mmi%nummat) :: rhom, rhoem
real*8  :: ucons(g_tdof,g_neqns,1), alm, vel

associate (nummat=>g_mmi%nummat)

  !        get primitive variables
  do imat = 1,nummat
    alm         = ucons(1,imat,1)
    rhom(imat)  = ucons(1,g_mmi%irmin+imat-1,1)/alm
    rhoem(imat) = ucons(1,g_mmi%iemin+imat-1,1)/alm
  end do !imat
  vel  = ucons(1,g_mmi%imome,1)/sum( ucons(1,g_mmi%irmin:g_mmi%irmax,1) )

  do imat = 1,nummat
  !        Volume fraction: Keep limiter function the same
    ucons(2,imat,1) = theta_al * ucons(2,imat,1)
  !        Continuity:
    ucons(2,g_mmi%irmin+imat-1,1) = rhom(imat) * ucons(2,imat,1)! &
  !        Energy:
    ucons(2,g_mmi%iemin+imat-1,1) = rhoem(imat) * ucons(2,imat,1)! &
  end do !imat

  !        Momentum:
  ucons(2,g_mmi%imome,1) = vel * sum( ucons(2,g_mmi%irmin:g_mmi%irmax,1) )

end associate

end subroutine intfac_limiting

!-------------------------------------------------------------------------------
!----- Function that consistently applies limiter to P2 dofs for
!----- near-interface cell
!-------------------------------------------------------------------------------

subroutine intfac_limiting_p2(ucons, theta1_al, theta2_al)

real*8,  intent(in) :: theta1_al, theta2_al

integer :: imat
real*8, dimension(g_mmi%nummat) :: rhom, rhoem
real*8  :: ucons(g_tdof,g_neqns,1), alm, vel

associate (nummat=>g_mmi%nummat)

  !        get primitive variables
  do imat = 1,nummat
    alm         = ucons(1,imat,1)
    rhom(imat)  = ucons(1,g_mmi%irmin+imat-1,1)/alm
    rhoem(imat) = ucons(1,g_mmi%iemin+imat-1,1)/alm
  end do !imat
  vel  = ucons(1,g_mmi%imome,1)/sum( ucons(1,g_mmi%irmin:g_mmi%irmax,1) )

  do imat = 1,nummat
  !        Volume fraction: Keep limiter function the same
    ucons(2,imat,1) = max(theta1_al, theta2_al) * ucons(2,imat,1)
    ucons(3,imat,1) = theta2_al * ucons(3,imat,1)
  !        Continuity:
    ucons(2:3,g_mmi%irmin+imat-1,1) = rhom(imat) * ucons(2:3,imat,1)
  !        Energy:
    ucons(2:3,g_mmi%iemin+imat-1,1) = rhoem(imat) * ucons(2:3,imat,1)
  end do !imat

  !        Momentum:
  ucons(2,g_mmi%imome,1) = vel * sum( ucons(2,g_mmi%irmin:g_mmi%irmax,1) )
  ucons(3,g_mmi%imome,1) = vel * sum( ucons(3,g_mmi%irmin:g_mmi%irmax,1) )

end associate

end subroutine intfac_limiting_p2

!-------------------------------------------------------------------------------
!----- Interface-cell / mixed-cell indicator based on vol-frac gradient
!-------------------------------------------------------------------------------

logical function interface_cell(alcell, dalcell)

real*8, intent(in) :: alcell, dalcell
real*8  :: scale_al

logical :: al, dal

  scale_al = 1.0d-2
  dal = dabs(dalcell) .gt. scale_al

  scale_al = 1.0d-4 !10000.0*g_alphamin
  al = ( (alcell .gt. scale_al) .and. (alcell .lt. 1.0-scale_al) )

  if (al .or. dal) then
    interface_cell = .true.

  else
    interface_cell = .false.

  end if

end function interface_cell

!-------------------------------------------------------------------------------
!----- Troubled-cell indicator
!-------------------------------------------------------------------------------

logical function troubled_cell(beta, u, dx)

real*8, intent(in) :: beta, u(g_tdof,1,-1:1), dx

integer :: ifc
real*8  :: dplus, dminu, difc, a1, a2, a3, mm

  troubled_cell = .false.

  dplus = u(1,1,0)-u(1,1,-1)
  dminu = u(1,1,1)-u(1,1,0)

  do ifc = 1,2
    difc = ((-1.0)**ifc) * u(2,1,0)

    a1 = difc
    a2 = dplus
    a3 = dminu

    ! modified minmod
    if (dabs(a1) < beta*dx*dx) then
      mm = a1

    else
      if (a1*a2 > 0.0 .and. a2*a3 > 0.0) then
        mm = dsign(1.0, a1) * min(dabs(a1), dabs(a2), dabs(a3))
      else
        mm = 0.0
      end if

    end if

    if (dabs(mm-difc) > 1.0e-14) then
      troubled_cell = troubled_cell .or. .true.
    else
      troubled_cell = troubled_cell .or. .false.
    end if
  end do !ifc

end function troubled_cell

!-------------------------------------------------------------------------------

END MODULE reconstruction_mm6eq
