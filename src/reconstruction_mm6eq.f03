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

subroutine reconstruction_p0(ucons, uprim, ndof_el)

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

end subroutine reconstruction_p0

!-------------------------------------------------------------------------------
!----- P0P1 reconstruction:
!-------------------------------------------------------------------------------

subroutine reconstruction_p0p1(ucons, uprim, ndof_el)

integer, intent(in) :: ndof_el(2,0:imax+1)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  !--- central difference reconstruction (least-squares for uniform meshes)
  do ie = 1,imax
    ucons(2,:,ie) = 0.25 * (ucons(1,:,ie+1) - ucons(1,:,ie-1))
  end do !ie

  !--- central difference reconstruction for primitive quantities
  do ie = 1,imax
    uprim(2,:,ie) = 0.25 * (uprim(1,:,ie+1) - uprim(1,:,ie-1))
  end do !ie

  !--- limit second-order solution
  call limiting_p1(ucons, uprim)

end subroutine reconstruction_p0p1

!-------------------------------------------------------------------------------
!----- P1 reconstruction (only limiting):
!-------------------------------------------------------------------------------

subroutine reconstruction_p1(ucons, uprim, ndof_el)

integer, intent(in) :: ndof_el(2,0:imax+1)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  !--- if interface reconstruction is active, check if any cells are FV2
  if (g_intreco > 0) then
    ! central difference reconstruction (least-squares for uniform meshes)
    ! for FV2 cells
    do ie = 1,imax
      if (ndof_el(1,ie) == 1) then
        ucons(2,:,ie) = 0.25 * (ucons(1,:,ie+1) - ucons(1,:,ie-1))
        uprim(2,:,ie) = 0.25 * (uprim(1,:,ie+1) - uprim(1,:,ie-1))
      end if
    end do !ie
  end if

  !--- limit second-order solution
  call limiting_p1(ucons, uprim)

end subroutine reconstruction_p1

!-------------------------------------------------------------------------------
!----- P1P2 reconstruction:
!-------------------------------------------------------------------------------

subroutine reconstruction_p1p2(ucons, uprim, ndof_el)

integer, intent(in) :: ndof_el(2,0:imax+1)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  !--- if interface reconstruction is active, check if any cells are FV2
  if (g_intreco > 0) then
    ! central difference reconstruction (least-squares for uniform meshes)
    ! for FV2 cells
    do ie = 1,imax
      if (ndof_el(1,ie) == 1) then
        ucons(2,:,ie) = 0.25 * (ucons(1,:,ie+1) - ucons(1,:,ie-1))
        uprim(2,:,ie) = 0.25 * (uprim(1,:,ie+1) - uprim(1,:,ie-1))
      end if
    end do !ie
  end if

  !--- reconstruct third-order conserved quantities
  call leastsquares_p1p2(g_neqns, ucons)

  !--- reconstruct third-order primitive quantities
  call leastsquares_p1p2(g_nprim, uprim)

  !--- limit third-order solution
  call limiting_p2(ucons, uprim)

end subroutine reconstruction_p1p2

!-------------------------------------------------------------------------------
!----- high-order reconstruction given the basis functions at required point
!-------------------------------------------------------------------------------

subroutine ho_reconstruction(neq, udof, basis, uho)

integer, intent(in) :: neq
real*8, intent(in) :: udof(g_tdof,neq), basis(g_tdof)

real*8, intent(out) :: uho(neq)

integer :: ieqn

  uho = 0.0

  !--- dgp0
  if (g_nsdiscr .eq. 0) then
    do ieqn = 1,neq
      uho(ieqn) = udof(1,ieqn)
    end do !ieqn

  !--- rdgp0p1 or dgp1
  elseif ((g_nsdiscr .eq. 11) .or. (g_nsdiscr .eq. 1)) then
    do ieqn = 1,neq
      uho(ieqn) = udof(1,ieqn) + basis(2)*udof(2,ieqn)
    end do !ieqn

  !--- rdgp1p2
  elseif (g_nsdiscr .eq. 12) then
    do ieqn = 1,neq
      uho(ieqn) = udof(1,ieqn) + basis(2)*udof(2,ieqn) &
        + basis(3)*udof(3,ieqn)
    end do !ieqn

  end if

end subroutine ho_reconstruction

!-------------------------------------------------------------------------------
!----- Least-squares reconstruction for P1P2:
!-------------------------------------------------------------------------------

subroutine leastsquares_p1p2(neq, usol)

integer, intent(in) :: neq

integer :: ig, j, ie, je, ieqn, ngauss
data       ngauss/2/

real*8  :: dxi, dxj, xci, xcj, xg, wi, &
           carea(2), weight(2), &
           b2, b3, b2t, b3t, &
           r1, r2, rhs, lhs
real*8  :: usol(g_tdof,neq,0:imax+1)

  call rutope(1, ngauss, carea, weight)

  do ieqn = 1,neq

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
        r1  = usol(1,ieqn,je) - ( usol(1,ieqn,ie) + usol(2,ieqn,ie)*b2t )
        r2  = usol(2,ieqn,je)*dxi/dxj - usol(2,ieqn,ie)
        rhs = rhs + b3t*(r1) + b2*(r2)

      end if
    end do !j

    ! get p2 dof
    usol(3,ieqn,ie) = rhs/lhs

  end do !ie

  !ie = 1
  !dxi = 0.5 * (coord(ie+1)-coord(ie))
  !dxj = 0.5 * (coord(ie+2)-coord(ie+1))
  !usol(3,ieqn,ie) = (dxi*dxi)*(usol(2,ieqn,ie+1)/dxj - usol(2,ieqn,ie)/dxi) / (dxj+dxi)

  !ie = imax
  !dxi = 0.5 * (coord(ie+1)-coord(ie))
  !dxj = 0.5 * (coord(ie)-coord(ie-1))
  !usol(3,ieqn,ie) = (dxi*dxi)*(usol(2,ieqn,ie)/dxi - usol(2,ieqn,ie-1)/dxj) / (dxi+dxj)

  end do !ieqn

end subroutine leastsquares_p1p2

!-------------------------------------------------------------------------------
!----- Reconstructing primitives from DG solution:
!-------------------------------------------------------------------------------

subroutine weak_recons_primitives(ucons, uprim, ndof_el)

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: ig, ie, ieqn, imat, ngauss

real*8  :: dxi, xci, xg, wi, &
           carea(3), weight(3), &
           b2, b3, rhs(g_gdof,g_nprim), lhs(g_gdof)

real*8  :: rhomat, rhoemat, upg(g_nprim), ug(g_neqns), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  if (g_pureco == 1) then

    ngauss = get_numqpoints(g_nsdiscr)
    call rutope(1, ngauss, carea, weight)

    do ie = 0,imax+1

      dxi = coord(ie+1)-coord(ie)
      xci = 0.5*(coord(ie+1)+coord(ie))

      ! get lhs
      lhs(1) = dxi
      if (g_nsdiscr .gt. 1) lhs(2) = lhs(1)/3.0
      !if (g_nsdiscr .gt. 11) lhs(3) = lhs(1)/45.0

      ! quadrature
      rhs = 0.0
      do ig = 1,ngauss
        xg = carea(ig) * 0.5*dxi + xci
        wi = 0.5 * weight(ig)
        b2 = p1basis(xg, xci, dxi)
        !b3 = p2basis(xg, xci, dxi)

        !--- dgp0 or rdgp0p1
        if ((g_nsdiscr.eq.0) .or. (g_nsdiscr.eq.1)) then
          do ieqn = 1,g_neqns
            ug(ieqn) = ucons(1,ieqn,ie)
          end do !ieqn

        !--- dgp1 or rdgp1p2
        elseif ((g_nsdiscr .ge. 11)) then
          do ieqn = 1,g_neqns
            ug(ieqn) = ucons(1,ieqn,ie) + b2*ucons(2,ieqn,ie)
          end do !ieqn

        !!--- rdgp1p2
        !elseif (g_nsdiscr .eq. 12) then
        !  do ieqn = 1,g_neqns
        !    ug(ieqn) = ucons(1,ieqn,ie) + b2*ucons(2,ieqn,ie) &
        !      + b3*ucons(3,ieqn,ie)
        !  end do !ieqn

        end if

        ! primitives at quadrature point
        upg = 0.0
        upg(vel_idx(nummat, 0)) = ug(g_mmi%imome)/sum(ug(g_mmi%irmin:g_mmi%irmax))
        do imat = 1,nummat
          rhomat = ug(g_mmi%irmin+imat-1)
          rhoemat = ug(g_mmi%iemin+imat-1)
          upg(apr_idx(nummat, imat)) = &
            eos3_alphapr(g_gam(imat), g_pc(imat), ug(imat), rhomat, rhoemat, &
              upg(vel_idx(nummat, 0)))
        end do !imat

        ! get rhs
        do ieqn = 1,g_nprim
          rhs(1,ieqn) = rhs(1,ieqn) + wi*dxi*upg(ieqn)
          if (g_nsdiscr .gt. 1) rhs(2,ieqn) = rhs(2,ieqn) + wi*dxi*upg(ieqn)*b2
          !if (g_nsdiscr .gt. 11) rhs(3,ieqn) = rhs(3,ieqn) + wi*dxi*upg(ieqn)*b3
        end do !ieqn

      end do !ig

      ! get primitive variable dofs
      do ieqn = 1,g_nprim
        uprim(1,ieqn,ie) = rhs(1,ieqn)/lhs(1)
        if (g_nsdiscr .gt. 1) uprim(2,ieqn,ie) = rhs(2,ieqn)/lhs(2)
        !if (g_nsdiscr .gt. 11) uprim(3,ieqn,ie) = rhs(3,ieqn)/lhs(3)
      end do !ieqn

    end do !ie

    if (g_nsdiscr .eq. 1) then
      uprim(2,:,0) = 0.0
      uprim(2,:,imax+1) = 0.0

    else if (g_nsdiscr .eq. 12) then
      uprim(3,:,0) = 0.0
      uprim(3,:,imax+1) = 0.0

    end if

  end if

end associate

end subroutine weak_recons_primitives

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
    call min_vertexbased(ucons, uprim)

  case(2)
    call vertexbased_p1(ucons, uprim)

  case(3)
    call oververtexbased(ucons, uprim)

  case(4)
    call weno_p1(ucons, uprim)

  case(6)
    call thincvertexbased_p1(ucons, uprim)

  case default
    write(*,*) "Error: incorrect p1-limiter index in control file: ", g_nlim
    call exit

  end select

  if (g_intreco == 0) call boundpreserve_alpha_p1(ucons)

end subroutine limiting_p1

!-------------------------------------------------------------------------------
!----- P2 limiting:
!-------------------------------------------------------------------------------

subroutine limiting_p2(ucons, uprim)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  select case (g_nlim)

  case(0)

  case(2)
    call vertexbased_p2(ucons, uprim)

  case(4)
    call weno_p2(ucons, uprim)

  case(5)
    call vertexbasedweno_p2(ucons, uprim)

  case(6)
    call thincvertexbased_p2(ucons, uprim)

  case(7)
    call thincvertexbasedweno_p2(ucons, uprim)

  case default
    write(*,*) "Error: incorrect p2-limiter index in control file: ", g_nlim
    call exit

  end select

  if (g_intreco == 0) call boundpreserve_alpha_p2(ucons)

end subroutine limiting_p2

!-------------------------------------------------------------------------------
!----- Bound preserving limiter for p1:
!-------------------------------------------------------------------------------

subroutine boundpreserve_alpha_p1(ucons)

integer :: ie, ig, ngauss, imat, mmax
real*8  :: careap(2), &
           xc, dx, xg, basis(g_tdof), &
           alm(1), thal(g_mmi%nummat), thal1, thal2, &
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
      xc = 0.5*(coord(ie+1)+coord(ie))
      dx = coord(ie+1)-coord(ie)
      xg = careap(ig) * 0.5*dx + xc
      call get_basisfns(xg, xc, dx, basis)
      call ho_reconstruction(1, ucons(:,imat,ie), basis, alm)
      !arhom = ucons(1,g_mmi%irmin+imat-1,ie) + careap(ig) * ucons(2,g_mmi%irmin+imat-1,ie)

      alm_min = min(alm(1), alm_min)
      alm_max = max(alm(1), alm_max)
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

  mmax = maxloc(ucons(1,1:nummat,ie), 1)
  thal(1:nummat) = minval(thal)

  ucons(2,1:nummat,ie) = thal(1:nummat) * ucons(2,1:nummat,ie)

  if (g_nmatint .eq. 1) then
    do imat = 1,nummat
      ucons(2,g_mmi%irmin+imat-1,ie) = thal(imat) * ucons(2,g_mmi%irmin+imat-1,ie)
      ucons(2,g_mmi%iemin+imat-1,ie) = thal(imat) * ucons(2,g_mmi%iemin+imat-1,ie)
    end do !imat
  end if

  end do !ie

end associate

end subroutine boundpreserve_alpha_p1

!-------------------------------------------------------------------------------
!----- Bound preserving limiter for p2:
!-------------------------------------------------------------------------------

subroutine boundpreserve_alpha_p2(ucons)

integer :: ie, ig, ngauss, imat, mmax
real*8  :: carea(2), weight(2), careap(4), &
           xc, dx, xg, basis(g_tdof), &
           alm(1), thal(g_mmi%nummat), thal1, thal2, &
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

      ! reconstructed volume fraction
      xc = 0.5*(coord(ie+1)+coord(ie))
      dx = coord(ie+1)-coord(ie)
      xg = careap(ig) * 0.5*dx + xc
      call get_basisfns(xg, xc, dx, basis)
      call ho_reconstruction(1, ucons(:,imat,ie), basis, alm)
      !arhom = ucons(1,g_mmi%irmin+imat-1,ie) &
      !        + careap(ig) * ucons(2,g_mmi%irmin+imat-1,ie) &
      !        + b3 * ucons(3,g_mmi%irmin+imat-1,ie)

      alm_min = min(alm(1), alm_min)
      alm_max = max(alm(1), alm_max)
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

  mmax = maxloc(ucons(1,1:nummat,ie), 1)
  thal(1:nummat) = minval(thal)

  do imat = 1,nummat
    ucons(2:3,imat,ie) = thal(imat) * ucons(2:3,imat,ie)
  end do !imat

  if (g_nmatint .eq. 1) then
    do imat = 1,nummat
      ucons(2:3,g_mmi%irmin+imat-1,ie) = thal(imat) * ucons(2:3,g_mmi%irmin+imat-1,ie)
      ucons(2:3,g_mmi%iemin+imat-1,ie) = thal(imat) * ucons(2:3,g_mmi%iemin+imat-1,ie)
    end do !imat
  end if

  end do !ie

end associate

end subroutine boundpreserve_alpha_p2

!-------------------------------------------------------------------------------
!----- vertexbased limiter:
!-------------------------------------------------------------------------------

subroutine min_vertexbased(ucons, uprim)

integer :: ie, mmax
real*8  :: theta(g_neqns), thetap(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    ! 1. compute limiter function
    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta)

    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    theta(1:nummat) = theta(mmax)

    ! min for all equations
    theta = minval(theta)

    ! primitive quantities
    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap)

    ! 2. limit 2nd dofs
    ucons(2,:,ie) = theta(:) * ucons(2,:,ie)

    uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = &
      minval(thetap(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
      * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
    uprim(2,vel_idx(nummat, 0),ie) = thetap(vel_idx(nummat, 0)) &
      * uprim(2,vel_idx(nummat, 0),ie)

  end do !ie

end associate

end subroutine min_vertexbased

!-------------------------------------------------------------------------------
!----- system-consistent vertexbased limiter for P1 dofs:
!-------------------------------------------------------------------------------

subroutine vertexbased_p1(ucons, uprim)

integer :: ie, ieqn, mmax, imat
real*8  :: almax, dalmax, theta(g_neqns), thetap(g_nprim), thetac
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    !--- 1. detect interface/single-material cell
    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    dalmax = maxval(dabs(ucons(2,1:nummat,ie))/(0.5 * (coord(ie+1)-coord(ie))))

    !--- 2. obtain limiter function for individual unknowns
    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta)

    ! primitive quantities
    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap)

    ! use common limiter function for all volume-fractions
    theta(1:nummat) = theta(mmax) !minval(theta(1:nummat))

    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then
    !--- 3a. Obtain consistent limiter functions for equation system at interface

      ! consistent limiting
      call intfac_limiting(ucons(:,:,ie), theta, theta(mmax))

      !! consistent limiting of primitives
      !uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = 0.0
      !uprim(2,vel_idx(nummat, 0),ie) = 0.0

      ! monotonicity of primitives
      uprim(2,:,ie) = thetap(:) * uprim(2,:,ie)

      !--- common for all equations
      !thetac = minval(theta)
      !ucons(2,:,ie) = thetac * ucons(2,:,ie)

      !uprim(2,1:nummat,ie) = minval(thetap(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
      !  * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      !uprim(2,vel_idx(nummat, 0),ie) = min(theta(mmax), minval(theta(nummat+1:))) &
      !  * uprim(2,vel_idx(nummat, 0),ie)
      !---

      !--- separate-per-equation
      !ucons(2,1:nummat,ie) = theta(mmax) * ucons(2,1:nummat,ie)
      !thetac = minval(theta(g_mmi%irmin:g_mmi%irmax))
      !ucons(2,g_mmi%irmin:g_mmi%irmax,ie) = thetac &
      !  * ucons(2,g_mmi%irmin:g_mmi%irmax,ie)
      !thetac = minval(theta(g_mmi%iemin:g_mmi%iemax))
      !ucons(2,g_mmi%iemin:g_mmi%iemax,ie) = thetac &
      !  * ucons(2,g_mmi%iemin:g_mmi%iemax,ie)
      !thetac = theta(g_mmi%imome)
      !ucons(2,g_mmi%imome,ie) = thetac * ucons(2,g_mmi%imome,ie)

      !uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie) = &
      !  minval(thetap(apr_idx(nummat,1):apr_idx(nummat,nummat))) &
      !  * uprim(2,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
      !uprim(2,vel_idx(nummat, 0),ie) = thetap(vel_idx(nummat, 0)) &
      !  * uprim(2,vel_idx(nummat, 0),ie)
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
      !  call vertexbased_fn(g_neqns, uneigh, theta)
      !  theta(g_mmi%irmin:g_mmi%irmax) = minval(theta(g_mmi%irmin:g_mmi%irmax))
      !  theta(g_mmi%iemin:g_mmi%iemax) = minval(theta(g_mmi%iemin:g_mmi%iemax))

      !  do ieqn = g_mmi%irmin,g_mmi%iemax
      !    ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      !  end do !ieqn
      !end if
      !---

    else
    !--- 3b. Obtain limiter functions for equation system in single-material cell

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

      uprim(2,:,ie) = thetap(:) * uprim(2,:,ie)

    end if

  end do !ie

end associate

end subroutine vertexbased_p1

!-------------------------------------------------------------------------------
!----- system-consistent vertexbased limiter for P1 dofs:
!-------------------------------------------------------------------------------

subroutine thincvertexbased_p1(ucons, uprim)

integer :: ie, ieqn, mmax, imat, matint(g_mmi%nummat)
real*8  :: almax, theta(g_neqns), thetap(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)
logical :: is_intcell

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    !--- 1. obtain limiter function for individual unknowns
    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta)

    ! primitive quantities
    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap)

    !--- 2. detect interface/single-material cell
    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)

    is_intcell = intrecons_cell(ucons(1,1:nummat,ie), matint)

    !--- 3. determine limiting appropriately
    if (.not.is_intcell) then

      ! use common limiter function for all volume-fractions
      theta(1:nummat) = theta(mmax) !minval(theta(1:nummat))

    else

      ! no limiter for volume fractions in interface cells
      do imat = 1,nummat
        if (matint(imat) == 1) theta(imat) = 1.0
      end do !imat

    end if

    !--- 4. apply limiting
    do ieqn = 1,g_neqns
      ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
    end do !ieqn

    uprim(2,:,ie) = thetap(:) * uprim(2,:,ie)

  end do !ie

end associate

end subroutine thincvertexbased_p1

!-------------------------------------------------------------------------------
!----- system-consistent weno limiter for P1:
!-------------------------------------------------------------------------------

subroutine weno_p1(ucons, uprim)

integer :: ie, ieqn, mmax!, imat
real*8  :: dx2, almax, dalmax, theta(g_neqns)!, thetap(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), &
           uxlim(g_neqns,0:imax+1), &
           pxlim(g_nprim,0:imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)
logical :: trcell

associate (nummat=>g_mmi%nummat)

  uxlim = 0.0
  pxlim = 0.0

  !--- 1. WENO reconstructed dofs
  do ie = 1,imax

    dx2 = 0.5 * (coord(ie+1)-coord(ie))

    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    trcell = .false.
    if ((ie.eq.1) .or. (ie.eq.imax)) then
      trcell = .true.
    else
      do ieqn = g_mmi%irmin, g_mmi%irmax
        trcell = trcell &
          .or. troubled_cell(2.0, uneigh(:,ieqn,:), coord(ie+1)-coord(ie))
      end do !ieqn
    end if

    if (trcell) then
      if (ie.eq.1) then
        uneigh(2,:,-1) = 0.0
      else
        uneigh(2,:,-1) = (ucons(1,:,ie)-ucons(1,:,ie-2))/(4.0*dx2)
      end if

      if (ie.eq.imax) then
        uneigh(2,:,1) = 0.0
      else
        uneigh(2,:,1) = (ucons(1,:,ie+2)-ucons(1,:,ie))/(4.0*dx2)
      end if

      uneigh(2,:,0) = uneigh(2,:,0)/dx2
      call weno_fn(g_neqns, 100.0, uneigh)
      uxlim(:,ie) = uneigh(2,:,0) * dx2

    else
      uxlim(:,ie) = ucons(2,:,ie)

    end if

    ! primitive quantities
    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    trcell = .false.
    if ((ie.eq.1) .or. (ie.eq.imax)) then
      trcell = .true.
    else
      do ieqn = 1,g_nprim
        trcell = trcell &
          .or. troubled_cell(2.0, uneigh(:,ieqn,:), coord(ie+1)-coord(ie))
      end do !ieqn
    end if

    if (trcell) then
      if (ie.eq.1) then
        uneigh(2,1:g_nprim,-1) = 0.0
      else
        uneigh(2,1:g_nprim,-1) = (uprim(1,:,ie)-uprim(1,:,ie-2))/(4.0*dx2)
      end if

      if (ie.eq.imax) then
        uneigh(2,1:g_nprim,1) = 0.0
      else
        uneigh(2,1:g_nprim,1) = (uprim(1,:,ie+2)-uprim(1,:,ie))/(4.0*dx2)
      end if

      uneigh(2,1:g_nprim,0) = uneigh(2,1:g_nprim,0)/dx2
      call weno_fn(g_nprim, 100.0, uneigh(:,1:g_nprim,:))
      pxlim(:,ie) = uneigh(2,1:g_nprim,0) * dx2

    else
      pxlim(:,ie) = uprim(2,:,ie)

    end if

  end do !ie

  !--- 2. replace dofs with WENO dofs
  ucons(2,:,:) = uxlim(:,:)
  uprim(2,:,:) = pxlim(:,:)

  !--- 3. Obtain consistent limiter at material interface cells
  do ie = 1,imax

    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    dalmax = maxval(ucons(2,1:nummat,ie)/dx2)

    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      theta = 1.0
      call intfac_limiting(ucons(:,:,ie), theta, 1.0)

    end if

  end do !ie

end associate

end subroutine weno_p1

!-------------------------------------------------------------------------------
!----- system-consistent vertexbased limiter for P2 dofs:
!-------------------------------------------------------------------------------

subroutine vertexbased_p2(ucons, uprim)

integer :: ie, ieqn, mmax, imat
real*8  :: dx2, almax, dalmax, theta2(g_neqns), theta1(g_neqns), &
           theta1c, theta2c, &
           thetap2(g_nprim), thetap1(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), alneigh(2,1,-1:1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)
logical :: is_weno

associate (nummat=>g_mmi%nummat)

  is_weno = .false.

  do ie = 1,imax

    dx2 = 0.5 * (coord(ie+1)-coord(ie))

    !--- 1. detect interface/single-material cell
    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    dalmax = max( maxval(dabs(ucons(2,1:nummat,ie) + ucons(3,1:nummat,ie))), &
                  maxval(dabs(ucons(2,1:nummat,ie) - ucons(3,1:nummat,ie))) ) /dx2

    !--- 2. obtain limiter function for individual unknowns
    ! i. P2 derivative limiting

    ! conserved quantities
    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / dx2

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,:,0)  = ucons(2:3,:,ie) / dx2

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / dx2

    call vertexbased_fn(g_neqns, uneigh, theta2)

    ! primitive quantities
    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,1:g_nprim,-1) = uprim(2:3,:,ie-1) / dx2

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,1:g_nprim,0)  = uprim(2:3,:,ie) / dx2

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,1:g_nprim,1)  = uprim(2:3,:,ie+1) / dx2

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap2)

    ! ii. P1 derivative limiting

    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta1)

    ! primitive quantities
    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap1)

    ! use common limiter function for all volume-fractions
    theta1(1:nummat) = theta1(mmax) !minval(theta1(1:nummat))
    theta2(1:nummat) = minval(theta2(1:nummat))

    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then
    !--- 3a. Obtain consistent limiter functions for equation system at interface

      !--- separate-per-equation
      !ucons(2,1:nummat,ie) = max(theta1(mmax), theta2(mmax)) &
      !  * ucons(2,1:nummat,ie)
      !ucons(3,1:nummat,ie) = theta2(mmax) * ucons(3,1:nummat,ie)

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

      call intfac_limiting_p2(ucons(:,:,ie), is_weno, theta1(mmax), &
        theta2(mmax))

      !! consistent limiting of primitives
      !do imat = 1,nummat
      !  uprim(2:3,apr_idx(nummat, imat),ie) = 0.0
      !end do !imat
      !uprim(2:3,vel_idx(nummat, 0),ie) = 0.0

      ! monotonicity of primitives
      do ieqn = 1,g_nprim
        uprim(3,ieqn,ie) = thetap2(ieqn) * uprim(3,ieqn,ie)
        uprim(2,ieqn,ie) = max(thetap1(ieqn), thetap2(ieqn)) * uprim(2,ieqn,ie)
      end do !ieqn

    else
    !--- 3b. Obtain limiter functions for equation system in single-material cell

      do ieqn = 1,g_neqns
        ucons(3,ieqn,ie) = theta2(ieqn) * ucons(3,ieqn,ie)
        ucons(2,ieqn,ie) = max(theta1(ieqn), theta2(ieqn)) * ucons(2,ieqn,ie)
      end do !ieqn

      do ieqn = 1,g_nprim
        uprim(3,ieqn,ie) = thetap2(ieqn) * uprim(3,ieqn,ie)
        uprim(2,ieqn,ie) = max(thetap1(ieqn), thetap2(ieqn)) * uprim(2,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end associate

end subroutine vertexbased_p2

!-------------------------------------------------------------------------------
!----- system-consistent THINC+vertexbased limiter for P2 dofs:
!-------------------------------------------------------------------------------

subroutine thincvertexbased_p2(ucons, uprim)

integer :: ie, ieqn, mmax, imat, matint(g_mmi%nummat)
real*8  :: dx2, almax, theta2(g_neqns), theta1(g_neqns), &
           thetap2(g_nprim), thetap1(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), alneigh(2,1,-1:1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)
logical :: is_intcell

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    dx2 = 0.5 * (coord(ie+1)-coord(ie))

    !--- i. obtain limiter function for individual unknowns

    ! 1. P2 derivative limiting
    ! conserved quantities
    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / dx2

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,:,0)  = ucons(2:3,:,ie) / dx2

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / dx2

    call vertexbased_fn(g_neqns, uneigh, theta2)

    ! primitive quantities
    dx2 = 0.5 * (coord(ie)-coord(ie-1))
    uneigh(1:2,1:g_nprim,-1) = uprim(2:3,:,ie-1) / dx2

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    uneigh(1:2,1:g_nprim,0)  = uprim(2:3,:,ie) / dx2

    dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
    uneigh(1:2,1:g_nprim,1)  = uprim(2:3,:,ie+1) / dx2

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap2)

    ! 2. P1 derivative limiting
    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta1)

    ! primitive quantities
    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap1)

    !--- ii. detect interface/single-material cell
    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    is_intcell = intrecons_cell(ucons(1,1:nummat,ie), matint)

    !--- iii. determine limiting appropriately
    if (.not.is_intcell) then

      ! use common limiter function for all volume-fractions
      theta1(1:nummat) = theta1(mmax) !minval(theta1(1:nummat))
      theta2(1:nummat) = minval(theta2(1:nummat))

    else

      ! no limiter for volume fractions in interface cells
      do imat = 1,nummat
        if (matint(imat) == 1) then
          theta1(imat) = 1.0
          theta2(imat) = 1.0
        end if
      end do !imat

    end if

    !--- iv. apply limiting
    do ieqn = 1,g_neqns
      ucons(3,ieqn,ie) = theta2(ieqn) * ucons(3,ieqn,ie)
      ucons(2,ieqn,ie) = max(theta1(ieqn), theta2(ieqn)) * ucons(2,ieqn,ie)
    end do !ieqn

    do ieqn = 1,g_nprim
      uprim(3,ieqn,ie) = thetap2(ieqn) * uprim(3,ieqn,ie)
      uprim(2,ieqn,ie) = max(thetap1(ieqn), thetap2(ieqn)) * uprim(2,ieqn,ie)
    end do !ieqn

  end do !ie

end associate

end subroutine thincvertexbased_p2

!-------------------------------------------------------------------------------
!----- system-consistent weno limiter for P2:
!-------------------------------------------------------------------------------

subroutine weno_p2(ucons, uprim)

integer :: ie, ieqn, mmax, imat
real*8  :: dx2, almax, dalmax, theta2(g_neqns), theta1(g_neqns), &
           thetap(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), alneigh(2,1,-1:1), &
           uxxlim(g_neqns,0:imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)
logical :: is_weno

associate (nummat=>g_mmi%nummat)

  is_weno = .true.
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

    call weno_fn(g_neqns, 100.0, uneigh)
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

    call weno_fn(g_neqns, 100.0, uneigh)
    uxxlim(:,ie) = uneigh(2,:,0) * dx2

  end do !ie

  ucons(2,:,:) = uxxlim(:,:)

  do ie = 1,imax

    dx2 = 0.5 * (coord(ie+1)-coord(ie))
    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    dalmax = maxval(ucons(2,1:nummat,ie)/dx2)

    !--- 3. Obtain consistent limiter functions for the equation system
    !       with interface detection
    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      call intfac_limiting_p2(ucons(:,:,ie), is_weno, theta1(mmax), &
        theta2(mmax))

    end if

  end do !ie

end associate

end subroutine weno_p2

!-------------------------------------------------------------------------------
!----- system-consistent vertexbased+weno limiter for P1P2:
!-------------------------------------------------------------------------------

subroutine vertexbasedweno_p2(ucons, uprim)

integer :: ie, ieqn, mmax, imat
real*8  :: dx2, almax, dalmax, theta1(g_neqns), thetap1(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), upneigh(2,g_nprim,-1:1), &
           ueq(2,-1:1), &
           uxxlim(g_neqns,0:imax+1), &
           pxxlim(g_nprim,0:imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)
logical :: trcell, is_weno

associate (nummat=>g_mmi%nummat)

  is_weno = .true.
  uxxlim = 0.0
  pxxlim = 0.0

  !--- 1. P2 derivative limiting
  do ie = 1,imax

    trcell = .false.
    if ((ie.eq.1) .or. (ie.eq.imax)) then
      trcell = .true.
    else
      do ieqn = g_mmi%irmin, g_mmi%irmax
        ueq = ucons(1:2,ieqn,ie-1:ie+1)
        trcell = trcell &
          .or. troubled_cell(2.0, ueq, coord(ie+1)-coord(ie))
      end do !ieqn
      do ieqn = 1,g_nprim
        ueq = uprim(1:2,ieqn,ie-1:ie+1)
        trcell = trcell &
          .or. troubled_cell(2.0, ueq, coord(ie+1)-coord(ie))
      end do !ieqn
    end if

    if (trcell) then
      dx2 = 0.5 * (coord(ie)-coord(ie-1))
      uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+1)-coord(ie))
      uneigh(1:2,:,0)  = ucons(2:3,:,ie) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
      uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / (dx2*dx2)

      call weno_fn(g_neqns, 200.0, uneigh)
      uxxlim(:,ie) = uneigh(2,:,0) * dx2 * dx2

    else
      uxxlim(:,ie) = ucons(3,:,ie)

    end if

    ! primitive quantities
    if (trcell) then
      dx2 = 0.5 * (coord(ie)-coord(ie-1))
      upneigh(1:2,:,-1) = uprim(2:3,:,ie-1) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+1)-coord(ie))
      upneigh(1:2,:,0)  = uprim(2:3,:,ie) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
      upneigh(1:2,:,1)  = uprim(2:3,:,ie+1) / (dx2*dx2)

      call weno_fn(g_nprim, 200.0, upneigh)
      pxxlim(:,ie) = upneigh(2,:,0) * dx2 * dx2

    else
      pxxlim(:,ie) = uprim(3,:,ie)

    end if

  end do !ie

  ucons(3,:,:) = uxxlim(:,:)
  uprim(3,:,:) = pxxlim(:,:)

  do ie = 1,imax

    !--- 2. P1 derivative limiting
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta1)

    ! primitive quantities
    upneigh(1:2,:,-1) = uprim(1:2,:,ie-1)
    upneigh(1:2,:,0)  = uprim(1:2,:,ie)
    upneigh(1:2,:,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, upneigh, thetap1)

    do ieqn = 1,g_nprim
      uprim(2,ieqn,ie) = thetap1(ieqn) * uprim(2,ieqn,ie)
    end do !ieqn

    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    dalmax = maxval(ucons(2,1:nummat,ie)/dx2)

    ! use common limiter function for all volume-fractions
    theta1(1:nummat) = theta1(mmax)

    !--- 3. Obtain consistent limiter functions for the equation system
    !       with interface detection
    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      call intfac_limiting_p2(ucons(:,:,ie), is_weno, theta1(mmax), 1.0)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta1(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

  end do !ie

end associate

end subroutine vertexbasedweno_p2

!-------------------------------------------------------------------------------
!----- system-consistent vertexbased+weno+thinc limiter for P1P2:
!-------------------------------------------------------------------------------

subroutine thincvertexbasedweno_p2(ucons, uprim)

integer :: ie, ieqn, mmax, imat, matint(g_mmi%nummat)
real*8  :: dx2, almax, theta1(g_neqns), thetap1(g_nprim)
real*8  :: uneigh(2,g_neqns,-1:1), upneigh(2,g_nprim,-1:1), &
           ueq(2,-1:1), &
           uxxlim(g_neqns,0:imax+1), &
           pxxlim(g_nprim,0:imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)
logical :: trcell, is_intcell

associate (nummat=>g_mmi%nummat)

  uxxlim = 0.0
  pxxlim = 0.0

  !--- 1. P2 derivative limiting
  do ie = 1,imax

    !--- i. detect troubled-cell
    trcell = .false.
    if ((ie.eq.1) .or. (ie.eq.imax)) then
      trcell = .true.
    else
      do ieqn = g_mmi%irmin, g_mmi%irmax
        ueq = ucons(1:2,ieqn,ie-1:ie+1)
        trcell = trcell &
          .or. troubled_cell(2.0, ueq, coord(ie+1)-coord(ie))
      end do !ieqn
      do ieqn = 1,g_nprim
        ueq = uprim(1:2,ieqn,ie-1:ie+1)
        trcell = trcell &
          .or. troubled_cell(2.0, ueq, coord(ie+1)-coord(ie))
      end do !ieqn
    end if

    !--- ii. determine limiting appropriately
    ! conserved quantities
    if (trcell) then
      dx2 = 0.5 * (coord(ie)-coord(ie-1))
      uneigh(1:2,:,-1) = ucons(2:3,:,ie-1) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+1)-coord(ie))
      uneigh(1:2,:,0)  = ucons(2:3,:,ie) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
      uneigh(1:2,:,1)  = ucons(2:3,:,ie+1) / (dx2*dx2)

      call weno_fn(g_neqns, 200.0, uneigh)
      uxxlim(:,ie) = uneigh(2,:,0) * dx2 * dx2

    else
      uxxlim(:,ie) = ucons(3,:,ie)

    end if

    ! primitive quantities
    if (trcell) then
      dx2 = 0.5 * (coord(ie)-coord(ie-1))
      upneigh(1:2,:,-1) = uprim(2:3,:,ie-1) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+1)-coord(ie))
      upneigh(1:2,:,0)  = uprim(2:3,:,ie) / (dx2*dx2)

      dx2 = 0.5 * (coord(ie+2)-coord(ie+1))
      upneigh(1:2,:,1)  = uprim(2:3,:,ie+1) / (dx2*dx2)

      call weno_fn(g_nprim, 200.0, upneigh)
      pxxlim(:,ie) = upneigh(2,:,0) * dx2 * dx2

    else
      pxxlim(:,ie) = uprim(3,:,ie)

    end if

  end do !ie

  ucons(3,:,:) = uxxlim(:,:)
  uprim(3,:,:) = pxxlim(:,:)

  !--- 2. P1 derivative limiting
  do ie = 1,imax

    !--- i. obtain limiter function for individial unknowns
    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta1)

    ! primitive quantities
    upneigh(1:2,:,-1) = uprim(1:2,:,ie-1)
    upneigh(1:2,:,0)  = uprim(1:2,:,ie)
    upneigh(1:2,:,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, upneigh, thetap1)

    !--- ii. detect interface/single-material cell
    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    is_intcell = intrecons_cell(ucons(1,1:nummat,ie), matint)

    !--- iii. determine limiting appropriately
    if (.not.is_intcell) then

      ! use common limiter function for all volume-fractions
      theta1(1:nummat) = theta1(mmax) !minval(theta1(1:nummat))

    else

      ! no limiter for volume fractions in interface cells
      do imat = 1,nummat
        if (matint(imat) == 1) then
          theta1(imat) = 1.0
        end if
      end do !imat

    end if

    !--- iv. apply limiting
    do ieqn = 1,g_neqns
      ucons(2,ieqn,ie) = theta1(ieqn) * ucons(2,ieqn,ie)
    end do !ieqn

    do ieqn = 1,g_nprim
      uprim(2,ieqn,ie) = thetap1(ieqn) * uprim(2,ieqn,ie)
    end do !ieqn

  end do !ie

end associate

end subroutine thincvertexbasedweno_p2

!-------------------------------------------------------------------------------
!----- system-consistent overbee+vertexbased limiter:
!-------------------------------------------------------------------------------

subroutine oververtexbased(ucons, uprim)

integer :: ie, ieqn, mmax
real*8  :: theta(g_neqns), thetap(g_nprim), theta_al, thrho(2)
real*8  :: almax, dalmax, rho1, rho2, vel, rhoe1, rhoe2
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    !--- 1. detect interface/single-material cell
    mmax = maxloc(ucons(1,1:nummat,ie), 1)
    almax = ucons(1,mmax,ie)
    dalmax = maxval(ucons(2,1:nummat,ie)/(0.5 * (coord(ie+1)-coord(ie))))

    !--- 2. obtain limiter function for individual unknowns
    ! conserved quantities
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, uneigh, theta)

    ! compressive limiting for volume fraction
    call overbee_fn(uneigh(:,mmax,:), theta_al)

    ! primitive quantities
    uneigh(1:2,1:g_nprim,-1) = uprim(1:2,:,ie-1)
    uneigh(1:2,1:g_nprim,0)  = uprim(1:2,:,ie)
    uneigh(1:2,1:g_nprim,1)  = uprim(1:2,:,ie+1)

    call vertexbased_fn(g_nprim, uneigh(:,1:g_nprim,:), thetap)

    ! use common limiter function for all volume-fractions
    theta(1:nummat) = theta_al !theta(mmax)

    ! 3. Obtain consistent limiter functions for the equation system
    !    Interface detection
    if ( (g_nmatint .eq. 1) .and. &
         interface_cell(almax, dalmax) ) then

      call intfac_limiting(ucons(:,:,ie), theta, theta_al)

    else

      do ieqn = 1,g_neqns
        ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
      end do !ieqn

    end if

    uprim(2,:,ie) = thetap(:) * uprim(2,:,ie)

  end do !ie

end associate

end subroutine oververtexbased

!-------------------------------------------------------------------------------
!----- Vertex-based limiter for n equations individually
!----- this sub actually calculates the limiter function according to vertexbased.
!----- the input to this function can be modified to limit P1 or P2 dofs, refer
!----- to vertexbased_p1 and vertexbased_p2 respectively
!-------------------------------------------------------------------------------

subroutine vertexbased_fn(neq, ucons, theta)

integer, intent(in) :: neq
real*8,  intent(in) :: ucons(2,neq,-1:1)

integer :: ieqn, ifc, ne
real*8  :: ui, ug, umin, umax, diff, phi, theta(neq)

  theta(:) = 1.0

  do ieqn = 1,neq

    ui = ucons(1,ieqn,0)

    do ifc = 1,2

      ! find min and max in neighborhood
      if (ifc == 1) then
        ne = -1
      else if (ifc == 2) then
        ne = 1
      end if
      umax = max(ui, ucons(1,ieqn,ne))
      umin = min(ui, ucons(1,ieqn,ne))

      ! unlimited 2nd order solution
      ug = ucons(1,ieqn,0) + ((-1.0)**ifc) * ucons(2,ieqn,0)

      ! bounds
      diff = ug-ui
      if (diff > 1.0d-16) then
        phi = min(1.0, (umax-ui)/diff)

      else if (diff < -1.0d-16) then
        phi = min(1.0, (umin-ui)/diff)

      else
        phi = 1.0

      end if

      theta(ieqn) = min(phi, theta(ieqn))

    end do !ifc

  end do !ieqn

end subroutine vertexbased_fn

!-------------------------------------------------------------------------------
!----- WENO limiter function
!-------------------------------------------------------------------------------

subroutine weno_fn(neq, wenocp1, ucons)

integer, intent(in) :: neq
real*8, intent(in) :: wenocp1

integer :: ieqn, is, nsten, imat
real*8  :: wi,epsweno,wt, &
           weight(3), &
           gradv(3),osc(3),gradu, &
           ulim(2,neq), &
           ucons(2,neq,-1:1)

associate (nummat=>g_mmi%nummat)

  epsweno = 1.d-10
  nsten = 3

  do ieqn = 1,neq

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

    ulim(2,ieqn) = gradu

  end do !ieqn

  ucons(2,:,0) = ulim(2,:)

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
!----- LINC reconstruction for volume fraction
!-------------------------------------------------------------------------------

subroutine linc_reconstruction(udof, pdof, basis, uho, pho, dx, xc)

real*8, intent(in) :: udof(g_tdof,g_neqns), pdof(g_tdof,g_nprim), basis(g_tdof), &
  dx, xc

real*8, intent(out) :: uho(g_neqns), pho(g_nprim)

integer :: imat, mmax, matint(g_mmi%nummat)
real*8 :: beta_linc, lolim, hilim, almax, alm, al_reco, xt, nx, x, sig, vel, &
  alsum, velavg, almat(g_mmi%nummat)

associate (nummat=>g_mmi%nummat)

  !--- TVD reconstruction
  call ho_reconstruction(g_neqns, udof, basis, uho)
  if (g_pureco == 1) then
    call ho_reconstruction(g_nprim, pdof, basis, pho)
  else
    call get_uprim_mm6eq(uho, pho)
  end if

  beta_linc = 2.5

  lolim = 2.0*1e-8
  hilim = 1.0 - lolim

  almat = udof(1,g_mmi%iamin:g_mmi%iamax)
  mmax = maxloc(almat, 1)
  almax = maxval(almat)

  !--- algebraic interface reconstruction, if required and appropriate
  if ((g_intreco==1) .and. intrecons_cell(udof(1,g_mmi%iamin:g_mmi%iamax), matint)) then

    x = (0.5*dx*basis(2)) + xc

    alsum = 0.0

    do imat = 1,nummat
      if (matint(imat) == 1) then
      ! THINC reconstruction for materials with interface in cell
      ! TVD reconstruction for materials without interface in cell
        alm = udof(1,imat)
        sig = dsign(1.0,udof(2,imat))

        nx = udof(2,imat)/(dabs(udof(2,imat))+1e-12)

        !! LINC reconstruction
        !xt = (1.0 + 0.5*nx*beta_linc - 2.0*alm) / (nx*beta_linc)
        !al_reco = 0.5 * (1.0 + nx * beta_linc * ((x-(xc-0.5*dx))/dx - xt))

        !if (al_reco > hilim) then
        !  al_reco = hilim
        !else if (al_reco < lolim) then
        !  al_reco = lolim
        !end if

        ! THINC reconstruction
        !xt = dlog( dexp(beta_linc*(1.0+sig-2.0*alm)/sig) &
        !  / (1.0-dexp(beta_linc*(1.0-sig-2.0*alm)/sig)) ) &
        !  / (2.0*beta_linc)
        xt = dlog( (dexp(2.0*sig*beta_linc*alm) - 1.0) &
          / (dexp(2.0*sig*beta_linc) - dexp(2.0*sig*beta_linc*alm)) ) &
          / (2.0*beta_linc)
        al_reco = 0.5 * (1.0 + dtanh( beta_linc*(sig*(x-(xc-0.5*dx))/dx+xt) ))

        ! volfrac
        uho(imat) = min(1.0-1e-14, max(1e-14, al_reco))

      end if

      alsum = alsum + uho(imat)
    end do !imat

    ! enforce unit sum
    uho(mmax) = uho(mmax) + (1.0 - alsum)
    alsum = 1.0

    velavg = udof(1,g_mmi%imome)/sum(udof(1,g_mmi%irmin:g_mmi%irmax))

    ! consistent reconstruction
    do imat = 1, nummat
      if (matint(imat) == 1) then
        alm = udof(1,imat)
        ! density
        uho(g_mmi%irmin+imat-1) = udof(1,g_mmi%irmin+imat-1)/alm * uho(imat)
        ! energy
        uho(g_mmi%iemin+imat-1) = udof(1,g_mmi%iemin+imat-1)/alm * uho(imat)
        ! pressure
        if (g_pureco == 1) then
          pho(apr_idx(nummat, imat)) = pdof(1,apr_idx(nummat, imat))/alm &
            * uho(imat)
        else
          pho(apr_idx(nummat, imat)) = eos3_alphapr(g_gam(imat), g_pc(imat), &
            alm, udof(1,g_mmi%irmin+imat-1), udof(1,g_mmi%iemin+imat-1), &
            velavg)/alm &
            * uho(imat)
        end if
      end if
    end do !imat

    vel = udof(1,g_mmi%imome)/sum( udof(1,g_mmi%irmin:g_mmi%irmax) )
    ! bulk-momentum
    uho(g_mmi%imome) = vel * sum( uho(g_mmi%irmin:g_mmi%irmax) )
    ! velocity
    if (g_pureco == 1) then
      pho(vel_idx(nummat, 0)) = pdof(1,vel_idx(nummat, 0))
    else
      pho(vel_idx(nummat, 0)) = velavg
    end if

  end if

end associate

end subroutine linc_reconstruction

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

subroutine intfac_limiting_p2(ucons, is_weno, theta1_al, theta2_al)

real*8,  intent(in) :: theta1_al, theta2_al
logical, intent(in) :: is_weno

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
    if (is_weno) then
      ucons(2,imat,1) = theta1_al * ucons(2,imat,1)
    else
      ucons(2,imat,1) = max(theta1_al, theta2_al) * ucons(2,imat,1)
    end if
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
!----- Determine reconstruction stencil based on interface and material
!-------------------------------------------------------------------------------

subroutine recons_stencil(ie, matint_cell, lefte, righte, skew)

integer, intent(in) :: ie
logical, intent(in) :: matint_cell(-2:2)

integer, intent(out) :: lefte, righte, skew

  ! 0. default central stencil
  lefte = ie-1
  righte = ie+1
  skew = 0
  ! only for interior cells, try to modify the limiting stencil according to
  ! material interfaces
  if ((ie.gt.1) .and. (ie.lt.imax)) then
    ! 1. for material interface cells
    if (matint_cell(0)) then
      ! if this is the only material interface cell in the vicinity, diffuse
      ! the interface
      if (.not.matint_cell(-1) .and. .not.matint_cell(1)) then
        lefte = ie-1
        righte = ie+1
        skew = 0
      ! if the left cell is not an interface cell, extend stencil to right,
      ! but only if the right cell is an interface cell
      else if (.not.matint_cell(-1) .and. matint_cell(2)) then
        lefte = ie+2
        skew = 1
      ! if the right cell is not an interface cell, extend stencil to left,
      ! but only if the left cell is an interface cell
      else if (.not.matint_cell(1) .and. matint_cell(-2)) then
        righte = ie-2
        skew = -1
      end if
    ! 2. for 'pure' cells
    else
      if (matint_cell(-1) .and. matint_cell(1)) then
        lefte = ie-1
        righte = ie+1
        skew = 0
      else if (matint_cell(-1) .and. .not.matint_cell(2)) then
        lefte = ie+2
        skew = 1
      else if (matint_cell(1) .and. .not.matint_cell(-2)) then
        righte = ie-2
        skew = -1
      end if
    end if
  end if

end subroutine recons_stencil

!-------------------------------------------------------------------------------
!----- Interface-cell / mixed-cell indicator based on vol-frac gradient
!-------------------------------------------------------------------------------

logical function interface_cell(alcell, dalcell)

real*8, intent(in) :: alcell, dalcell
real*8  :: scale_al

logical :: al, dal

  scale_al = 0.0001
  dal = dabs(dalcell) .gt. scale_al

  scale_al = 0.0001 !10000.0*g_alphamin
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

real*8, intent(in) :: beta, u(g_tdof,-1:1), dx

integer :: ifc
real*8  :: dplus, dminu, difc, a1, a2, a3, mm

  troubled_cell = .false.

  dplus = u(1,0)-u(1,-1)
  dminu = u(1,1)-u(1,0)

  do ifc = 1,2
    difc = ((-1.0)**ifc) * u(2,0)

    a1 = difc
    a2 = dplus
    a3 = dminu

    ! modified minmod
    if (dabs(a1) <= beta*dx*dx) then
      mm = a1

    else
      if (a1*a2 > 0.0 .and. a2*a3 > 0.0) then
        mm = dsign(1.0, a1) * min(dabs(a1), dabs(a2), dabs(a3))
      else
        mm = 0.0
      end if

    end if

    if (dabs(mm-difc) > 1.0e-12) then
      troubled_cell = troubled_cell .or. .true.
    else
      troubled_cell = troubled_cell .or. .false.
    end if
  end do !ifc

end function troubled_cell

!-------------------------------------------------------------------------------
!----- Detect cells for algebraic interface reconstruction
!-------------------------------------------------------------------------------

logical function intrecons_cell(almat, matint)

real*8, intent(in) :: almat(g_mmi%nummat)
integer, intent(out) :: matint(g_mmi%nummat)

integer :: imat
real*8 :: lolim, hilim, almax

  lolim = 2.0*1e-8
  hilim = 1.0 - lolim

  almax = maxval(almat)
  matint = 0

  do imat = 1, g_mmi%nummat
    if ((almat(imat) > lolim) .and. (almat(imat) < hilim)) matint(imat) = 1
  end do !imat

  if ((almax > lolim) .and. (almax < hilim)) then
    intrecons_cell = .true.
  else
    intrecons_cell = .false.
  end if

end function intrecons_cell

!-------------------------------------------------------------------------------
!----- Fill array of cells indicating material interface cells
!-------------------------------------------------------------------------------

subroutine fill_matintel(ucons, matint_el)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: matint_el(0:imax+1), matint(g_mmi%nummat), ie

  matint_el = 0

  do ie = 0, imax+1
    ! use indicator for algebraic reconstruction to determine interface cells
    if (intrecons_cell(ucons(1,g_mmi%iamin:g_mmi%iamax,ie), matint)) &
      matint_el(ie) = 1
  end do !ie

end subroutine fill_matintel

!-------------------------------------------------------------------------------

END MODULE reconstruction_mm6eq
