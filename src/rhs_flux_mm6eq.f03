!!------------------------------------------------------------------------------
!!----- Pressure non-equilibrium multi-material (Flux computation module)
!!----- by
!!----- Aditya K Pandare
!!------------------------------------------------------------------------------

MODULE rhs_flux_mm6eq

USE reconstruction_mm6eq

implicit none

CONTAINS

!-------------------------------------------------------------------------------
!----- P0 RHS:
!-------------------------------------------------------------------------------

subroutine rhs_p0_mm6eq(ucons, ulim, rhsel)

real*8  :: rhsel(g_gdof,g_neqns,imax), &
           ulim(g_tdof,g_neqns,0:imax+1)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  call flux_p0_mm6eq(ucons, ulim, rhsel)

  if (g_nprelx .eq. 1) then
    call relaxpressure_p0(ulim, rhsel)
  end if

end subroutine rhs_p0_mm6eq

!-------------------------------------------------------------------------------
!----- P0P1 RHS:
!-------------------------------------------------------------------------------

subroutine rhs_p0p1_mm6eq(ucons, ulim, rhsel)

real*8  :: rhsel(g_gdof,g_neqns,imax), &
           ulim(g_tdof,g_neqns,0:imax+1)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  call boundpreserve_alpha_p1(ucons)

  call flux_p0p1_mm6eq(ucons, ulim, rhsel)

  if (g_nprelx .eq. 1) then
    call relaxpressure_p0(ulim, rhsel)
  end if

end subroutine rhs_p0p1_mm6eq

!-------------------------------------------------------------------------------
!----- P1 RHS:
!-------------------------------------------------------------------------------

subroutine rhs_p1_mm6eq(ucons, ulim, rhsel)

real*8  :: rhsel(g_gdof,g_neqns,imax), &
           ulim(g_tdof,g_neqns,0:imax+1)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  call boundpreserve_alpha_p1(ucons)

  call flux_p1_mm6eq(ucons, ulim, rhsel)

  if (g_nprelx .eq. 1) then
    call relaxpressure_p1(ulim, rhsel)
  end if

end subroutine rhs_p1_mm6eq

!-------------------------------------------------------------------------------
!----- P1P2 RHS:
!-------------------------------------------------------------------------------

subroutine rhs_p1p2_mm6eq(ucons, ulim, rhsel)

real*8  :: rhsel(g_gdof,g_neqns,imax), &
           ulim(g_tdof,g_neqns,0:imax+1)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  call flux_p1p2_mm6eq(ucons, ulim, rhsel)

  if (g_nprelx .eq. 1) then
    call relaxpressure_p1p2(ulim, rhsel)
  end if

end subroutine rhs_p1p2_mm6eq

!-------------------------------------------------------------------------------
!----- P0 Advective-flux contribution to RHS:
!-------------------------------------------------------------------------------

subroutine flux_p0_mm6eq(ucons, ulim, rhsel)

integer :: ifc, iel, ier, ieqn
real*8  :: ul(g_neqns), ur(g_neqns), &
           ncnflux(g_neqns,2), intflux(g_neqns), lplus, lminu, lmag, &
           rhsel(g_gdof,g_neqns,imax), &
           ulim(g_tdof,g_neqns,0:imax+1)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

  ulim = ucons

  do ifc = 1,imax+1

  intflux = 0.0
  ncnflux = 0.0

  iel = ifc - 1
  ier = ifc

  ul(:) = ucons(1,:,iel)
  ur(:) = ucons(1,:,ier)

  !--- fluxes

  if (i_flux .eq. 1) then
     call llf_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
     call llf_nonconserv(ul, ur, ul, ur, lplus, lminu, lmag, ncnflux)
  else if (i_flux .eq. 2) then
     call ausmplus_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
     call ausmplus_nonconserv(ul, ur, ul, ur, lplus, lminu, lmag, ncnflux)
  else
     write(*,*) "Invalid flux scheme."
     stop
  endif

  if (iel .gt. 0) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - ncnflux(ieqn,1)
    end do !ieqn
  end if

  if (ier .lt. (imax+1)) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + ncnflux(ieqn,2)
    end do !ieqn
  end if

  end do !ifc

end subroutine flux_p0_mm6eq

!-------------------------------------------------------------------------------
!----- P0P1 Advective-flux contribution to RHS:
!-------------------------------------------------------------------------------

subroutine flux_p0p1_mm6eq(ucons, ulim, rhsel)

integer :: ifc, iel, ier, ieqn
real*8  :: ul(g_neqns), ur(g_neqns), uavgl(g_neqns), uavgr(g_neqns), &
           ncnflux(g_neqns,2), intflux(g_neqns), rhsel(g_gdof,g_neqns,imax), &
           lplus, lminu, lmag

real*8  :: ulim(g_tdof,g_neqns,0:imax+1), ucons(g_tdof,g_neqns,0:imax+1)

  call reconstruction_p0p1(ucons)
  ulim = ucons

  do ifc = 1,imax+1

  intflux = 0.0
  ncnflux = 0.0

  iel = ifc - 1
  ier = ifc

  do ieqn = 1,g_neqns
    ul(ieqn) = ucons(1,ieqn,iel) + ucons(2,ieqn,iel)
    ur(ieqn) = ucons(1,ieqn,ier) - ucons(2,ieqn,ier)

    uavgl(ieqn) = ucons(1,ieqn,iel)
    uavgr(ieqn) = ucons(1,ieqn,ier)
  end do !ieqn

  !--- fluxes

  if (i_flux .eq. 1) then
     call llf_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
     call llf_nonconserv(ul, ur, uavgl, uavgr, &
                         lplus, lminu, lmag, ncnflux)
  else if (i_flux .eq. 2) then
     call ausmplus_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
     call ausmplus_nonconserv(ul, ur, uavgl, uavgr, &
                              lplus, lminu, lmag, ncnflux)
  else
     write(*,*) "Invalid flux scheme."
     stop
  endif

  if (iel .gt. 0) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - ncnflux(ieqn,1)
    end do !ieqn
  end if

  if (ier .lt. (imax+1)) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + ncnflux(ieqn,2)
    end do !ieqn
  end if

  end do !ifc

end subroutine flux_p0p1_mm6eq

!-------------------------------------------------------------------------------
!----- P1 Advective-flux contribution to RHS:
!-------------------------------------------------------------------------------

subroutine flux_p1_mm6eq(ucons, ulim, rhsel)

integer :: ie
real*8, intent(in)  :: ucons(g_tdof,g_neqns,0:imax+1)

real*8  :: riemanngrad(g_mmi%nummat+1,imax), &
           vriemann(g_mmi%nummat+1,imax+1), &
           ulim(g_tdof,g_neqns,0:imax+1), &
           rhsel(g_gdof,g_neqns,imax)

  !--- limiting
  ulim = ucons
  call limiting_p1(ulim)

  riemanngrad = 0.0
  vriemann = 0.0

  !--- surface integration
  call surfaceint_p1(ulim, riemanngrad, vriemann, rhsel)
  !--- volume integration
  call volumeint_p1(ulim, riemanngrad, vriemann, rhsel)

end subroutine flux_p1_mm6eq

!-------------------------------------------------------------------------------
!----- P1P2 Advective-flux contribution to RHS:
!-------------------------------------------------------------------------------

subroutine flux_p1p2_mm6eq(ucons, ulim, rhsel)

real*8, intent(in)  :: ucons(g_tdof,g_neqns,0:imax+1)

real*8  :: riemanngrad(g_mmi%nummat+1,imax), &
           ulim(g_tdof,g_neqns,0:imax+1), &
           rhsel(g_gdof,g_neqns,imax)

  !--- limiting
  ulim = ucons
  call reconstruction_p1p2(ulim)

  riemanngrad = 0.0

  !--- surface integration
  call surfaceint_p1p2(ulim, riemanngrad, rhsel)
  !--- volume integration
  call volumeint_p1p2(ulim, riemanngrad, rhsel)

end subroutine flux_p1p2_mm6eq

!-------------------------------------------------------------------------------
!----- P1 surface contribution to RHS:
!-------------------------------------------------------------------------------

subroutine surfaceint_p1(ucons, rgrad, vriem, rhsel)

integer :: ifc, iel, ier, ieqn, imat
real*8  :: ul(g_neqns), ur(g_neqns), uavgl(g_neqns), uavgr(g_neqns), &
           intflux(g_neqns), rhsel(g_gdof,g_neqns,imax), &
           lplus, lminu, lmag, &
           alpha_star, u_star

real*8  :: rgrad(g_mmi%nummat+1,imax), vriem(g_mmi%nummat+1,imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ifc = 1,imax+1

  intflux = 0.0

  iel = ifc - 1
  ier = ifc

  do ieqn = 1,g_neqns
    ul(ieqn) = ucons(1,ieqn,iel) + ucons(2,ieqn,iel)
    ur(ieqn) = ucons(1,ieqn,ier) - ucons(2,ieqn,ier)

    uavgl(ieqn) = ucons(1,ieqn,iel)
    uavgr(ieqn) = ucons(1,ieqn,ier)
  end do !ieqn

  !--- fluxes

  if (i_flux .eq. 2) then
     call ausmplus_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
  else
     write(*,*) "Invalid flux scheme."
     stop
  endif

  do imat = 1,nummat
    alpha_star = dabs(lplus) * ul(imat) + dabs(lminu) * ur(imat)
    vriem(imat,ifc) = alpha_star
  end do !imat
  u_star = lmag*(lplus+lminu)
  vriem(nummat+1,ifc) = u_star

  if (iel .gt. 0) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
          rhsel(2,ieqn,iel) = rhsel(2,ieqn,iel) - intflux(ieqn)
    end do !ieqn

    !--- compute gradients of volume fractions and velocity for the
    !--- non-conservative terms from Riemann reconstructed values
    do imat = 1,nummat
      alpha_star = dabs(lplus) * ul(imat) + dabs(lminu) * ur(imat)
      rgrad(imat,iel) = rgrad(imat,iel) + alpha_star
    end do !imat
    rgrad(nummat+1,iel) = rgrad(nummat+1,iel) + u_star

  end if

  if (ier .lt. (imax+1)) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
          rhsel(2,ieqn,ier) = rhsel(2,ieqn,ier) - intflux(ieqn)
    end do !ieqn

    !--- compute gradients of volume fractions and velocity for the
    !--- non-conservative terms from Riemann reconstructed values
    do imat = 1,nummat
      alpha_star = dabs(lplus) * ul(imat) + dabs(lminu) * ur(imat)
      rgrad(imat,ier) = rgrad(imat,ier) - alpha_star
    end do !imat
    rgrad(nummat+1,ier) = rgrad(nummat+1,ier) - u_star

  end if

  end do !ifc

end associate

end subroutine surfaceint_p1

!-------------------------------------------------------------------------------
!----- P1 volume contribution to RHS:
!-------------------------------------------------------------------------------

subroutine volumeint_p1(ucons, rgrad, vriem, rhsel)

integer :: ig, ie, ieqn, ngauss, imat
data       ngauss/2/

real*8  :: dx2, p, hmat, viriem, &
           u(g_neqns), up(g_neqns), &
           carea(2), weight(2), &
           cflux(g_neqns), &
           nflux(g_gdof,g_neqns), &
           rhsel(g_gdof,g_neqns,imax)

real*8, intent(in) :: rgrad(g_mmi%nummat+1,imax), vriem(g_mmi%nummat+1,imax+1), &
                      ucons(g_tdof,g_neqns,0:imax+1)

associate (nummat=>g_mmi%nummat)

  ngauss = 2

  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax
  cflux = 0.0
  nflux = 0.0
  do ig = 1,ngauss

    dx2 = weight(ig) ! <-- 2.0/dx * weight(ig)/2.0 * dx

    ! basis function

    do ieqn = 1,g_neqns
      u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie)
    end do !ieqn
    viriem = 0.5* (vriem(nummat+1,ie) + vriem(nummat+1,ie+1)) + carea(ig) * rgrad(nummat+1,ie)/2.0

    call get_uprim_mm6eq(u, up)

    ! bulk pressure
    p = 0.0
    do imat = 1,nummat
    p = p + u(imat)*up(g_mmi%irmin+imat-1)
    end do !imat

    !--- flux terms
    ! momentum flux
    cflux(g_mmi%imome) = up(g_mmi%imome) * u(g_mmi%imome) + p
    nflux(1,g_mmi%imome) = 0.0
    nflux(2,g_mmi%imome) = 0.0
    do imat = 1,nummat
      hmat = u(g_mmi%iemin+imat-1) + u(imat)*up(g_mmi%irmin+imat-1)
      ! other conservative fluxes
      cflux(imat) = 0.0
      cflux(g_mmi%irmin+imat-1) = up(g_mmi%imome) * u(g_mmi%irmin+imat-1)
      cflux(g_mmi%iemin+imat-1) = up(g_mmi%imome) * hmat!(u(g_mmi%iemin+imat-1) )
      ! non-conservative fluxes
      nflux(1,imat) = u(imat) * rgrad(nummat+1,ie)
      nflux(1,g_mmi%irmin+imat-1) = 0.0
      nflux(1,g_mmi%iemin+imat-1) = p * up(g_mmi%imome) * rgrad(imat,ie)
      nflux(2,imat) = nflux(1,imat) * carea(ig) + u(imat) * viriem * 2.0 !up(g_mmi%imome) * 2.0
      nflux(2,g_mmi%irmin+imat-1) = 0.0
      nflux(2,g_mmi%iemin+imat-1) = nflux(1,g_mmi%iemin+imat-1) * carea(ig)! &
                                    !+ p * up(g_mmi%imome) * u(imat) * 2.0
    end do !imat

    do ieqn = 1,g_neqns
      rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + dx2 * cflux(ieqn)
    end do !ieqn

    do ieqn = 1,g_neqns
      rhsel(1,ieqn,ie) = rhsel(1,ieqn,ie) + 0.5 * dx2 * nflux(1,ieqn)
      rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + 0.5 * dx2 * nflux(2,ieqn)
    end do !ieqn

  end do !ig
  end do !ie

end associate

end subroutine volumeint_p1

!-------------------------------------------------------------------------------
!----- P1P2 surface contribution to RHS:
!-------------------------------------------------------------------------------

subroutine surfaceint_p1p2(ucons, rgrad, rhsel)

integer :: ifc, iel, ier, ieqn
real*8  :: ul(g_neqns), ur(g_neqns), uavgl(g_neqns), uavgr(g_neqns), &
           intflux(g_neqns), rhsel(g_gdof,g_neqns,imax), &
           lplus, lminu, lmag, &
           alpha_star, u_star

real*8  :: rgrad(2,imax), ucons(g_tdof,g_neqns,0:imax+1)

  do ifc = 1,imax+1

  intflux = 0.0

  iel = ifc - 1
  ier = ifc

  do ieqn = 1,g_neqns
    ul(ieqn) = ucons(1,ieqn,iel) + ucons(2,ieqn,iel) + 1.0/3.0*ucons(3,ieqn,iel)
    ur(ieqn) = ucons(1,ieqn,ier) - ucons(2,ieqn,ier) + 1.0/3.0*ucons(3,ieqn,ier)

    uavgl(ieqn) = ucons(1,ieqn,iel)
    uavgr(ieqn) = ucons(1,ieqn,ier)
  end do !ieqn

  !--- fluxes

  if (i_flux .eq. 2) then
     call ausmplus_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
  else
     write(*,*) "Invalid flux scheme."
     stop
  endif

  !--- compute gradients of volume fraction and velocity for the
  !--- non-conservative terms from Riemann reconstructed values
  alpha_star = dabs(lplus) * ul(1) + dabs(lminu) * ur(1)
  u_star = lmag*(lplus+lminu)

  if (iel .gt. 0) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
          rhsel(2,ieqn,iel) = rhsel(2,ieqn,iel) - intflux(ieqn)
    end do !ieqn

    rgrad(1,iel) = rgrad(1,iel) + alpha_star
    rgrad(2,iel) = rgrad(2,iel) + u_star

  end if

  if (ier .lt. (imax+1)) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
          rhsel(2,ieqn,ier) = rhsel(2,ieqn,ier) - intflux(ieqn)
    end do !ieqn

    rgrad(1,ier) = rgrad(1,ier) - alpha_star
    rgrad(2,ier) = rgrad(2,ier) - u_star

  end if

  end do !ifc

end subroutine surfaceint_p1p2

!-------------------------------------------------------------------------------
!----- P1P2 volume contribution to RHS:
!-------------------------------------------------------------------------------

subroutine volumeint_p1p2(ucons, rgrad, rhsel)

integer :: ig, ie, ieqn, ngauss
data       ngauss/2/

real*8  :: dx2, b3, p, &
           u(g_neqns), up(g_neqns), &
           carea(2), weight(2), &
           cflux(g_neqns), &
           nflux(g_gdof,g_neqns), &
           rhsel(g_gdof,g_neqns,imax)

real*8, intent(in) :: rgrad(2,imax), ucons(g_tdof,g_neqns,0:imax+1)

  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax
  do ig = 1,ngauss

    dx2 = weight(ig) ! <-- 2.0/dx * weight(ig)/2.0 * dx

    ! basis function

    b3 = 0.5*carea(ig)*carea(ig) - 1.0/6.0
    do ieqn = 1,g_neqns
      u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie) + b3 * ucons(3,ieqn,ie)
    end do !ieqn

    call get_uprim_mm6eq(u, up)

    p = up(1)*up(2) + (1.0-up(1))*up(3)

    ! conservative fluxes
    cflux(1) = 0.0
    cflux(2) = up(4) * u(2)
    cflux(3) = up(4) * u(3)
    cflux(4) = up(4) * u(4) + p
    cflux(5) = up(4) * (u(5) )
    cflux(6) = up(4) * (u(6) )

    ! non-conservative fluxes
    nflux(1,1) = up(1) * rgrad(2,ie)
    nflux(1,2) = 0.0
    nflux(1,3) = 0.0
    nflux(1,4) = 0.0
    nflux(1,5) = p * up(4) * rgrad(1,ie)
    nflux(1,6) = p * up(4) * (-rgrad(1,ie))

    nflux(2,1) = nflux(1,1) * carea(ig) + up(1) * up(4) * 2.0
    nflux(2,2) = 0.0
    nflux(2,3) = 0.0
    nflux(2,4) = 0.0
    nflux(2,5) = nflux(1,5) * carea(ig) + p * up(4) * up(1) * 2.0
    nflux(2,6) = nflux(1,6) * carea(ig) + p * up(4) * (1.0-up(1)) * 2.0

    do ieqn = 1,g_neqns
      rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + dx2 * cflux(ieqn)
    end do !ieqn

    do ieqn = 1,g_neqns
      rhsel(1,ieqn,ie) = rhsel(1,ieqn,ie) + 0.5 * dx2 * nflux(1,ieqn)
      rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + 0.5 * dx2 * nflux(2,ieqn)
    end do !ieqn

  end do !ig
  end do !ie

end subroutine volumeint_p1p2

!-------------------------------------------------------------------------------
!----- 2fluid Lax-Friedrichs flux:
!-------------------------------------------------------------------------------

subroutine llf_mm6eq(ul, ur, flux, lplus, lminu, lmag)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns)

integer :: imat
real*8 :: flux(g_neqns), lplus, lminu, lmag
real*8 :: up_l(g_neqns), up_r(g_neqns)
real*8 :: ffunc_l(g_neqns), ffunc_r(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,rhom_l,em_l,am_l,hm_l,pm_l, &
                                   arhom_r,rhom_r,em_r,am_r,hm_r,pm_r, &
                                   am_12,rhom_12,al_12

real*8 :: rhou_l, u_l, rho_l, p_l, &
          rhou_r, u_r, rho_r, p_r
real*8 :: rho_12,ac_12
  
real*8 :: lambda

associate (nummat=>g_mmi%nummat)

  flux(:) = 0.0

  do imat = 1,nummat
  ! ul
    al_l(imat)    = ul(imat)
    arhom_l(imat) = ul(g_mmi%irmin+imat-1)
    em_l(imat)    = ul(g_mmi%iemin+imat-1)
  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r(imat)    = ur(g_mmi%iemin+imat-1)
  end do !imat
  rhou_l  = ul(g_mmi%imome)
  rhou_r  = ur(g_mmi%imome)

  call get_uprim_mm6eq(ul, up_l)
  call get_uprim_mm6eq(ur, up_r)

  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l    = rhou_l / rho_l
  u_r    = rhou_r / rho_r

  do imat = 1,nummat
    rhom_l(imat) = arhom_l(imat) / al_l(imat)
    pm_l(imat)   = up_l(g_mmi%irmin+imat-1)
    am_l(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_l(imat), pm_l(imat))
    hm_l(imat)   = em_l(imat) + al_l(imat)*pm_l(imat)

    rhom_r(imat) = arhom_r(imat) / al_r(imat)
    pm_r(imat)   = up_r(g_mmi%irmin+imat-1)
    am_r(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_r(imat), pm_r(imat))
    hm_r(imat)   = em_r(imat) + al_r(imat)*pm_r(imat)
  end do !imat
  p_l = sum(pm_l)
  p_r = sum(pm_r)

  ! average states
  rho_12  = 0.5*(rho_l + rho_r)
  do imat = 1,nummat
    rhom_12(imat) = 0.5*(rhom_l(imat) + rhom_r(imat))
    am_12(imat)   = 0.5*(am_l(imat) + am_r(imat))
    al_12(imat)   = 0.5*(al_l(imat) + al_r(imat))
  end do !imat

  ! numerical speed of sound choice:
  ac_12 = 0.0
  do imat = 1,nummat
    ac_12 = ac_12 + ( al_12(imat)*rhom_12(imat)*am_12(imat)*am_12(imat) )
  end do !imat
  ac_12 = dsqrt( ac_12 / rho_12 )

  lambda = ac_12 + max(dabs(u_l),dabs(u_r));

  ! flux functions
  do imat = 1,nummat
    ffunc_l(imat) = u_l * al_l(imat)
    ffunc_l(g_mmi%irmin+imat-1) = u_l * arhom_l(imat)
    ffunc_l(g_mmi%iemin+imat-1) = u_l * hm_l(imat)

    ffunc_r(imat) = u_r * al_r(imat)
    ffunc_r(g_mmi%irmin+imat-1) = u_r * arhom_r(imat)
    ffunc_r(g_mmi%iemin+imat-1) = u_r * hm_r(imat)
  end do !imat
  ffunc_l(g_mmi%imome) = u_l * rhou_l + p_l
  ffunc_r(g_mmi%imome) = u_r * rhou_r + p_r

  flux = 0.5 * ( ffunc_l+ffunc_r - lambda*(ur-ul) )

  lplus = 0.5
  lminu = 0.5

  lmag = lambda

end associate

end subroutine llf_mm6eq

!-------------------------------------------------------------------------------

subroutine llf_nonconserv(ul, ur, uavgl, uavgr, lplus, lminu, lmag, ncnflux)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), &
                      uavgl(g_neqns), uavgr(g_neqns), &
                      lplus, lminu, lmag

integer :: imat
real*8  :: ncnflux(g_neqns,2), &
           uprim_avg(g_neqns), &
           u_conv_l, u_conv_r, uf, p, &
           alpf(g_mmi%nummat), ncnfl, ncnfr

associate (nummat=>g_mmi%nummat)

  do imat = 1,nummat
    alpf(imat) = lplus*ul(imat) + lminu*ur(imat)
  end do !imat

  ! left element
  call get_uprim_mm6eq(uavgl, uprim_avg)
  u_conv_l = uprim_avg(g_mmi%imome)
  p = 0.0
  do imat = 1,nummat
    p = p + uavgl(imat)*uprim_avg(g_mmi%irmin+imat-1)
  end do !imat
  ncnfl = p * u_conv_l

  ! right element
  call get_uprim_mm6eq(uavgr, uprim_avg)
  u_conv_r = uprim_avg(g_mmi%imome)
  p = 0.0
  do imat = 1,nummat
    p = p + uavgr(imat)*uprim_avg(g_mmi%irmin+imat-1)
  end do !imat
  ncnfr = p * u_conv_r

  uf = 0.5 * (u_conv_r-u_conv_l)

  do imat = 1,nummat
    ncnflux(imat,1)               = - ur(imat) * uf
    ncnflux(g_mmi%iemin+imat-1,1) = - alpf(imat) * ncnfl

    ncnflux(imat,2)               = + ul(imat) * uf
    ncnflux(g_mmi%iemin+imat-1,2) = - alpf(imat) * ncnfr
  end do !imat

end associate

end subroutine llf_nonconserv

!-------------------------------------------------------------------------------
!----- 2fluid AUSM+UP:
!-------------------------------------------------------------------------------

subroutine ausmplus_mm6eq(ul, ur, flux, lambda_plus, lambda_minu, lambda_mag)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns)

integer :: imat
real*8 :: flux(g_neqns)
real*8 :: up_l(g_neqns), up_r(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,rhom_l,em_l,am_l,hm_l, &
                                   arhom_r,rhom_r,em_r,am_r,hm_r, &
                                   am_12,rhom_12,al_12

real*8 :: rhou_l, u_l, m_l, rho_l, pi_l, p_l, &
          rhou_r, u_r, m_r, rho_r, pi_r, p_r
real*8 :: rho_12, ac_12, &
          f_a, m_12, p_12, m_p, p_u
real*8 :: msplus_l(3),msplus_r(3),msminu_l(3),msminu_r(3)
real*8 :: psplus_l,psplus_r,psminu_l,psminu_r
real*8 :: temp

real*8 :: lambda,lambda_plus, lambda_minu, lambda_mag

real*8 :: k_p, k_u

associate (nummat=>g_mmi%nummat)

  k_p = 1.0;
  k_u = 0.1;

  flux(:) = 0.0

  call get_uprim_mm6eq(ul, up_l)
  call get_uprim_mm6eq(ur, up_r)

  ! material states and bulk pressure
  p_l = 0.0
  p_r = 0.0
  do imat = 1,nummat
  ! ul
    al_l(imat)    = ul(imat)
    arhom_l(imat) = ul(g_mmi%irmin+imat-1)
    em_l(imat)    = ul(g_mmi%iemin+imat-1)

    rhom_l(imat) = arhom_l(imat) / al_l(imat)
    pi_l         = up_l(g_mmi%irmin+imat-1)
    am_l(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_l(imat), pi_l)
    hm_l(imat)   = em_l(imat) + al_l(imat)*pi_l
    p_l = p_l + al_l(imat)*pi_l

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r(imat)    = ur(g_mmi%iemin+imat-1)

    rhom_r(imat) = arhom_r(imat) / al_r(imat)
    pi_r         = up_r(g_mmi%irmin+imat-1)
    am_r(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_r(imat), pi_r)
    hm_r(imat)   = em_r(imat) + al_r(imat)*pi_r
    p_r = p_r + al_r(imat)*pi_r
  end do !imat

  ! bulk state
  rhou_l  = ul(g_mmi%imome)
  rhou_r  = ur(g_mmi%imome)

  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l    = rhou_l / rho_l
  u_r    = rhou_r / rho_r

  ! average states
  rho_12  = 0.5*(rho_l + rho_r)
  do imat = 1,nummat
    rhom_12(imat) = 0.5*(rhom_l(imat) + rhom_r(imat))
    am_12(imat)   = 0.5*(am_l(imat) + am_r(imat))
    al_12(imat)   = 0.5*(al_l(imat) + al_r(imat))
  end do !imat

  ! numerical speed of sound choice:
  ac_12 = 0.0
  do imat = 1,nummat
    ac_12 = ac_12 + ( al_12(imat)*rhom_12(imat)*am_12(imat)*am_12(imat) )
  end do !imat
  ac_12 = dsqrt( ac_12 / rho_12 )

  m_l = u_l/ac_12
  m_r = u_r/ac_12

  ! all-speed scaling:
  !mbar2 = c05 * (vn_l*vn_l + vn_r*vn_r)/(ac_12*ac_12)
  !umag_fs = sqrt(g_freestream.u*g_freestream.u + g_freestream.v*g_freestream.v)
  !m_0   = sqrt(min(1.0, max(mbar2, umag_fs/ac_12)))
  f_a = 1.0 !m_0 * (2.0 - m_0)

  ! split mach number functions
  call splitmach_as(f_a,m_l,msplus_l,msminu_l,psplus_l,psminu_l)
  call splitmach_as(f_a,m_r,msplus_r,msminu_r,psplus_r,psminu_r)

  ! "p"
  temp  = 1.0 - (0.5*(u_l*u_l + u_r*u_r)/(ac_12*ac_12))
  m_p   = -k_p* (max(temp,0.0))* (p_r-p_l) / (f_a* rho_12*ac_12*ac_12)
  m_12 = msplus_l(3) + msminu_r(3) + m_p

  ! "u"
  p_u   = -k_u* psplus_l* psminu_r* f_a* rho_12* ac_12* (u_r-u_l)
  p_12 = psplus_l*p_l + psminu_r*p_r + p_u

  lambda = ac_12 * m_12

  ! flux vector splitting
  lambda_plus = 0.5 * (lambda + dabs(lambda))
  lambda_minu = 0.5 * (lambda - dabs(lambda))
  
  do imat = 1,nummat
    flux(imat)               = lambda_plus*al_l(imat)    + lambda_minu*al_r(imat)
    flux(g_mmi%irmin+imat-1) = lambda_plus*arhom_l(imat) + lambda_minu*arhom_r(imat)
    flux(g_mmi%iemin+imat-1) = lambda_plus*hm_l(imat)    + lambda_minu*hm_r(imat)
  end do !imat
  flux(g_mmi%imome) = lambda_plus*rhou_l + lambda_minu*rhou_r + p_12

  lambda_mag = dabs(lambda) + 1.d-16

  lambda_plus = lambda_plus/(lambda_mag) 
  lambda_minu = lambda_minu/(lambda_mag)

end associate

end subroutine ausmplus_mm6eq

!-------------------------------------------------------------------------------

subroutine ausmplus_nonconserv(ul, ur, uavgl, uavgr, lplus, lminu, lmag, &
                               ncnflux)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), &
                      uavgl(g_neqns), uavgr(g_neqns), &
                      lplus, lminu, lmag

integer :: imat
real*8  :: ncnflux(g_neqns,2), &
           uprim_avg(g_neqns), &
           u_conv_l, u_conv_r, uf, p, &
           alpf(g_mmi%nummat), ncnfl, ncnfr

associate (nummat=>g_mmi%nummat)

  do imat = 1,nummat
    alpf(imat) = dabs(lplus)*ul(imat) + dabs(lminu)*ur(imat)
  end do !imat

  ! left element
  call get_uprim_mm6eq(uavgl, uprim_avg)
  u_conv_l = uprim_avg(g_mmi%imome)
  p = 0.0
  do imat = 1,nummat
    p = p + uavgl(imat)*uprim_avg(g_mmi%irmin+imat-1)
  end do !imat
  ncnfl = p * u_conv_l

  ! right element
  call get_uprim_mm6eq(uavgr, uprim_avg)
  u_conv_r = uprim_avg(g_mmi%imome)
  p = 0.0
  do imat = 1,nummat
    p = p + uavgr(imat)*uprim_avg(g_mmi%irmin+imat-1)
  end do !imat
  ncnfr = p * u_conv_r

  uf = lmag*(lplus+lminu)

  do imat = 1,nummat
    ncnflux(imat,1) = - uavgl(imat) * uf
    ncnflux(g_mmi%iemin+imat-1,1) = - alpf(imat) * ncnfl

    ncnflux(imat,2) = - uavgr(imat) * uf
    ncnflux(g_mmi%iemin+imat-1,2) = - alpf(imat) * ncnfr
  end do !imat

end associate

end subroutine ausmplus_nonconserv

!-------------------------------------------------------------------------------

subroutine ausmplus_nonconserv_p1(ul, ur, uavgl, uavgr, lplus, lminu, lmag, &
                                  ncnflux)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), &
                      uavgl(g_neqns), uavgr(g_neqns), &
                      lplus, lminu, lmag

real*8  :: ncnflux(g_neqns,2), &
           uprim(g_neqns), &
           u_conv_l, u_conv_r, uf, p, &
           alp1f, alp2f, ncnfl, ncnfr

  alp1f = dabs(lplus)*ul(1) + dabs(lminu)*ur(1)
  alp2f = dabs(lplus)*(1.0-ul(1)) + dabs(lminu)*(1.0-ur(1))

  ! left element
  call get_uprim_mm6eq(ul, uprim)
  u_conv_l = uprim(4)
  p = ul(1)*uprim(2) + (1.0-ul(1))*uprim(3)
  ncnfl = p * u_conv_l

  ! right element
  call get_uprim_mm6eq(ur, uprim)
  u_conv_r = uprim(4)
  p = ur(1)*uprim(2) + (1.0-ur(1))*uprim(3)
  ncnfr = p * u_conv_r

  uf = lmag*(lplus+lminu)

  ncnflux(1,1) = - uavgl(1) * uf
  ncnflux(5,1) = - (alp1f-uavgl(1)) * ncnfl
  ncnflux(6,1) = - (alp2f-(1.0-uavgl(1))) * ncnfl

  ncnflux(1,2) = - uavgr(1) * uf
  ncnflux(5,2) = - (alp1f-uavgr(1)) * ncnfr
  ncnflux(6,2) = - (alp2f-(1.0-uavgr(1))) * ncnfr

end subroutine ausmplus_nonconserv_p1

!-------------------------------------------------------------------------------
!----- Split Mach polynomials for AUSM+UP:
!-------------------------------------------------------------------------------
subroutine splitmach_as(fa, mach, msplus, msminu, psplus, psminu)

real*8, intent(in) :: fa, mach

real*8             :: alph_fa, msplus(3), msminu(3), psplus, psminu

    msplus(1) = 0.5d0*(mach + dabs(mach))
    msminu(1) = 0.5d0*(mach - dabs(mach))

    msplus(2) = +0.25d0*(mach + 1.0d0)*(mach + 1.0d0)
    msminu(2) = -0.25d0*(mach - 1.0d0)*(mach - 1.0d0)

    alph_fa = (3.d0/16.d0) * (-4.d0 + 5.d0*fa*fa)

    if (dabs(mach) .ge. 1.d0) then
    
        msplus(3) = msplus(1)
        msminu(3) = msminu(1)
        psplus    = msplus(1)/mach
        psminu    = msminu(1)/mach
   
    else
    
        msplus(3) = msplus(2)* (1.d0 - 2.d0*msminu(2))
        msminu(3) = msminu(2)* (1.d0 + 2.d0*msplus(2))
        psplus    = msplus(2)* &
                    ((+2.d0 - mach) - (16.d0 * alph_fa)*mach*msminu(2))
        psminu    = msminu(2)* &
                    ((-2.d0 - mach) + (16.d0 * alph_fa)*mach*msplus(2))

    end if

end subroutine splitmach_as

!-------------------------------------------------------------------------------
!----- Boundary conditions:
!-------------------------------------------------------------------------------

subroutine get_bc_mm6eq(ucons)

integer :: imat
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)
real*8  :: pmat, u_conv

associate (nummat=>g_mmi%nummat)

  !----- left boundary

  if (g_lbflag .eq. -1) then
     !--- exact inlet
     ucons(1,:,0) = gaussian(coord(1),g_time*a_nd)

  else if (g_lbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     ucons(1,:,0) = ucons(1,:,1)

  else if (g_lbflag .eq. 1) then
     !--- supersonic inflow
     do imat = 1,nummat
        ucons(1,imat,0) = alpha_fs(imat)
        ucons(1,g_mmi%irmin+imat-1,0) = alpha_fs(imat) * rhomat_fs(imat)
        ucons(1,g_mmi%iemin+imat-1,0) = alpha_fs(imat) * &
          eos3_rhoe(g_gam(imat), g_pc(imat), pr_fs, rhomat_fs(imat), u_fs)
     end do !imat
     ucons(1,g_mmi%imome,0) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,0)) * u_fs

  else if (g_lbflag .eq. 2) then
     !--- subsonic inflow
     u_conv = ucons(1,g_mmi%imome,1)/sum(ucons(1,g_mmi%irmin:g_mmi%irmax,1))
     do imat = 1,nummat
        ucons(1,imat,0) = alpha_fs(imat)
        ucons(1,g_mmi%irmin+imat-1,0) = alpha_fs(imat) * rhomat_fs(imat)
        pmat = eos3_pr(g_gam(imat), g_pc(imat), &
                       ucons(1,g_mmi%irmin+imat-1,1)/ucons(1,imat,1), &
                       ucons(1,g_mmi%iemin+imat-1,1)/ucons(1,imat,1), u_conv)
        ucons(1,g_mmi%iemin+imat-1,0) = alpha_fs(imat) * &
          eos3_rhoe(g_gam(imat), g_pc(imat), pmat, rhomat_fs(imat), u_fs)
     end do !imat
     ucons(1,g_mmi%imome,0) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,0)) * u_fs

  else if (g_lbflag .eq. 3) then
     !--- periodic
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,0) = ucons(1,:,imax)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,0) = ucons(1,:,imax) + ucons(2,:,imax)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,0) = ucons(1,:,imax) + ucons(2,:,imax) + ucons(3,:,imax)/3.0
     end if

  else
     write(*,*) "BC-type not set for flag ", g_lbflag

  end if

  !----- right boundary

  if (g_rbflag .eq. -1) then
     !--- exact inlet
     ucons(1,:,imax+1) = gaussian(coord(imax+1),g_time*a_nd)

  else if (g_rbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     ucons(1,:,imax+1) = ucons(1,:,imax)

  else if (g_rbflag .eq. 1) then
     !--- supersonic inflow
     do imat = 1,nummat
        ucons(1,imat,imax+1) = alpha_fs(imat)
        ucons(1,g_mmi%irmin+imat-1,imax+1) = alpha_fs(imat) * rhomat_fs(imat)
        ucons(1,g_mmi%iemin+imat-1,imax+1) = alpha_fs(imat) * &
          eos3_rhoe(g_gam(imat), g_pc(imat), pr_fs, rhomat_fs(imat), u_fs)
     end do !imat
     ucons(1,g_mmi%imome,imax+1) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax+1)) * u_fs

  else if (g_rbflag .eq. 2) then
     !--- subsonic inflow
     u_conv = ucons(1,g_mmi%imome,imax)/sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax))
     do imat = 1,nummat
        ucons(1,imat,imax+1) = alpha_fs(imat)
        ucons(1,g_mmi%irmin+imat-1,imax+1) = alpha_fs(imat) * rhomat_fs(imat)
        pmat = eos3_pr(g_gam(imat), g_pc(imat), &
                       ucons(1,g_mmi%irmin+imat-1,imax)/ucons(1,imat,imax), &
                       ucons(1,g_mmi%iemin+imat-1,imax)/ucons(1,imat,imax), u_conv)
        ucons(1,g_mmi%iemin+imat-1,imax+1) = alpha_fs(imat) * &
          eos3_rhoe(g_gam(imat), g_pc(imat), pmat, rhomat_fs(imat), u_fs)
     end do !imat
     ucons(1,g_mmi%imome,imax+1) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax+1)) * u_fs

  else if (g_rbflag .eq. 3) then
     !--- periodic
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,imax+1) = ucons(1,:,1)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,imax+1) = ucons(1,:,1) - ucons(2,:,1)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,imax+1) = ucons(1,:,1) - ucons(2,:,1) + ucons(3,:,1)/3.0
     end if

  else
     write(*,*) "BC-type not set for flag ", g_rbflag

  end if

  !----- high-order dofs
  if (g_nsdiscr .ge. 1) then
    ucons(2,:,0) = 0.0
    ucons(2,:,imax+1) = 0.0
      if (g_nsdiscr .ge. 12) then
        ucons(3,:,0) = 0.0
        ucons(3,:,imax+1) = 0.0
      end if
  end if

end associate

end subroutine get_bc_mm6eq

!-------------------------------------------------------------------------------

subroutine relaxpressure_p0(ulim, rhsel)

real*8, intent(in) :: ulim(g_tdof,g_neqns,0:imax+1)

integer :: ie, imat
real*8  :: dx, p_star, rel_time, &
           rho, p, aimat, nume, deno
real*8  :: u(g_neqns), up(g_neqns), rhsel(g_gdof,g_neqns,imax)

real*8, dimension(g_mmi%nummat) :: rhom, pm, km, s_alp

associate (nummat=>g_mmi%nummat)

  do ie = 1,imax

    dx = coord(ie+1)-coord(ie)
    u(:) = ulim(1,:,ie)

    call get_uprim_mm6eq(u, up)

    rho  = sum(u(g_mmi%irmin:g_mmi%irmax))
    p = 0.0
    do imat = 1,nummat
      rhom(imat) = u(g_mmi%irmin+imat-1) / u(imat)
      pm(imat)   = up(g_mmi%irmin+imat-1)
      p = p + u(imat)*pm(imat)
    end do !imat

    ! relaxed pressure calculations
    rel_time = 0.0
    nume = 0.0
    deno = 0.0
    do imat = 1,nummat
      aimat = eos3_ss(g_gam(imat), g_pc(imat), rhom(imat), pm(imat))
      km(imat) = rhom(imat) * aimat*aimat
      rel_time = max( rel_time, g_prelct * dx/aimat )
      nume = nume + pm(imat)*u(imat)/km(imat)
      deno = deno +          u(imat)/km(imat)
    end do !imat
    p_star = nume/deno

    do imat = 1,nummat
      s_alp(imat) = 1.0/rel_time * (pm(imat)-p_star)*(u(imat)/km(imat))
    end do !imat

    do imat = 1,nummat
      rhsel(1,imat,ie) = rhsel(1,imat,ie) + dx * s_alp(imat)
      rhsel(1,g_mmi%iemin+imat-1,ie) = rhsel(1,g_mmi%iemin+imat-1,ie) - dx * p*s_alp(imat)
    end do !imat

  end do !ie

end associate

end subroutine relaxpressure_p0

!-------------------------------------------------------------------------------

subroutine relaxpressure_p1(ulim, rhsel)

real*8, intent(in) :: ulim(g_tdof,g_neqns,0:imax+1)

integer :: ig, ie, ieqn, ngauss, imat
data       ngauss/2/
real*8  :: dx, dx2, p_star, rel_time, &
           rho, p, aimat, nume, deno, &
           carea(2), weight(2)
real*8  :: u(g_neqns), up(g_neqns), rhsel(g_gdof,g_neqns,imax)

real*8, dimension(g_mmi%nummat) :: rhom, pm, km, s_alp

associate (nummat=>g_mmi%nummat)

  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax
  do ig = 1,ngauss

    dx = coord(ie+1)-coord(ie)
    dx2 = weight(ig)/2.0 * dx

    ! basis function

    do ieqn = 1,g_neqns
      u(ieqn) = ulim(1,ieqn,ie) + carea(ig) * ulim(2,ieqn,ie)
    end do !ieqn

    call get_uprim_mm6eq(u, up)

    rho  = sum(u(g_mmi%irmin:g_mmi%irmax))
    p = 0.0
    do imat = 1,nummat
      rhom(imat) = u(g_mmi%irmin+imat-1) / u(imat)
      pm(imat)   = up(g_mmi%irmin+imat-1)
      p = p + u(imat)*pm(imat)
    end do !imat

    ! relaxed pressure calculations
    rel_time = 0.0
    nume = 0.0
    deno = 0.0
    do imat = 1,nummat
      aimat = eos3_ss(g_gam(imat), g_pc(imat), rhom(imat), pm(imat))
      km(imat) = rhom(imat) * aimat*aimat
      rel_time = max( rel_time, g_prelct * dx/aimat )
      nume = nume + pm(imat)*u(imat)/km(imat)
      deno = deno +          u(imat)/km(imat)
    end do !imat
    p_star = nume/deno

    do imat = 1,nummat
      s_alp(imat) = 1.0/rel_time * (pm(imat)-p_star)*(u(imat)/km(imat))
    end do !imat

    do imat = 1,nummat
      rhsel(1,imat,ie) = rhsel(1,imat,ie) + dx2 * s_alp(imat)
      rhsel(1,g_mmi%iemin+imat-1,ie) = rhsel(1,g_mmi%iemin+imat-1,ie) - dx2 * p*s_alp(imat)

      rhsel(2,imat,ie) = rhsel(2,imat,ie) + dx2 * carea(ig) * s_alp(imat)
      rhsel(2,g_mmi%iemin+imat-1,ie) = rhsel(2,g_mmi%iemin+imat-1,ie) - dx2 * carea(ig) * p*s_alp(imat)
    end do !imat

  end do !ig
  end do !ie

end associate

end subroutine relaxpressure_p1

!-------------------------------------------------------------------------------

subroutine relaxpressure_p1p2(ulim, rhsel)

real*8, intent(in) :: ulim(g_tdof,g_neqns,0:imax+1)

integer :: ig, ie, ieqn, ngauss
data       ngauss/2/
real*8  :: dx, dx2, b3, p_star, rel_time, &
           al2, rho, p, rho1, rho2, p1, p2, a1, a2, k1, k2, &
           carea(2), weight(2), &
           s_alp1, s_alp2
real*8  :: u(g_neqns), up(g_neqns), rhsel(g_gdof,g_neqns,imax)

  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax
  do ig = 1,ngauss

    dx = coord(ie+1)-coord(ie)
    dx2 = weight(ig)/2.0 * dx

    ! basis function

    b3 = 0.5*carea(ig)*carea(ig) - 1.0/6.0
    do ieqn = 1,g_neqns
      u(ieqn) = ulim(1,ieqn,ie) + carea(ig) * ulim(2,ieqn,ie) + b3 * ulim(3,ieqn,ie)
    end do !ieqn

    call get_uprim_mm6eq(u, up)
    al2 = 1.0-u(1)

    rho  = u(2) + u(3)
    rho1 = u(2) / u(1)
    rho2 = u(3) / al2
    p1   = up(2)
    p2   = up(3)
    p    = u(1)*p1 + al2*p2

    ! relaxed pressure calculations
    a1 = eos3_ss(g_gam(1), g_pc(1), rho1, p1)
    a2 = eos3_ss(g_gam(2), g_pc(2), rho2, p2)
    k1 = rho1 * a1*a1
    k2 = rho2 * a2*a2
    p_star = ( p1*u(1)/k1 + p2*al2/k2 ) / ( u(1)/k1 + al2/k2 )
    rel_time = g_prelct * max(dx/a1, dx/a2)
    s_alp1 = 1.0/rel_time * (p1-p_star)*(u(1)/k1)
    s_alp2 = 1.0/rel_time * (p2-p_star)*(al2/k2)

    rhsel(1,1,ie) = rhsel(1,1,ie) + dx2 * s_alp1
    rhsel(1,5,ie) = rhsel(1,5,ie) - dx2 * p*s_alp1
    rhsel(1,6,ie) = rhsel(1,6,ie) - dx2 * p*s_alp2

    rhsel(2,1,ie) = rhsel(2,1,ie) + dx2 * carea(ig) * s_alp1
    rhsel(2,5,ie) = rhsel(2,5,ie) - dx2 * carea(ig) * p*s_alp1
    rhsel(2,6,ie) = rhsel(2,6,ie) - dx2 * carea(ig) * p*s_alp2

  end do !ig
  end do !ie

end subroutine relaxpressure_p1p2

!-------------------------------------------------------------------------------
!----- 2-material Gaussian function in volume-fraction
!-------------------------------------------------------------------------------

function gaussian(x,t)

real*8, intent(in)  :: x, t
real*8  :: xc, al1, rho1, rho2, gaussian(g_neqns)

  xc  = 0.25 + u_fs*t
  al1 = (1.0-alpha_fs(1)) * dexp( -(x-xc)*(x-xc)/(2.0 * 0.002) ) + alpha_fs(1)

  rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
  rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)
  gaussian(1) = al1
  gaussian(2) = 1.0-al1
  gaussian(3) = al1 * rho1
  gaussian(4) = (1.0-al1) * rho2
  gaussian(5) = (gaussian(2)+gaussian(3)) * u_fs
  gaussian(6) = al1 * eos3_rhoe(g_gam(1), g_pc(1), pr_fs, rho1, u_fs)
  gaussian(7) = (1.0-al1) * eos3_rhoe(g_gam(2), g_pc(2), pr_fs, rho2, u_fs)

end function

!!----- Old DGP1 surface and volume subroutines.
!!-------------------------------------------------------------------------------
!!----- P1 surface contribution to RHS:
!!-------------------------------------------------------------------------------
!
!subroutine surfaceint_p1(ucons, rgrad, rhsel)
!
!integer :: ifc, iel, ier, ieqn
!real*8  :: ul(g_neqns), ur(g_neqns), uavgl(g_neqns), uavgr(g_neqns), &
!           ncnflux(g_neqns,2), intflux(g_neqns), rhsel(g_gdof,g_neqns,imax), &
!           lplus, lminu, lmag
!
!real*8  :: ucons(g_tdof,g_neqns,0:imax+1)
!
!  do ifc = 1,imax+1
!
!  intflux = 0.0
!  ncnflux = 0.0
!
!  iel = ifc - 1
!  ier = ifc
!
!  do ieqn = 1,g_neqns
!    ul(ieqn) = ucons(1,ieqn,iel) + ucons(2,ieqn,iel)
!    ur(ieqn) = ucons(1,ieqn,ier) - ucons(2,ieqn,ier)
!
!    uavgl(ieqn) = ucons(1,ieqn,iel)
!    uavgr(ieqn) = ucons(1,ieqn,ier)
!  end do !ieqn
!
!  !--- fluxes
!
!  if (i_flux .eq. 2) then
!     call ausmplus_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
!     call ausmplus_nonconserv_p1(ul, ur, uavgl, uavgr, &
!                                 lplus, lminu, lmag, ncnflux)
!  else
!     write(*,*) "Invalid flux scheme."
!     stop
!  endif
!
!  if (iel .gt. 0) then
!    do ieqn = 1,g_neqns
!          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
!          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - ncnflux(ieqn,1)
!          rhsel(2,ieqn,iel) = rhsel(2,ieqn,iel) - intflux(ieqn)
!          rhsel(2,ieqn,iel) = rhsel(2,ieqn,iel) - ncnflux(ieqn,1)
!    end do !ieqn
!  end if
!
!  if (ier .lt. (imax+1)) then
!    do ieqn = 1,g_neqns
!          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
!          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + ncnflux(ieqn,2)
!          rhsel(2,ieqn,ier) = rhsel(2,ieqn,ier) - intflux(ieqn)
!          rhsel(2,ieqn,ier) = rhsel(2,ieqn,ier) - ncnflux(ieqn,2)
!    end do !ieqn
!  end if
!
!  end do !ifc
!
!end subroutine surfaceint_p1
!
!!-------------------------------------------------------------------------------
!!----- P1 volume contribution to RHS:
!!-------------------------------------------------------------------------------
!
!subroutine volumeint_p1(ucons, rgrad, rhsel)
!
!integer :: ig, ie, ieqn, ngauss
!data       ngauss/2/
!
!real*8  :: dx2, p, &
!           u(g_neqns), up(g_neqns), &
!           carea(2), weight(2), &
!           cflux(g_neqns), &
!           rhsel(g_gdof,g_neqns,imax)
!
!real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)
!
!  ngauss = 2
!
!  call rutope(1, ngauss, carea, weight)
!
!  do ie = 1,imax
!  do ig = 1,ngauss
!
!    dx2 = weight(ig) ! <-- 2.0/dx * weight(ig)/2.0 * dx
!
!    ! basis function
!
!    do ieqn = 1,g_neqns
!      u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie)
!    end do !ieqn
!
!    call get_uprim_mm6eq(u, up)
!
!    p = up(1)*up(2) + (1.0-up(1))*up(3)
!
!    ! conservative fluxes
!    cflux(1) = 0.0!up(4) * up(1)
!    cflux(2) = up(4) * u(2)
!    cflux(3) = up(4) * u(3)
!    cflux(4) = up(4) * u(4) + p
!    cflux(5) = up(4) * (u(5) + up(1)*up(2))
!    cflux(6) = up(4) * (u(6) + (1.0-up(1))*up(3))
!
!    do ieqn = 1,g_neqns
!      rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + dx2 * cflux(ieqn)
!    end do !ieqn
!
!  end do !ig
!  end do !ie
!
!end subroutine volumeint_p1

!-------------------------------------------------------------------------------

END MODULE rhs_flux_mm6eq
