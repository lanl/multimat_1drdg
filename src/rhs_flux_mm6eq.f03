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
!----- rDG RHS:
!-------------------------------------------------------------------------------

subroutine rhs_rdg_mm6eq(ucons, uprim, rhsel)

real*8  :: rhsel(g_gdof,g_neqns,imax)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

  call flux_rdg_mm6eq(ucons, uprim, rhsel)

  if (g_nprelx .eq. 1) then
    call relaxpressure_rdg(ucons, rhsel)
  end if

end subroutine rhs_rdg_mm6eq

!-------------------------------------------------------------------------------
!----- rDG Advective-flux contribution to RHS:
!-------------------------------------------------------------------------------

subroutine flux_rdg_mm6eq(ucons, uprim, rhsel)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

real*8  :: riemanngrad(g_mmi%nummat+1,imax), &
           vriemann(g_mmi%nummat+1,imax+1), &
           rhsel(g_gdof,g_neqns,imax)

  riemanngrad = 0.0
  vriemann = 0.0

  !--- surface integration
  call surfaceint_dg(ucons, uprim, riemanngrad, vriemann, rhsel)
  !--- volume integration
  call volumeint_dg(ucons, uprim, riemanngrad, vriemann, rhsel)

end subroutine flux_rdg_mm6eq

!-------------------------------------------------------------------------------
!----- DG surface contribution to RHS:
!-------------------------------------------------------------------------------

subroutine surfaceint_dg(ucons, uprim, rgrad, vriem, rhsel)

integer :: ifc, iel, ier, ieqn, imat
real*8  :: ul(g_neqns), ur(g_neqns), uavgl(g_neqns), uavgr(g_neqns), &
           up_l(g_neqns), up_r(g_neqns), &
           pp_l(g_nprim), pp_r(g_nprim), &
           intflux(g_neqns), rhsel(g_gdof,g_neqns,imax), &
           lplus, lminu, lmag, pplus, pminu, pstar

real*8  :: rgrad(g_mmi%nummat+1,imax), vriem(g_mmi%nummat+1,imax+1)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ifc = 1,imax+1

  intflux = 0.0

  iel = ifc - 1
  ier = ifc

  !--- reconstructed and limited values of conserved variables and
  !--- primitive quantities i.e. material partial pressures (alpha_k * p_k)
  !--- and bulk fluid velocity.

  !--- dgp0
  if (g_nsdiscr .eq. 0) then
    do ieqn = 1,g_neqns
      ul(ieqn) = ucons(1,ieqn,iel)
      ur(ieqn) = ucons(1,ieqn,ier)
    end do !ieqn
    do ieqn = 1,g_nprim
      pp_l(ieqn) = uprim(1,ieqn,iel)
      pp_r(ieqn) = uprim(1,ieqn,ier)
    end do !ieqn

  !--- rdgp0p1 or dgp1
  elseif ((g_nsdiscr .eq. 11) .or. (g_nsdiscr .eq. 1)) then
    do ieqn = 1,g_neqns
      ul(ieqn) = ucons(1,ieqn,iel) + ucons(2,ieqn,iel)
      ur(ieqn) = ucons(1,ieqn,ier) - ucons(2,ieqn,ier)
    end do !ieqn
    do ieqn = 1,g_nprim
      pp_l(ieqn) = uprim(1,ieqn,iel) + uprim(2,ieqn,iel)
      pp_r(ieqn) = uprim(1,ieqn,ier) - uprim(2,ieqn,ier)
    end do !ieqn

  !--- rdgp1p2
  elseif (g_nsdiscr .eq. 12) then
    do ieqn = 1,g_neqns
      ul(ieqn) = ucons(1,ieqn,iel) + ucons(2,ieqn,iel) + 1.0/3.0*ucons(3,ieqn,iel)
      ur(ieqn) = ucons(1,ieqn,ier) - ucons(2,ieqn,ier) + 1.0/3.0*ucons(3,ieqn,ier)
    end do !ieqn
    do ieqn = 1,g_nprim
      pp_l(ieqn) = uprim(1,ieqn,iel) + uprim(2,ieqn,iel) + 1.0/3.0*uprim(3,ieqn,iel)
      pp_r(ieqn) = uprim(1,ieqn,ier) - uprim(2,ieqn,ier) + 1.0/3.0*uprim(3,ieqn,ier)
    end do !ieqn

  end if

  call check_volfrac(iel, ul)
  call check_volfrac(ier, ur)

  uavgl(:) = ucons(1,:,iel)
  uavgr(:) = ucons(1,:,ier)

  call get_uprim_mm6eq(ul, up_l)
  call get_uprim_mm6eq(ur, up_r)

  !--- fluxes

  if (i_flux .eq. 1) then
     call llf_mm6eq(ul, ur, up_l, up_r, pp_l, pp_r, &
                    intflux, lplus, lminu, lmag, pplus, pminu)
  else if (i_flux .eq. 2) then
     call ausmplus_mm6eq(ul, ur, up_l, up_r, pp_l, pp_r, &
                         intflux, lplus, lminu, lmag, pplus, pminu)
  !else if (i_flux .eq. 3) then
  !   call hllc_mm6eq(ul, ur, up_l, up_r, intflux, lplus, lminu, lmag)
  !   pplus = lplus
  !   pminu = lminu
  else
     write(*,*) "Invalid flux scheme."
     stop
  endif

  !--- compute gradients of volume fractions and velocity for the
  !--- non-conservative terms from Riemann reconstructed values
  do imat = 1,nummat
    vriem(imat,ifc) = al_star(lplus, lminu, &
                              ul(imat)*up_l(g_mmi%irmin+imat-1), &
                              ur(imat)*up_r(g_mmi%irmin+imat-1))
                              !pp_l(apr_idx(nummat, imat)), &
                              !pp_r(apr_idx(nummat, imat)))
    !vriem(imat,ifc) = pplus*ul(imat)*up_l(g_mmi%irmin+imat-1) &
    !  + pminu*ur(imat)*up_r(g_mmi%irmin+imat-1)
    !vriem(imat,ifc) = pplus*pp_l(apr_idx(nummat, imat)) &
    !  + pminu*pp_r(apr_idx(nummat, imat))
  end do !imat
  vriem(nummat+1,ifc) = lmag*(lplus+lminu)

  if (iel .gt. 0) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
          if (g_nsdiscr .ge. 11) rhsel(2,ieqn,iel) = rhsel(2,ieqn,iel) - intflux(ieqn)
    end do !ieqn

    do imat = 1,nummat
      rgrad(imat,iel) = rgrad(imat,iel) + vriem(imat,ifc)
    end do !imat
    rgrad(nummat+1,iel) = rgrad(nummat+1,iel) + vriem(nummat+1,ifc)

  end if

  if (ier .lt. (imax+1)) then
    do ieqn = 1,g_neqns
          rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
          if (g_nsdiscr .ge. 11) rhsel(2,ieqn,ier) = rhsel(2,ieqn,ier) - intflux(ieqn)
    end do !ieqn

    do imat = 1,nummat
      rgrad(imat,ier) = rgrad(imat,ier) - vriem(imat,ifc)
    end do !imat
    rgrad(nummat+1,ier) = rgrad(nummat+1,ier) - vriem(nummat+1,ifc)

  end if

  end do !ifc

end associate

end subroutine surfaceint_dg

!-------------------------------------------------------------------------------
!----- DG volume contribution to RHS:
!-------------------------------------------------------------------------------

subroutine volumeint_dg(ucons, uprim, rgrad, vriem, rhsel)

integer :: ig, ie, ieqn, ngauss, imat
data       ngauss/2/

real*8  :: dx2, b3, p, hmat, viriem, &
           u(g_neqns), up(g_neqns), pp(g_nprim), &
           rhob, y(g_mmi%nummat), dapdx, &
           carea(2), weight(2), &
           cflux(g_neqns), &
           nflux(g_gdof,g_neqns), &
           rhsel(g_gdof,g_neqns,imax)

real*8, intent(in) :: rgrad(g_mmi%nummat+1,imax), vriem(g_mmi%nummat+1,imax+1), &
                      ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax
  cflux = 0.0
  nflux = 0.0
  do ig = 1,ngauss

    dx2 = weight(ig) ! <-- 2.0/dx * weight(ig)/2.0 * dx

    !--- dgp0
    if (g_nsdiscr .eq. 0) then
      do ieqn = 1,g_neqns
        u(ieqn) = ucons(1,ieqn,ie)
      end do !ieqn
      do ieqn = 1,g_nprim
        pp(ieqn) = uprim(1,ieqn,ie)
      end do !ieqn

    !--- rdgp0p1 or dgp1
    elseif ((g_nsdiscr .eq. 11) .or. (g_nsdiscr .eq. 1)) then
      do ieqn = 1,g_neqns
        u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie)
      end do !ieqn
      do ieqn = 1,g_nprim
        pp(ieqn) = uprim(1,ieqn,ie) + carea(ig) * uprim(2,ieqn,ie)
      end do !ieqn

    !--- rdgp1p2
    elseif (g_nsdiscr .eq. 12) then
      b3 = 0.5*carea(ig)*carea(ig) - 1.0/6.0
      do ieqn = 1,g_neqns
        u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie) + b3 * ucons(3,ieqn,ie)
      end do !ieqn
      do ieqn = 1,g_nprim
        pp(ieqn) = uprim(1,ieqn,ie) + carea(ig) * uprim(2,ieqn,ie) + b3 * uprim(3,ieqn,ie)
      end do !ieqn

    end if

    call check_volfrac(ie, u)

    viriem = 0.5* (vriem(nummat+1,ie) + vriem(nummat+1,ie+1)) + carea(ig) * rgrad(nummat+1,ie)/2.0

    call get_uprim_mm6eq(u, up)

    p = 0.0
    dapdx = 0.0
    rhob = sum(u(g_mmi%irmin:g_mmi%irmax))
    do imat = 1,nummat
      p = p + u(imat)*up(g_mmi%irmin+imat-1)
      !p = p + pp(apr_idx(nummat, imat))
      dapdx = dapdx + rgrad(imat,ie)
      y(imat) = u(g_mmi%irmin+imat-1) / rhob
    end do !imat

    !--- flux terms
    ! momentum flux
    cflux(g_mmi%imome) = pp(vel_idx(nummat, 0)) &
      * sum(pp(mmom_idx(nummat,1):mmom_idx(nummat,nummat))) + p
    !cflux(g_mmi%imome) = up(g_mmi%imome) &
    !  * sum(pp(mmom_idx(nummat,1):mmom_idx(nummat,nummat))) + p
    nflux(1,g_mmi%imome) = 0.0
    if (g_nsdiscr .ge. 11) nflux(2,g_mmi%imome) = 0.0
    do imat = 1,nummat
      !hmat = u(g_mmi%iemin+imat-1) + pp(apr_idx(nummat, imat))
      hmat = u(g_mmi%iemin+imat-1) + u(imat)*up(g_mmi%irmin+imat-1)
      ! other conservative fluxes
      cflux(imat) = 0.0
      cflux(g_mmi%irmin+imat-1) = pp(vel_idx(nummat, 0)) * u(g_mmi%irmin+imat-1)
      cflux(g_mmi%iemin+imat-1) = pp(vel_idx(nummat, 0)) * hmat
      !cflux(g_mmi%irmin+imat-1) = up(g_mmi%imome) * u(g_mmi%irmin+imat-1)
      !cflux(g_mmi%iemin+imat-1) = up(g_mmi%imome) * hmat

      ! non-conservative fluxes
      nflux(1,imat) = u(imat) * rgrad(nummat+1,ie)
      nflux(1,g_mmi%irmin+imat-1) = 0.0
      nflux(1,g_mmi%iemin+imat-1) = - pp(vel_idx(nummat, 0)) * ( y(imat) * dapdx &
                                                         - rgrad(imat,ie) )
      !nflux(1,g_mmi%iemin+imat-1) = - up(g_mmi%imome) * ( y(imat) * dapdx &
      !                                                   - rgrad(imat,ie) )
      if (g_nsdiscr .ge. 11) then
        nflux(2,imat) = nflux(1,imat) * carea(ig) + u(imat) * viriem * 2.0 !pp(vel_idx(nummat, 0)) * 2.0
        nflux(2,g_mmi%irmin+imat-1) = 0.0
        nflux(2,g_mmi%iemin+imat-1) = nflux(1,g_mmi%iemin+imat-1) * carea(ig)
      end if
    end do !imat

    do ieqn = 1,g_neqns
      if (g_nsdiscr .ge. 11) rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + dx2 * cflux(ieqn)
    end do !ieqn

    do ieqn = 1,g_neqns
      rhsel(1,ieqn,ie) = rhsel(1,ieqn,ie) + 0.5 * dx2 * nflux(1,ieqn)
      if (g_nsdiscr .ge. 11) rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + 0.5 * dx2 * nflux(2,ieqn)
    end do !ieqn

  end do !ig
  end do !ie

end associate

end subroutine volumeint_dg

!-------------------------------------------------------------------------------
!----- 2fluid Lax-Friedrichs flux:
!-------------------------------------------------------------------------------

subroutine llf_mm6eq(ul, ur, up_l, up_r, pp_l, pp_r, &
                     flux, lplus, lminu, lmag, pplus, pminu)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), up_l(g_neqns), up_r(g_neqns), &
                      pp_l(g_nprim), pp_r(g_nprim)

integer :: imat
real*8 :: flux(g_neqns), lplus, lminu, lmag
real*8 :: ffunc_l(g_neqns), ffunc_r(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,rhom_l,em_l,am_l,hm_l,pm_l, &
                                   arhom_r,rhom_r,em_r,am_r,hm_r,pm_r, &
                                   am_12,rhom_12,al_12

real*8 :: rhou_l, u_l, rho_l, pi_l, p_l, &
          rhou_r, u_r, rho_r, pi_r, p_r
real*8 :: rho_12,ac_12,pplus,pminu
  
real*8 :: lambda

associate (nummat=>g_mmi%nummat)

  flux(:) = 0.0

  ! material and bulk states
  p_l = 0.0
  p_r = 0.0
  rhou_l = 0.0
  rhou_r = 0.0
  do imat = 1,nummat
  ! ul
    al_l(imat)    = ul(imat)
    arhom_l(imat) = ul(g_mmi%irmin+imat-1)
    em_l(imat)    = ul(g_mmi%iemin+imat-1)

    rhom_l(imat) = arhom_l(imat) / al_l(imat)
    pi_l         = up_l(g_mmi%irmin+imat-1)
    !pi_l         = pp_l(apr_idx(nummat, imat)) / al_l(imat)
    am_l(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_l(imat), al_l(imat), pi_l)
    hm_l(imat)   = em_l(imat) + al_l(imat)*pi_l
    p_l = p_l + al_l(imat)*pi_l
    rhou_l = rhou_l + pp_l(mmom_idx(nummat, imat))
    pm_l(imat)   = pi_l

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r(imat)    = ur(g_mmi%iemin+imat-1)

    rhom_r(imat) = arhom_r(imat) / al_r(imat)
    pi_r         = up_r(g_mmi%irmin+imat-1)
    !pi_r         = pp_r(apr_idx(nummat, imat)) / al_r(imat)
    am_r(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_r(imat), al_r(imat), pi_r)
    hm_r(imat)   = em_r(imat) + al_r(imat)*pi_r
    p_r = p_r + al_r(imat)*pi_r
    rhou_r = rhou_r + pp_r(mmom_idx(nummat, imat))
    pm_r(imat)   = pi_r
  end do !imat

  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l = pp_l(vel_idx(nummat, 0))
  u_r = pp_r(vel_idx(nummat, 0))

  ! average states
  rho_12  = 0.5*(rho_l + rho_r)
  do imat = 1,nummat
    rhom_12(imat) = 0.5*(rhom_l(imat) + rhom_r(imat))
    am_12(imat)   = 0.5*(am_l(imat) + am_r(imat))
    al_12(imat)   = 0.5*(al_l(imat) + al_r(imat))
  end do !imat

  ! numerical speed of sound choice:

  ! Kapila
  ac_12 = 0.0
  do imat = 1,nummat
    ac_12 = ac_12 + ( al_12(imat)*rhom_12(imat)*am_12(imat)*am_12(imat) )
  end do !imat
  ac_12 = dsqrt( ac_12 / rho_12 )

  lambda = ac_12 + max(dabs(u_l),dabs(u_r))

  ! flux functions
  ffunc_l = 0.0
  ffunc_r = 0.0
  do imat = 1,nummat
    ffunc_l(imat) = u_l * al_l(imat)
    ffunc_l(g_mmi%irmin+imat-1) = u_l * arhom_l(imat)
    ffunc_l(g_mmi%iemin+imat-1) = u_l * hm_l(imat)
    ffunc_l(g_mmi%imome) = ffunc_l(g_mmi%imome) &
      + u_l * pp_l(mmom_idx(nummat, imat)) + al_l(imat)*pm_l(imat)

    ffunc_r(imat) = u_r * al_r(imat)
    ffunc_r(g_mmi%irmin+imat-1) = u_r * arhom_r(imat)
    ffunc_r(g_mmi%iemin+imat-1) = u_r * hm_r(imat)
    ffunc_r(g_mmi%imome) = ffunc_r(g_mmi%imome) &
      + u_r * pp_r(mmom_idx(nummat, imat)) + al_r(imat)*pm_r(imat)
  end do !imat

  flux = 0.5 * ( ffunc_l+ffunc_r - lambda*(ur-ul) )

  lplus = 0.5 * (u_l + lambda)
  lminu = 0.5 * (u_r - lambda)
  pplus = 0.5
  pminu = 0.5

  lmag = dabs(lplus+lminu) + 1.d-16
  lplus = lplus/lmag
  lminu = lminu/lmag

end associate

end subroutine llf_mm6eq

!-------------------------------------------------------------------------------
!----- 2fluid AUSM+UP:
!-------------------------------------------------------------------------------

subroutine ausmplus_mm6eq(ul, ur, up_l, up_r, pp_l, pp_r, &
  flux, lambda_plus, lambda_minu, lambda_mag, psplus_l, psminu_r)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), up_l(g_neqns), up_r(g_neqns), &
                      pp_l(g_nprim), pp_r(g_nprim)

integer :: imat
real*8 :: flux(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,rhom_l,em_l,am_l,hm_l,pm_l, &
                                   arhom_r,rhom_r,em_r,am_r,hm_r,pm_r, &
                                   am_12,rhom_12,al_12

real*8 :: rhou_l, u_l, m_l, rho_l, pi_l, p_l, &
          rhou_r, u_r, m_r, rho_r, pi_r, p_r
real*8 :: rho_12, ac_12, &
          f_a, m_12, p_12, m_p, p_u(g_mmi%nummat)
real*8 :: msplus_l(3),msplus_r(3),msminu_l(3),msminu_r(3)
real*8 :: psplus_l,psplus_r,psminu_l,psminu_r,pplus,pminu
real*8 :: temp

real*8 :: lambda,lambda_plus, lambda_minu, lambda_mag

real*8 :: k_p, k_u
!real*8 :: mbar2, umag_fs, m_0

associate (nummat=>g_mmi%nummat)

  k_p = 0.5;
  k_u = 0.5;

  flux(:) = 0.0

  ! material and bulk states
  p_l = 0.0
  p_r = 0.0
  rhou_l = 0.0
  rhou_r = 0.0
  do imat = 1,nummat
  ! ul
    al_l(imat)    = ul(imat)
    arhom_l(imat) = ul(g_mmi%irmin+imat-1)
    em_l(imat)    = ul(g_mmi%iemin+imat-1)

    rhom_l(imat) = arhom_l(imat) / al_l(imat)
    pi_l         = up_l(g_mmi%irmin+imat-1)
    !pi_l         = pp_l(apr_idx(nummat, imat)) / al_l(imat)
    am_l(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_l(imat), al_l(imat), pi_l)
    hm_l(imat)   = em_l(imat) + al_l(imat)*pi_l
    p_l = p_l + al_l(imat)*pi_l
    rhou_l = rhou_l + pp_l(mmom_idx(nummat, imat))
    pm_l(imat)   = pi_l

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r(imat)    = ur(g_mmi%iemin+imat-1)

    rhom_r(imat) = arhom_r(imat) / al_r(imat)
    pi_r         = up_r(g_mmi%irmin+imat-1)
    !pi_r         = pp_r(apr_idx(nummat, imat)) / al_r(imat)
    am_r(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_r(imat), al_r(imat), pi_r)
    hm_r(imat)   = em_r(imat) + al_r(imat)*pi_r
    p_r = p_r + al_r(imat)*pi_r
    rhou_r = rhou_r + pp_r(mmom_idx(nummat, imat))
    pm_r(imat)   = pi_r
  end do !imat

  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l = pp_l(vel_idx(nummat, 0))
  u_r = pp_r(vel_idx(nummat, 0))
  !u_l = rhou_l/rho_l
  !u_r = rhou_r/rho_r

  ! average states
  rho_12  = 0.5*(rho_l + rho_r)
  do imat = 1,nummat
    rhom_12(imat) = 0.5*(rhom_l(imat) + rhom_r(imat))
    am_12(imat)   = 0.5*(am_l(imat) + am_r(imat))
    al_12(imat)   = 0.5*(al_l(imat) + al_r(imat))
  end do !imat

  ! numerical speed of sound choice:

  ! Kapila
  ac_12 = 0.0
  do imat = 1,nummat
    ac_12 = ac_12 + ( al_12(imat)*rhom_12(imat)*am_12(imat)*am_12(imat) )
  end do !imat
  ac_12 = dsqrt( ac_12 / rho_12 )

  ! Wood
  !ac_12 = 0.0
  !do imat = 1,nummat
  !  ac_12 = ac_12 + al_12(imat) / ( rho_12*rhom_12(imat)*am_12(imat)*am_12(imat) )
  !end do !imat
  !ac_12 = dsqrt( 1.0 / ac_12 ) / rho_12

  m_l = u_l/ac_12
  m_r = u_r/ac_12

  ! all-speed scaling:
  !mbar2 = 0.5 * (u_l*u_l + u_r*u_r)/(ac_12*ac_12)
  !umag_fs = u_fs*u_fs/(ac_12*ac_12)
  !m_0   = dsqrt(min(1.0, max(mbar2, umag_fs)))
  f_a = 1.0 !m_0 * (2.0 - m_0)

  ! split mach number functions
  call splitmach_as(f_a,m_l,msplus_l,msminu_l,psplus_l,psminu_l)
  call splitmach_as(f_a,m_r,msplus_r,msminu_r,psplus_r,psminu_r)

  ! "p"
  temp  = 1.0 - (0.5*(u_l*u_l + u_r*u_r)/(ac_12*ac_12))
  m_p   = -k_p* (max(temp,0.0))* (p_r-p_l) / (f_a* rho_12*ac_12*ac_12)
  m_12 = msplus_l(3) + msminu_r(3) + m_p

  ! "u"
  !p_u   = -k_u* psplus_l* psminu_r* f_a* rho_12* ac_12* (u_r-u_l)
  !p_12 = psplus_l*p_l + psminu_r*p_r + p_u
  do imat = 1,nummat
    p_u(imat) = -k_u* psplus_l* psminu_r* f_a* 0.5*(arhom_l(imat)+arhom_r(imat)) &
      * ac_12* (u_r-u_l)
  end do !imat
  p_12 = psplus_l*p_l + psminu_r*p_r + sum(p_u)

  lambda = ac_12 * m_12

  ! flux vector splitting
  lambda_plus = 0.5 * (lambda + dabs(lambda))
  lambda_minu = 0.5 * (lambda - dabs(lambda))

  pplus = psplus_l - k_u* psplus_l* psminu_r* f_a* (u_r-u_l) / ac_12 !0.5
  pminu = psminu_r - k_u* psplus_l* psminu_r* f_a* (u_r-u_l) / ac_12 !0.5

  !psplus_l = pplus
  !psminu_r = pminu

  do imat = 1,nummat
    flux(imat)               = lambda_plus*al_l(imat)    + lambda_minu*al_r(imat)
    flux(g_mmi%irmin+imat-1) = lambda_plus*arhom_l(imat) + lambda_minu*arhom_r(imat)
    flux(g_mmi%iemin+imat-1) = lambda_plus*hm_l(imat)    + lambda_minu*hm_r(imat)
    !flux(g_mmi%iemin+imat-1) = lambda_plus*em_l(imat) + lambda_minu*em_r(imat) &
    !  + psplus_l*al_l(imat)*pm_l(imat)*u_l &
    !  + psminu_r*al_r(imat)*pm_r(imat)*u_r
    flux(g_mmi%imome) = flux(g_mmi%imome) &
      + lambda_plus*pp_l(mmom_idx(nummat, imat)) &
      + lambda_minu*pp_r(mmom_idx(nummat, imat)) &
      + psplus_l*al_l(imat)*pm_l(imat) &
      + psminu_r*al_r(imat)*pm_r(imat) + p_u(imat)
  end do !imat
  !flux(g_mmi%imome) = lambda_plus*rhou_l + lambda_minu*rhou_r + p_12

  lambda_mag = dabs(lambda) + 1.d-16

  lambda_plus = lambda_plus/(lambda_mag)
  lambda_minu = lambda_minu/(lambda_mag)

end associate

end subroutine ausmplus_mm6eq

!------------------------------------------------------------------------------
!----- Numerical flux by HLLC:
!------------------------------------------------------------------------------

subroutine hllc_mm6eq(ul, ur, up_l, up_r, flux, lambda_plus, lambda_minu, lambda_mag)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), up_l(g_neqns), up_r(g_neqns)

integer :: imat
real*8 :: flux(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,rhom_l,em_l,am_l,hm_l, &
                                   arhom_r,rhom_r,em_r,am_r,hm_r

real*8 :: rhou_l, u_l, rho_l, pi_l, p_l, ac_l, &
          rhou_r, u_r, rho_r, pi_r, p_r, ac_r

real*8  :: rij, rij1, vroe, croe, &
           si,sj,sm,pstar, &
           frac,UstarI(g_neqns),UstarJ(g_neqns), &
           tmp1,tmp2

real*8 :: lambda,lambda_plus, lambda_minu, lambda_mag

associate (nummat=>g_mmi%nummat)

  flux(:) = 0.0

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
    !pi_l         = pp_l(apr_idx(nummat, imat)) / al_l(imat)
    am_l(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_l(imat), al_l(imat), pi_l)
    hm_l(imat)   = em_l(imat) + al_l(imat)*pi_l
    p_l = p_l + al_l(imat)*pi_l

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r(imat)    = ur(g_mmi%iemin+imat-1)

    rhom_r(imat) = arhom_r(imat) / al_r(imat)
    pi_r         = up_r(g_mmi%irmin+imat-1)
    !pi_r         = pp_r(apr_idx(nummat, imat)) / al_r(imat)
    am_r(imat)   = eos3_ss(g_gam(imat), g_pc(imat), rhom_r(imat), al_r(imat), pi_r)
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

  ! numerical speed of sound choice:
  ac_l = 0.0
  ac_r = 0.0
  do imat = 1,nummat
    ac_l = ac_l + ( al_l(imat)*rhom_l(imat)*am_l(imat)*am_l(imat) )
    ac_r = ac_r + ( al_r(imat)*rhom_r(imat)*am_r(imat)*am_r(imat) )
  end do !imat
  ac_l = dsqrt( ac_l / rho_l )
  ac_r = dsqrt( ac_r / rho_r )

  ! Roe-averaged variables
  rij = dsqrt(rho_r/rho_l)
  rij1 = rij + 1.d0
  vroe = (rij*u_r + u_l)/rij1
  croe = (rij*ac_r + ac_l)/rij1

  ! signal velocities
  si = min((u_l - ac_l), (vroe - croe))
  sj = max((u_r + ac_r), (vroe + croe))
  tmp1 = rho_r*u_r*(sj-u_r) - rho_l*u_l*(si-u_l) + p_l - p_r
  tmp2 = rho_r*(sj-u_r) - rho_l*(si-u_l)
  sm = tmp1/tmp2
  pstar = rho_l*(u_l-si)*(u_l-sm) + p_l

  ! UstarI
  frac = 1.d0 / (si - sm)
  do imat = 1,nummat
  UstarI(imat)               = frac* (si-u_l) * al_l(imat)
  UstarI(g_mmi%irmin+imat-1) = frac* (si-u_l) * arhom_l(imat)
  UstarI(g_mmi%iemin+imat-1) = frac* ((si-u_l) * em_l(imat) &
                                      - al_l(imat)*up_l(g_mmi%irmin+imat-1)*u_l &
                                      + al_l(imat)*pstar*sm )
  end do !imat
  UstarI(g_mmi%imome) = frac* ((si-u_l)* (rhou_l) + (pstar - p_l))

  ! UstarJ
  frac = 1.d0 / (sj - sm)
  do imat = 1,nummat
  UstarJ(imat)               = frac* (sj-u_r) * al_r(imat)
  UstarJ(g_mmi%irmin+imat-1) = frac* (sj-u_r) * arhom_r(imat)
  UstarJ(g_mmi%iemin+imat-1) = frac* ((sj-u_r) * em_r(imat) &
                                      - al_r(imat)*up_r(g_mmi%irmin+imat-1)*u_r &
                                      + al_r(imat)*pstar*sm )
  end do !imat
  UstarJ(g_mmi%imome) = frac* ((sj-u_r)* (rhou_r) + (pstar - p_r))

  ! Flux
  if (si > 0.d0) then
    do imat = 1,nummat
    flux(imat)               = u_l * al_l(imat)
    flux(g_mmi%irmin+imat-1) = u_l * arhom_l(imat)
    flux(g_mmi%iemin+imat-1) = u_l * hm_l(imat)
    end do !imat
    flux(g_mmi%imome) = u_l * rhou_l + p_l

    lambda = u_l

  else if ( (si <= 0.d0) .and. (sm > 0.d0)) then
    do imat = 1,nummat
    flux(imat)               = sm * UstarI(imat)
    flux(g_mmi%irmin+imat-1) = sm * UstarI(g_mmi%irmin+imat-1)
    flux(g_mmi%iemin+imat-1) = sm * UstarI(g_mmi%iemin+imat-1) + sm*pstar
    end do !imat
    flux(g_mmi%imome) = sm * UstarI(g_mmi%imome) + pstar

    lambda = sm

  else if ((sm <= 0.d0) .and. (sj >= 0.d0)) then
    do imat = 1,nummat
    flux(imat)               = sm * UstarJ(imat)
    flux(g_mmi%irmin+imat-1) = sm * UstarJ(g_mmi%irmin+imat-1)
    flux(g_mmi%iemin+imat-1) = sm * UstarJ(g_mmi%iemin+imat-1) + sm*pstar
    end do !imat
    flux(g_mmi%imome) = sm * UstarJ(g_mmi%imome) + pstar

    lambda = sm

  else ! sj.lt.0.d0
    do imat = 1,nummat
    flux(imat)               = u_r * al_r(imat)
    flux(g_mmi%irmin+imat-1) = u_r * arhom_r(imat)
    flux(g_mmi%iemin+imat-1) = u_r * hm_r(imat)
    end do !imat
    flux(g_mmi%imome) = u_r * rhou_r + p_r

    lambda = u_r

  end if

  !--- these Riemann velocity estimates need to be improved so that
  !    well-balancedness is maintained
  lambda_plus = 0.5 * (lambda + dabs(lambda))
  lambda_minu = 0.5 * (lambda - dabs(lambda))

  lambda_mag = dabs(lambda) + 1.d-16

  lambda_plus = lambda_plus/(lambda_mag)
  lambda_minu = lambda_minu/(lambda_mag)

end associate

end subroutine hllc_mm6eq

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

subroutine get_bc_mm6eq(ucons, uprim)

integer :: imat
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)
real*8  :: pmat, u_conv

associate (nummat=>g_mmi%nummat)

  !----- left boundary

  if (g_lbflag .eq. -1) then
     !--- exact inlet
     ucons(1,:,0) = gaussian(coord(1),g_time*a_nd)
     uprim(1,apr_idx(nummat,1):apr_idx(nummat,nummat),0) = ucons(1,1:nummat,0) * pr_fs
     uprim(1,mmom_idx(nummat,1):mmom_idx(nummat,nummat),0) = &
       ucons(1,g_mmi%irmin:g_mmi%irmax,0) * u_fs
     uprim(1,vel_idx(nummat, 0),0) = u_fs

  else if (g_lbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,0) = ucons(1,:,1)
       uprim(1,:,0) = uprim(1,:,1)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,0) = ucons(1,:,1) - ucons(2,:,1)
       uprim(1,:,0) = uprim(1,:,1) - uprim(2,:,1)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,0) = ucons(1,:,1) - ucons(2,:,1) + ucons(3,:,1)/3.0
       uprim(1,:,0) = uprim(1,:,1) - uprim(2,:,1) + uprim(3,:,1)/3.0
     end if

  else if (g_lbflag .eq. 1) then
     !--- supersonic inflow
     do imat = 1,nummat
        ucons(1,imat,0) = alpha_fs(imat)
        ucons(1,g_mmi%irmin+imat-1,0) = alpha_fs(imat) * rhomat_fs(imat)
        ucons(1,g_mmi%iemin+imat-1,0) = alpha_fs(imat) * &
          eos3_rhoe(g_gam(imat), g_pc(imat), pr_fs, rhomat_fs(imat), u_fs)
        uprim(1,apr_idx(nummat, imat),0) = alpha_fs(imat) * pr_fs
        uprim(1,mmom_idx(nummat, imat),0) = ucons(1,g_mmi%irmin+imat-1,0) &
          * u_fs
     end do !imat
     uprim(1,vel_idx(nummat, 0),0) = u_fs
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
        uprim(1,apr_idx(nummat, imat),0) = alpha_fs(imat) * pmat
        uprim(1,mmom_idx(nummat, imat),0) = ucons(1,g_mmi%irmin+imat-1,0) &
          * u_fs
     end do !imat
     uprim(1,vel_idx(nummat, 0),0) = u_fs
     ucons(1,g_mmi%imome,0) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,0)) * u_fs

  else if (g_lbflag .eq. 3) then
     !--- periodic
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,0) = ucons(1,:,imax)
       uprim(1,:,0) = uprim(1,:,imax)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,0) = ucons(1,:,imax) + ucons(2,:,imax)
       uprim(1,:,0) = uprim(1,:,imax) + uprim(2,:,imax)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,0) = ucons(1,:,imax) + ucons(2,:,imax) + ucons(3,:,imax)/3.0
       uprim(1,:,0) = uprim(1,:,imax) + uprim(2,:,imax) + uprim(3,:,imax)/3.0
     end if

  else
     write(*,*) "BC-type not set for flag ", g_lbflag

  end if

  !----- right boundary

  if (g_rbflag .eq. -1) then
     !--- exact inlet
     ucons(1,:,imax+1) = gaussian(coord(imax+1),g_time*a_nd)
     uprim(1,apr_idx(nummat,1):apr_idx(nummat,nummat),imax+1) = &
       ucons(1,1:nummat,imax+1) * pr_fs
     uprim(1,mmom_idx(nummat,1):mmom_idx(nummat,nummat),imax+1) = &
       ucons(1,g_mmi%irmin:g_mmi%irmax,imax+1) * u_fs
     uprim(1,vel_idx(nummat, 0),imax+1) = u_fs

  else if (g_rbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,imax+1) = ucons(1,:,imax)
       uprim(1,:,imax+1) = uprim(1,:,imax)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,imax+1) = ucons(1,:,imax) - ucons(2,:,imax)
       uprim(1,:,imax+1) = uprim(1,:,imax) - uprim(2,:,imax)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,imax+1) = ucons(1,:,imax) - ucons(2,:,imax) + ucons(3,:,imax)/3.0
       uprim(1,:,imax+1) = uprim(1,:,imax) - uprim(2,:,imax) + uprim(3,:,imax)/3.0
     end if

  else if (g_rbflag .eq. 1) then
     !--- supersonic inflow
     do imat = 1,nummat
        ucons(1,imat,imax+1) = alpha_fs(imat)
        ucons(1,g_mmi%irmin+imat-1,imax+1) = alpha_fs(imat) * rhomat_fs(imat)
        ucons(1,g_mmi%iemin+imat-1,imax+1) = alpha_fs(imat) * &
          eos3_rhoe(g_gam(imat), g_pc(imat), pr_fs, rhomat_fs(imat), u_fs)
        uprim(1,apr_idx(nummat, imat),imax+1) = alpha_fs(imat) * pr_fs
        uprim(1,mmom_idx(nummat, imat),imax+1) = ucons(1,g_mmi%irmin+imat-1,imax+1) &
          * u_fs
     end do !imat
     uprim(1,vel_idx(nummat, 0),imax+1) = u_fs
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
        uprim(1,apr_idx(nummat, imat),imax+1) = alpha_fs(imat) * pmat
        uprim(1,mmom_idx(nummat, imat),imax+1) = ucons(1,g_mmi%irmin+imat-1,imax+1) &
          * u_fs
     end do !imat
     uprim(1,vel_idx(nummat, 0),imax+1) = u_fs
     ucons(1,g_mmi%imome,imax+1) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax+1)) * u_fs

  else if (g_rbflag .eq. 3) then
     !--- periodic
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,imax+1) = ucons(1,:,1)
       uprim(1,:,imax+1) = uprim(1,:,1)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,imax+1) = ucons(1,:,1) - ucons(2,:,1)
       uprim(1,:,imax+1) = uprim(1,:,1) - uprim(2,:,1)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,imax+1) = ucons(1,:,1) - ucons(2,:,1) + ucons(3,:,1)/3.0
       uprim(1,:,imax+1) = uprim(1,:,1) - uprim(2,:,1) + uprim(3,:,1)/3.0
     end if

  else
     write(*,*) "BC-type not set for flag ", g_rbflag

  end if

  !----- high-order dofs
  if (g_nsdiscr .ge. 1) then
    ucons(2,:,0) = 0.0
    ucons(2,:,imax+1) = 0.0
    uprim(2,:,0) = 0.0
    uprim(2,:,imax+1) = 0.0
      if (g_nsdiscr .ge. 12) then
        ucons(3,:,0) = 0.0
        ucons(3,:,imax+1) = 0.0
        uprim(3,:,0) = 0.0
        uprim(3,:,imax+1) = 0.0
      end if
  end if

end associate

end subroutine get_bc_mm6eq

!-------------------------------------------------------------------------------

subroutine relaxpressure_rdg(ucons, rhsel)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: ig, ie, ieqn, ngauss, imat
data       ngauss/2/
real*8  :: dx, dx2, b3, p_star, rel_time, &
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

    !--- dgp0 or rdgp0p1
    if ((g_nsdiscr .eq. 0) .or. (g_nsdiscr .eq. 1)) then
      do ieqn = 1,g_neqns
        u(ieqn) = ucons(1,ieqn,ie)
      end do !ieqn

    !--- dgp1
    elseif (g_nsdiscr .eq. 11) then
      do ieqn = 1,g_neqns
        u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie)
      end do !ieqn

    !--- rdgp1p2
    elseif (g_nsdiscr .eq. 12) then
      b3 = 0.5*carea(ig)*carea(ig) - 1.0/6.0
      do ieqn = 1,g_neqns
        u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie) + b3 * ucons(3,ieqn,ie)
      end do !ieqn

    end if

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
      aimat = eos3_ss(g_gam(imat), g_pc(imat), rhom(imat), u(imat), pm(imat))
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

      if (g_nsdiscr .ge. 11) then
        rhsel(2,imat,ie) = rhsel(2,imat,ie) + dx2 * carea(ig) * s_alp(imat)
        rhsel(2,g_mmi%iemin+imat-1,ie) = rhsel(2,g_mmi%iemin+imat-1,ie) - dx2 * carea(ig) * p*s_alp(imat)
      end if
    end do !imat

  end do !ig
  end do !ie

end associate

end subroutine relaxpressure_rdg

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
  gaussian(5) = (gaussian(3)+gaussian(4)) * u_fs
  gaussian(6) = al1 * eos3_rhoe(g_gam(1), g_pc(1), pr_fs, rho1, u_fs)
  gaussian(7) = (1.0-al1) * eos3_rhoe(g_gam(2), g_pc(2), pr_fs, rho2, u_fs)

end function

!-------------------------------------------------------------------------------
!----- 2-material contact in volume-fraction
!-------------------------------------------------------------------------------

function contact(x,t)

real*8, intent(in)  :: x, t
real*8  :: xc, al1, rho1, rho2, contact(g_neqns)

  xc  = 0.5 + u_fs*t
  if (x .le. xc) then
    al1 = alpha_fs(1)
  else
    al1 = 1.0-alpha_fs(1)
  end if

  rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
  rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)
  contact(1) = al1
  contact(2) = 1.0-al1
  contact(3) = al1 * rho1
  contact(4) = (1.0-al1) * rho2
  contact(5) = (contact(3)+contact(4)) * u_fs
  contact(6) = al1 * eos3_rhoe(g_gam(1), g_pc(1), pr_fs, rho1, u_fs)
  contact(7) = (1.0-al1) * eos3_rhoe(g_gam(2), g_pc(2), pr_fs, rho2, u_fs)

end function

!-------------------------------------------------------------------------------

function al_star(lplus, lminu, al, ar)
  real*8, intent(in) :: al, ar, lplus, lminu
  real*8 :: al_star

  if (dabs(lplus) .gt. 1.0e-10) then
    al_star = al
  elseif (dabs(lminu) .gt. 1.0e-10) then
    al_star = ar
  else
    al_star = 0.5 * (al+ar)
  end if

end function al_star

!-------------------------------------------------------------------------------

END MODULE rhs_flux_mm6eq
