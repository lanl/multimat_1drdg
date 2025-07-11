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

subroutine rhs_rdg_mm6eq(ucons, uprim, rhsel, matint_el, ndof_el)

real*8  :: rhsel(g_gdof,g_neqns,imax)

integer, intent(in) :: matint_el(0:imax+1), ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

  call flux_rdg_mm6eq(ucons, uprim, ndof_el, rhsel)

  if (g_nprelx .eq. 1) then
    call relaxpressure_rdg(ucons, uprim, ndof_el, rhsel)
  end if

  if (iprob == -2) then
    call src_mms_tanh_rdg(ndof_el, rhsel)
  else if (iprob == -3) then
    call src_mms_nleg_rdg(ndof_el, rhsel)
  end if

end subroutine rhs_rdg_mm6eq

!-------------------------------------------------------------------------------
!----- rDG Advective-flux contribution to RHS:
!-------------------------------------------------------------------------------

subroutine flux_rdg_mm6eq(ucons, uprim, ndof_el, rhsel)

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

real*8  :: riemanngrad(g_mmi%nummat+1,imax), &
           vriemann(g_mmi%nummat+1,imax+1), &
           rhsel(g_gdof,g_neqns,imax)

  riemanngrad = 0.0
  vriemann = 0.0

  !--- surface integration
  call surfaceint_dg(ucons, uprim, ndof_el, riemanngrad, vriemann, rhsel)
  !--- volume integration
  call volumeint_dg(ucons, uprim, ndof_el, riemanngrad, vriemann, rhsel)

end subroutine flux_rdg_mm6eq

!-------------------------------------------------------------------------------
!----- DG surface contribution to RHS:
!-------------------------------------------------------------------------------

subroutine surfaceint_dg(ucons, uprim, ndof_el, rgrad, vriem, rhsel)

integer :: ifc, iel, ier, ieqn, imat
real*8  :: ul(g_neqns), ur(g_neqns), &
           pp_l(g_nprim), pp_r(g_nprim), &
           xc, dx, basis(g_tdof), &
           intflux(g_neqns), rhsel(g_gdof,g_neqns,imax), &
           ac_l, ac_r, lmda, pplus, pminu, ap_flux(g_mmi%nummat)

real*8  :: rgrad(g_mmi%nummat+1,imax), vriem(g_mmi%nummat+1,imax+1)

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  do ifc = 1,imax+1

  intflux = 0.0

  iel = ifc - 1
  ier = ifc

  !--- reconstructed and limited values of conserved variables and
  !--- primitive quantities: material pressures (p_k) and bulk fluid velocity.

  !--- left element
  xc = 0.5*(coord(iel+1)+coord(iel))
  dx = coord(iel+1)-coord(iel)
  call get_basisfns(xc, xc, dx, basis)
  call linc_reconstruction(ucons(:,:,iel), uprim(:,:,iel), basis, ul, pp_l, dx, xc)
  call get_multimatsoundspeed(ul, pp_l, ac_l)
  call get_basisfns(coord(ifc), xc, dx, basis)
  call linc_reconstruction(ucons(:,:,iel), uprim(:,:,iel), basis, ul, pp_l, dx, xc)

  !--- right element
  xc = 0.5*(coord(ier+1)+coord(ier))
  dx = coord(ier+1)-coord(ier)
  call get_basisfns(xc, xc, dx, basis)
  call linc_reconstruction(ucons(:,:,ier), uprim(:,:,ier), basis, ur, pp_r, dx, xc)
  call get_multimatsoundspeed(ur, pp_r, ac_r)
  call get_basisfns(coord(ifc), xc, dx, basis)
  call linc_reconstruction(ucons(:,:,ier), uprim(:,:,ier), basis, ur, pp_r, dx, xc)

  !if (dabs(sum(ul(g_mmi%iemin:g_mmi%iemax))-3.0000000000705) > 1e-14 .or. &
  !  dabs(sum(ur(g_mmi%iemin:g_mmi%iemax))-3.0000000000705) > 1e-14) then
  !  print*, "surfint-e: ", ifc, sum(ul(g_mmi%iemin:g_mmi%iemax)), &
  !    sum(ur(g_mmi%iemin:g_mmi%iemax))
  !end if
  !if (dabs(sum(ul(1:nummat))-1.0) > 1e-14 .or. dabs(sum(ur(1:nummat))-1.0) > 1e-14) then
  !  print*, "surfint-a: ", ifc, sum(ul(1:nummat)), sum(ur(1:nummat))
  !end if
  !if (dabs(sum(pp_l(apr_idx(nummat,1):apr_idx(nummat,nummat)))-1.0) > 1e-14 .or. &
  !  dabs(sum(pp_r(apr_idx(nummat,1):apr_idx(nummat,nummat)))-1.0) > 1e-14) then
  !  print*, "surfint-p: ", ifc, sum(pp_l(apr_idx(nummat,1):apr_idx(nummat,nummat))), &
  !  sum(pp_r(apr_idx(nummat,1):apr_idx(nummat,nummat)))
  !end if

  !--- fluxes

  if (i_flux .eq. 1) then
     call llf_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, &
                    intflux, lmda, pplus, pminu)
  else if (i_flux .eq. 2) then
     call ausmplus_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, &
                         intflux, lmda, pplus, pminu, iel, ier)
  else if (i_flux .eq. 3) then
     call hll_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, intflux, &
                    lmda, pplus, pminu)
  else if (i_flux .eq. 4) then
     call hllc_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, intflux, &
                     lmda, ap_flux)
  else
     write(*,*) "Invalid flux scheme."
     stop
  endif

  !--- compute gradients of volume fractions and velocity for the
  !--- non-conservative terms from Riemann reconstructed values
  if (i_flux .eq. 4) then
    do imat = 1,nummat
      vriem(imat,ifc) = ap_flux(imat)
    end do !imat
  else
    do imat = 1,nummat
      vriem(imat,ifc) = pplus*pp_l(apr_idx(nummat, imat)) &
        + pminu*pp_r(apr_idx(nummat, imat))
    end do !imat
  end if
  vriem(nummat+1,ifc) = lmda

  if (iel .gt. 0) then
    do ieqn = 1,g_neqns
      rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,iel) >= 2)) &
        rhsel(2,ieqn,iel) = rhsel(2,ieqn,iel) - intflux(ieqn)
      if ((g_nsdiscr .ge. 22) .and. (ndof_el(1,iel) >= 3)) &
        rhsel(3,ieqn,iel) = rhsel(3,ieqn,iel) - 1.0/3.0*intflux(ieqn)
    end do !ieqn

    do imat = 1,nummat
      rgrad(imat,iel) = rgrad(imat,iel) + vriem(imat,ifc)
    end do !imat
    rgrad(nummat+1,iel) = rgrad(nummat+1,iel) + vriem(nummat+1,ifc)

  end if

  if (ier .lt. (imax+1)) then
    do ieqn = 1,g_neqns
      rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ier) >= 2)) &
        rhsel(2,ieqn,ier) = rhsel(2,ieqn,ier) - intflux(ieqn)
      if ((g_nsdiscr .ge. 22) .and. (ndof_el(1,ier) >= 3)) &
        rhsel(3,ieqn,ier) = rhsel(3,ieqn,ier) + 1.0/3.0*intflux(ieqn)
    end do !ieqn

    do imat = 1,nummat
      rgrad(imat,ier) = rgrad(imat,ier) - vriem(imat,ifc)
    end do !imat
    rgrad(nummat+1,ier) = rgrad(nummat+1,ier) - vriem(nummat+1,ifc)

  end if

  end do !ifc

  !do iel = 1,imax
  !  if (dabs(g_fluxch(1,g_mmi%imome,iel)) > 1e-14) &
  !    print*, "WARNING: large pressure force:", iel, g_fluxch(1,g_mmi%imome,iel)
  !end do

end associate

end subroutine surfaceint_dg

!-------------------------------------------------------------------------------
!----- DG volume contribution to RHS:
!-------------------------------------------------------------------------------

subroutine volumeint_dg(ucons, uprim, ndof_el, rgrad, vriem, rhsel)

integer :: ig, ie, ieqn, ngauss, imat

real*8  :: dx2, db3dx, p, hmat, viriem, &
           xg, xc, dx, basis(g_tdof), &
           u(g_neqns), pp(g_nprim), &
           rhob, y(g_mmi%nummat), dapdx, &
           carea(5), weight(5), &
           cflux(g_neqns), &
           nflux(g_gdof,g_neqns), &
           rhsel(g_gdof,g_neqns,imax)

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8, intent(in) :: rgrad(g_mmi%nummat+1,imax), vriem(g_mmi%nummat+1,imax+1), &
                      ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  ngauss = get_numqpoints(g_nsdiscr)
  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax
  cflux = 0.0
  nflux = 0.0
  do ig = 1,ngauss

    dx2 = weight(ig) ! <-- 2.0/dx * weight(ig)/2.0 * dx

    !--- reconstruct high-order solution
    xc = 0.5*(coord(ie+1)+coord(ie))
    dx = coord(ie+1)-coord(ie)
    xg = carea(ig) * 0.5*dx + xc
    call get_basisfns(xg, xc, dx, basis)
    call linc_reconstruction(ucons(:,:,ie), uprim(:,:,ie), basis, u, pp, dx, xc)
    if (g_nsdiscr .ge. 22) db3dx = (xg-xc)/((dx/2.0)*(dx/2.0))

    !if (dabs(sum(u(g_mmi%iemin:g_mmi%iemax))-3.0000000000705) > 1e-14) then
    !  print*, "voluint-e: ", ie, sum(u(g_mmi%iemin:g_mmi%iemax))
    !end if
    !if (dabs(sum(u(1:nummat))-1.0) > 1e-14) then
    !  print*, "voluint-a: ", ie, sum(u(1:nummat))
    !end if
    !if (dabs(sum(pp(apr_idx(nummat,1):apr_idx(nummat,nummat)))-1.0) > 1e-14) then
    !  print*, "voluint-p: ", ie, sum(pp(apr_idx(nummat,1):apr_idx(nummat,nummat)))
    !end if

    viriem = 0.5* (vriem(nummat+1,ie) + vriem(nummat+1,ie+1)) + carea(ig) * rgrad(nummat+1,ie)/2.0

    p = 0.0
    dapdx = 0.0
    rhob = sum(u(g_mmi%irmin:g_mmi%irmax))
    do imat = 1,nummat
      p = p + pp(apr_idx(nummat, imat))
      dapdx = dapdx + rgrad(imat,ie)
      y(imat) = u(g_mmi%irmin+imat-1) / rhob
    end do !imat
    !if (dabs(dapdx) > 1e-14) print*, "WARNING: non-zero dpdx at: ", ie, dapdx

    !--- flux terms
    ! momentum flux
    cflux(g_mmi%imome) = pp(vel_idx(nummat, 0)) * u(g_mmi%imome) + p
    nflux(1,g_mmi%imome) = 0.0
    if (g_nsdiscr .ge. 11) nflux(2,g_mmi%imome) = 0.0
    if (g_nsdiscr .ge. 22) nflux(3,g_mmi%imome) = 0.0
    do imat = 1,nummat
      hmat = u(g_mmi%iemin+imat-1) + pp(apr_idx(nummat, imat))
      ! other conservative fluxes
      cflux(imat) = 0.0
      cflux(g_mmi%irmin+imat-1) = pp(vel_idx(nummat, 0)) * u(g_mmi%irmin+imat-1)
      cflux(g_mmi%iemin+imat-1) = pp(vel_idx(nummat, 0)) * hmat

      ! non-conservative fluxes
      nflux(1,imat) = u(imat) * rgrad(nummat+1,ie)
      nflux(1,g_mmi%irmin+imat-1) = 0.0
      nflux(1,g_mmi%iemin+imat-1) = - pp(vel_idx(nummat, 0)) * ( y(imat) * dapdx &
                                                         - rgrad(imat,ie) )
      if (g_nsdiscr .ge. 11) then
        nflux(2,imat) = nflux(1,imat) * carea(ig) + u(imat) * viriem * 2.0
        nflux(2,g_mmi%irmin+imat-1) = 0.0
        nflux(2,g_mmi%iemin+imat-1) = nflux(1,g_mmi%iemin+imat-1) * carea(ig)
      end if
      if (g_nsdiscr .ge. 22) then
        nflux(3,imat) = nflux(1,imat) * basis(3) + u(imat) * viriem * db3dx*dx
        nflux(3,g_mmi%irmin+imat-1) = 0.0
        nflux(3,g_mmi%iemin+imat-1) = nflux(1,g_mmi%iemin+imat-1) * basis(3)
      end if
    end do !imat

    do ieqn = 1,g_neqns
      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ie) >= 2)) &
        rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + dx2 * cflux(ieqn)
      if ((g_nsdiscr .ge. 22) .and. (ndof_el(1,ie) >= 3)) &
        rhsel(3,ieqn,ie) = rhsel(3,ieqn,ie) + db3dx*0.5*weight(ig)*dx * cflux(ieqn)
    end do !ieqn

    !do ieqn = 1,g_neqns
    !  g_fluxch(1,ieqn,ie) = g_fluxch(1,ieqn,ie) + 0.5*dx2*nflux(1,ieqn)
    !end do !ieqn

    do ieqn = 1,g_neqns
      rhsel(1,ieqn,ie) = rhsel(1,ieqn,ie) + 0.5 * dx2 * nflux(1,ieqn)
      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ie) >= 2)) &
        rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + 0.5 * dx2 * nflux(2,ieqn)
      if ((g_nsdiscr .ge. 22) .and. (ndof_el(1,ie) >= 3)) &
        rhsel(3,ieqn,ie) = rhsel(3,ieqn,ie) + 0.5 * dx2 * nflux(3,ieqn)
    end do !ieqn

  end do !ig

  !do imat = 1,nummat
  !  if (dabs(g_fluxch(1,g_mmi%iemin+imat-1,ie)) > 1e-14) &
  !    print*, "WARNING: large ncn term:", ie, imat, g_fluxch(1,g_mmi%iemin+imat-1,ie)
  !end do !imat

  end do !ie

end associate

end subroutine volumeint_dg

!-------------------------------------------------------------------------------
!----- multi-material sound-speed:
!-------------------------------------------------------------------------------

subroutine get_multimatsoundspeed(u, pp, ac)

real*8, intent(in) :: u(g_neqns), pp(g_nprim)

real*8, intent(out) :: ac

integer :: imat
real*8 :: al(g_mmi%nummat), rhom(g_mmi%nummat), am(g_mmi%nummat), &
          pr, rhob

associate (nummat=>g_mmi%nummat)

  do imat = 1,nummat
    al(imat) = u(imat)
    rhom(imat) = u(g_mmi%irmin+imat-1)
    pr = pp(apr_idx(nummat, imat))
    am(imat) = eos3_ss(g_gam(imat), g_pc(imat), rhom(imat), al(imat), pr)
  end do !imat

  rhob = sum(u(g_mmi%irmin:g_mmi%irmax))

  ! numerical speed of sound choice:

  ! Kapila
  ac = 0.0
  do imat = 1,nummat
    ac = ac + ( rhom(imat)*am(imat)*am(imat) )
  end do !imat
  ac = dsqrt( ac / rhob )

  ! Wood
  !ac_12 = 0.0
  !do imat = 1,nummat
  !  ac_12 = ac_12 + al_12(imat) / ( rho_12*rhom_12(imat)*am_12(imat)*am_12(imat) )
  !end do !imat
  !ac_12 = dsqrt( 1.0 / ac_12 ) / rho_12

end associate

end subroutine get_multimatsoundspeed

!-------------------------------------------------------------------------------
!----- 2fluid Lax-Friedrichs flux:
!-------------------------------------------------------------------------------

subroutine llf_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, &
                     flux, lambda, pplus, pminu)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), &
                      pp_l(g_nprim), pp_r(g_nprim), ac_l, ac_r

integer :: imat
real*8 :: flux(g_neqns), lplus, lminu, lambda, lmag
real*8 :: ffunc_l(g_neqns), ffunc_r(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,hm_l,pm_l, &
                                   arhom_r,hm_r,pm_r

real*8 :: rhou_l, em_l, u_l, rho_l, pi_l, p_l, &
          rhou_r, em_r, u_r, rho_r, pi_r, p_r
real*8 :: ac_12,pplus,pminu

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
    em_l          = ul(g_mmi%iemin+imat-1)

    pi_l         = pp_l(apr_idx(nummat, imat))
    hm_l(imat)   = em_l + pi_l
    p_l = p_l + pi_l
    pm_l(imat)   = pi_l

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r          = ur(g_mmi%iemin+imat-1)

    pi_r         = pp_r(apr_idx(nummat, imat))
    hm_r(imat)   = em_r + pi_r
    p_r = p_r + pi_r
    pm_r(imat)   = pi_r
  end do !imat

  rhou_l = ul(g_mmi%imome)
  rhou_r = ur(g_mmi%imome)
  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l = pp_l(vel_idx(nummat, 0))
  u_r = pp_r(vel_idx(nummat, 0))

  ! numerical speed of sound:
  ac_12 = 0.5 * (ac_l+ac_r)

  lambda = ac_12 + max(dabs(u_l),dabs(u_r))

  ! flux functions
  ffunc_l = 0.0
  ffunc_r = 0.0
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

  lplus = 0.5 * (u_l + lambda)
  lminu = 0.5 * (u_r - lambda)
  pplus = 0.5
  pminu = 0.5

  lmag = dabs(lplus+lminu) + 1.d-16
  lplus = lplus/lmag
  lminu = lminu/lmag
  lambda = lmag*(lplus+lminu)

end associate

end subroutine llf_mm6eq

!-------------------------------------------------------------------------------
!----- 2fluid AUSM+UP:
!-------------------------------------------------------------------------------

subroutine ausmplus_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, &
  flux, lambda, pplus, pminu, iel, ier)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), &
                      pp_l(g_nprim), pp_r(g_nprim), ac_l, ac_r
integer, intent(in) :: iel, ier

integer :: imat
real*8 :: flux(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,hm_l,pm_l, &
                                   arhom_r,hm_r,pm_r

real*8 :: rhou_l, em_l, u_l, m_l, rho_l, pi_l, p_l, &
          rhou_r, em_r, u_r, m_r, rho_r, pi_r, p_r
real*8 :: rho_12, ac_12, &
          f_a, m_12, p_12, m_p, p_u(g_mmi%nummat)
real*8 :: msplus_l(3),msplus_r(3),msminu_l(3),msminu_r(3)
real*8 :: psplus_l,psplus_r,psminu_l,psminu_r,pplus,pminu
real*8 :: temp

real*8 :: lambda,lambda_plus, lambda_minu, lambda_mag

real*8 :: k_p, k_u
!real*8 :: mbar2, umag_fs, m_0

associate (nummat=>g_mmi%nummat)

  k_p = 1.0;
  k_u = 1.0;

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
    em_l          = ul(g_mmi%iemin+imat-1)

    pi_l         = pp_l(apr_idx(nummat, imat))
    hm_l(imat)   = em_l + pi_l
    p_l = p_l + pi_l
    pm_l(imat)   = pi_l

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r          = ur(g_mmi%iemin+imat-1)

    pi_r         = pp_r(apr_idx(nummat, imat))
    hm_r(imat)   = em_r + pi_r
    p_r = p_r + pi_r
    pm_r(imat)   = pi_r
  end do !imat

  rhou_l = ul(g_mmi%imome)
  rhou_r = ur(g_mmi%imome)
  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l = pp_l(vel_idx(nummat, 0))
  u_r = pp_r(vel_idx(nummat, 0))

  ! average states
  rho_12  = 0.5*(rho_l + rho_r)

  ! numerical speed of sound:
  ac_12 = 0.5 * (ac_l+ac_r)

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

  do imat = 1,nummat
    flux(imat)               = lambda_plus*al_l(imat)    + lambda_minu*al_r(imat)
    flux(g_mmi%irmin+imat-1) = lambda_plus*arhom_l(imat) + lambda_minu*arhom_r(imat)
    flux(g_mmi%iemin+imat-1) = lambda_plus*hm_l(imat)    + lambda_minu*hm_r(imat)
    !g_fluxch(1,g_mmi%iemin+imat-1,iel) = g_fluxch(1,g_mmi%iemin+imat-1,iel) - &
    !  (lambda_plus*pp_l(apr_idx(nummat,imat)) + lambda_minu*pp_r(apr_idx(nummat,imat)))
    !g_fluxch(1,g_mmi%iemin+imat-1,ier) = g_fluxch(1,g_mmi%iemin+imat-1,ier) + &
    !  (lambda_plus*pp_l(apr_idx(nummat,imat)) + lambda_minu*pp_r(apr_idx(nummat,imat)))
  end do !imat
  flux(g_mmi%imome) = lambda_plus*rhou_l + lambda_minu*rhou_r + p_12
  !g_fluxch(1,g_mmi%imome,iel) = g_fluxch(1,g_mmi%imome,iel) - p_12
  !g_fluxch(1,g_mmi%imome,ier) = g_fluxch(1,g_mmi%imome,ier) + p_12

  lambda_mag = dabs(lambda) + 1.d-16

  pplus = lambda_plus/(lambda_mag)
  pminu = dabs(lambda_minu/(lambda_mag))

end associate

end subroutine ausmplus_mm6eq

!------------------------------------------------------------------------------
!----- Numerical flux by HLL:
!------------------------------------------------------------------------------

subroutine hll_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, &
  flux, lambda, pplus, pminu)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), &
                      pp_l(g_nprim), pp_r(g_nprim), ac_l, ac_r

integer :: imat
real*8 :: ffunc_l(g_neqns), ffunc_r(g_neqns), flux(g_neqns)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l,hm_l, &
                                   arhom_r,hm_r

real*8 :: rhou_l, em_l, u_l, rho_l, pi_l, p_l, &
          rhou_r, em_r, u_r, rho_r, pi_r, p_r

real*8 :: sl, sr

real*8 :: lambda_plus, lambda_minu, lambda_mag, lambda, pplus, pminu

associate (nummat=>g_mmi%nummat)

  flux(:) = 0.0

  ! material states and bulk pressure
  p_l = 0.0
  p_r = 0.0
  do imat = 1,nummat
  ! ul
    al_l(imat)    = ul(imat)
    arhom_l(imat) = ul(g_mmi%irmin+imat-1)
    em_l          = ul(g_mmi%iemin+imat-1)

    pi_l         = pp_l(apr_idx(nummat, imat))
    hm_l(imat)   = em_l + pi_l
    p_l = p_l + pi_l

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    em_r          = ur(g_mmi%iemin+imat-1)

    pi_r         = pp_r(apr_idx(nummat, imat))
    hm_r(imat)   = em_r + pi_r
    p_r = p_r + pi_r
  end do !imat

  ! bulk state
  rhou_l  = ul(g_mmi%imome)
  rhou_r  = ur(g_mmi%imome)

  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l = pp_l(vel_idx(nummat, 0))
  u_r = pp_r(vel_idx(nummat, 0))

  ! signal velocities
  sl = min((u_l-ac_l), (u_r-ac_r))
  sr = max((u_l+ac_l), (u_r+ac_r))

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

  ! Flux
  if (sl >= 0.0) then
    flux = ffunc_l
    lambda_plus = u_l
    lambda_minu = 0.0
    pplus = 1.0
    pminu = 0.0
  else if (sr <= 0.0) then
    flux = ffunc_r
    lambda_plus = 0.0
    lambda_minu = u_r
    pplus = 0.0
    pminu = 1.0
  else
    flux = ( sr*ffunc_l - sl*ffunc_r + sl*sr*(ur-ul) )&
      / (sr-sl)
    lambda_plus = (sr*u_l - sr*sl) / (sr-sl)
    lambda_minu = (sr*sl - sl*u_r) / (sr-sl)
    pplus = sr / (sr-sl)
    pminu = - sl / (sr-sl)
  end if

  lambda = lambda_plus + lambda_minu

end associate

end subroutine hll_mm6eq

!------------------------------------------------------------------------------
!----- Numerical flux by HLLC:
!------------------------------------------------------------------------------

subroutine hllc_mm6eq(ul, ur, pp_l, pp_r, ac_l, ac_r, &
  flux, lambda, ap_flux)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), &
                      pp_l(g_nprim), pp_r(g_nprim), ac_l, ac_r

integer :: imat
real*8 :: flux(g_neqns), &
          ustar_l(g_neqns), ustar_r(g_neqns), ap_star(g_mmi%nummat)
real*8, dimension(g_mmi%nummat) :: al_l, al_r, &
                                   arhom_l, &
                                   arhom_r, ap_flux

real*8 :: rhou_l, u_l, rho_l, p_l, &
          rhou_r, u_r, rho_r, p_r, p_star

real*8 :: sl, sr, sm

real*8 :: lambda

associate (nummat=>g_mmi%nummat)

  flux(:) = 0.0

  ! material states and bulk pressure
  p_l = 0.0
  p_r = 0.0
  do imat = 1,nummat
  ! ul
    al_l(imat)    = ul(imat)
    arhom_l(imat) = ul(g_mmi%irmin+imat-1)
    p_l = p_l + pp_l(apr_idx(nummat, imat))

  ! ur
    al_r(imat)    = ur(imat)
    arhom_r(imat) = ur(g_mmi%irmin+imat-1)
    p_r = p_r + pp_r(apr_idx(nummat, imat))
  end do !imat

  ! bulk state
  rhou_l  = ul(g_mmi%imome)
  rhou_r  = ur(g_mmi%imome)

  rho_l = sum(arhom_l)
  rho_r = sum(arhom_r)
  u_l = pp_l(vel_idx(nummat, 0))
  u_r = pp_r(vel_idx(nummat, 0))

  ! signal velocities
  sl = min((u_l-ac_l), (u_r-ac_r))
  sr = max((u_l+ac_l), (u_r+ac_r))
  sm = ( rho_r*u_r*(sr-u_r) - rho_l*u_l*(sl-u_l) + p_l-p_r )&
    /( rho_r*(sr-u_r) - rho_l*(sl-u_l) )

  ! middle-zone (star) variables
  do imat=1,nummat
    ap_star(imat) = pp_l(apr_idx(nummat, imat)) &
      + ul(g_mmi%irmin+imat-1)*(u_l-sl)*(u_l-sm)
    p_star = p_star + ap_star(imat)
  end do !imat

  ustar_l(g_mmi%imome) = ((sl-u_l)*rhou_l + (p_star-p_l))/(sl-sm)
  ustar_r(g_mmi%imome) = ((sr-u_r)*rhou_r + (p_star-p_r))/(sr-sm)

  do imat=1,nummat
    ! Left
    ustar_l(imat) = al_l(imat)
    ustar_l(g_mmi%irmin+imat-1) = (sl-u_l) * arhom_l(imat) / (sl-sm);
    ustar_l(g_mmi%iemin+imat-1) = ((sl-u_l) * ul(g_mmi%iemin+imat-1) &
      - pp_l(apr_idx(nummat, imat))*u_l + ap_star(imat)*sm) / (sl-sm);

    ! Right
    ustar_r(imat) = al_r(imat)
    ustar_r(g_mmi%irmin+imat-1) = (sr-u_r) * arhom_r(imat) / (sr-sm);
    ustar_r(g_mmi%iemin+imat-1) = ((sr-u_r) * ur(g_mmi%iemin+imat-1) &
      - pp_r(apr_idx(nummat, imat))*u_r + ap_star(imat)*sm) / (sr-sm);
  end do !imat

  ! Numerical flux
  if (sl >= 0.0) then
    do imat = 1,nummat
      flux(imat) = u_l * al_l(imat)
      flux(g_mmi%irmin+imat-1) = u_l * arhom_l(imat)
      flux(g_mmi%iemin+imat-1) = &
        u_l * (ul(g_mmi%iemin+imat-1) + pp_l(apr_idx(nummat, imat)))
    end do !imat
    flux(g_mmi%imome) = u_l * rhou_l + p_l

    lambda = u_l
    do imat = 1,nummat
      ap_flux(imat) = pp_l(apr_idx(nummat, imat))
    end do !imat

  else if (sl <= 0.0 .and. sm > 0.0) then
    do imat = 1,nummat
      flux(imat) = sm * ustar_l(imat)
      flux(g_mmi%irmin+imat-1) = sm * ustar_l(g_mmi%irmin+imat-1)
      flux(g_mmi%iemin+imat-1) = &
        sm * (ustar_l(g_mmi%iemin+imat-1) + ap_star(imat))
    end do !imat
    flux(g_mmi%imome) = sm * ustar_l(g_mmi%imome) + p_star

    lambda = sm
    do imat = 1,nummat
      ap_flux(imat) = ap_star(imat)
    end do !imat

  else if (sm <= 0.0 .and. sr >= 0.0) then
    do imat = 1,nummat
      flux(imat) = sm * ustar_r(imat)
      flux(g_mmi%irmin+imat-1) = sm * ustar_r(g_mmi%irmin+imat-1)
      flux(g_mmi%iemin+imat-1) = &
        sm * (ustar_r(g_mmi%iemin+imat-1) + ap_star(imat))
    end do !imat
    flux(g_mmi%imome) = sm * ustar_r(g_mmi%imome) + p_star

    lambda = sm
    do imat = 1,nummat
      ap_flux(imat) = ap_star(imat)
    end do !imat

  else !if (sr <= 0.0) then
    do imat = 1,nummat
      flux(imat) = u_r * al_r(imat)
      flux(g_mmi%irmin+imat-1) = u_r * arhom_r(imat)
      flux(g_mmi%iemin+imat-1) = &
        u_r * (ur(g_mmi%iemin+imat-1) + pp_r(apr_idx(nummat, imat)))
    end do !imat
    flux(g_mmi%imome) = u_r * rhou_r + p_r

    lambda = u_r
    do imat = 1,nummat
      ap_flux(imat) = pp_r(apr_idx(nummat, imat))
    end do !imat

  end if

end associate

end subroutine hllc_mm6eq

!-------------------------------------------------------------------------------
!----- Split Mach polynomials for AUSM+UP:
!-------------------------------------------------------------------------------
subroutine splitmach_as(fa, mach, msplus, msminu, psplus, psminu)

real*8, intent(in) :: fa, mach

real*8             :: alph_fa, beta, msplus(3), msminu(3), psplus, psminu

    msplus(1) = 0.5d0*(mach + dabs(mach))
    msminu(1) = 0.5d0*(mach - dabs(mach))

    msplus(2) = +0.25d0*(mach + 1.0d0)*(mach + 1.0d0)
    msminu(2) = -0.25d0*(mach - 1.0d0)*(mach - 1.0d0)

    alph_fa = (3.d0/16.d0) * (-4.d0 + 5.d0*fa*fa)
    beta = 1.d0/8.d0

    if (dabs(mach) .ge. 1.d0) then
    
        msplus(3) = msplus(1)
        msminu(3) = msminu(1)
        psplus    = msplus(1)/mach
        psminu    = msminu(1)/mach
   
    else
    
        msplus(3) = msplus(2)* (1.d0 - 16.d0*beta*msminu(2))
        msminu(3) = msminu(2)* (1.d0 + 16.d0*beta*msplus(2))
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
     if (iprob .eq. -1) then
       ucons(1,:,0) = gaussian(coord(1),g_time*a_nd)
     else if (iprob .eq. -2) then
       ucons(1,:,0) = mms_tanh(coord(1),g_time*a_nd + alpha_dt*dt)
     else if (iprob .eq. -3) then
       ucons(1,:,0) = mms_nleg(coord(1),g_time*a_nd + alpha_dt*dt)
     else if (iprob .eq. -4) then
       ucons(1,:,0) = gaussiantanh(coord(1),g_time*a_nd + alpha_dt*dt)
     else
       write(*,*) "Exact-BC not set for problem", iprob
       stop
     end if

  else if (g_lbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,0) = ucons(1,:,1)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,0) = ucons(1,:,1) - ucons(2,:,1)
     else if ((g_nsdiscr .eq. 12) .or. (g_nsdiscr .eq. 22)) then
       ucons(1,:,0) = ucons(1,:,1) - ucons(2,:,1) + ucons(3,:,1)/3.0
     end if

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
     else if ((g_nsdiscr .eq. 12) .or. (g_nsdiscr .eq. 22)) then
       ucons(1,:,0) = ucons(1,:,imax) + ucons(2,:,imax) + ucons(3,:,imax)/3.0
     end if

  else
     write(*,*) "BC-type not set for flag ", g_lbflag

  end if

  !----- right boundary

  if (g_rbflag .eq. -1) then
     !--- exact inlet
     if (iprob .eq. -1) then
       ucons(1,:,imax+1) = gaussian(coord(imax+1),g_time*a_nd)
     else if (iprob .eq. -2) then
       ucons(1,:,imax+1) = mms_tanh(coord(imax+1),g_time*a_nd + alpha_dt*dt)
     else if (iprob .eq. -3) then
       ucons(1,:,imax+1) = mms_nleg(coord(imax+1),g_time*a_nd + alpha_dt*dt)
     else if (iprob .eq. -4) then
       ucons(1,:,imax+1) = gaussiantanh(coord(imax+1),g_time*a_nd + alpha_dt*dt)
     else
       write(*,*) "Exact-BC not set for problem", iprob
       stop
     end if

  else if (g_rbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,imax+1) = ucons(1,:,imax)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,imax+1) = ucons(1,:,imax) + ucons(2,:,imax)
     else if ((g_nsdiscr .eq. 12) .or. (g_nsdiscr .eq. 22)) then
       ucons(1,:,imax+1) = ucons(1,:,imax) + ucons(2,:,imax) + ucons(3,:,imax)/3.0
     end if

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
     else if ((g_nsdiscr .eq. 12) .or. (g_nsdiscr .eq. 22)) then
       ucons(1,:,imax+1) = ucons(1,:,1) - ucons(2,:,1) + ucons(3,:,1)/3.0
     end if

  else if (g_rbflag .eq. 4) then
     !--- subsonic outflow
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,imax+1) = ucons(1,:,imax)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,imax+1) = ucons(1,:,imax) + ucons(2,:,imax)
     else if ((g_nsdiscr .eq. 12) .or. (g_nsdiscr .eq. 22)) then
       ucons(1,:,imax+1) = ucons(1,:,imax) + ucons(2,:,imax) + ucons(3,:,imax)/3.0
     end if
     u_conv = ucons(1,g_mmi%imome,imax+1)/sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax+1))
     do imat = 1,nummat
        ucons(1,g_mmi%iemin+imat-1,imax+1) = ucons(1,imat,imax+1) * &
          eos3_rhoe(g_gam(imat), g_pc(imat), pr_fs, &
          ucons(1,g_mmi%irmin+imat-1,imax+1)/ucons(1,imat,imax+1), u_conv)
     end do !imat

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

subroutine relaxpressure_rdg(ucons, uprim, ndof_el, rhsel)

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

integer :: ig, ie, ieqn, ngauss, imat
real*8  :: dx, dx2, b3, p_star, rel_time, rhorat, s_lim, &
           xc, xg, basis(g_tdof), &
           rho, p, aimat, nume, deno, &
           carea(5), weight(5)
real*8  :: u(g_neqns), pp(g_nprim), rhsel(g_gdof,g_neqns,imax)

real*8, dimension(g_mmi%nummat) :: rhom, pm, km, s_alp, s_max

associate (nummat=>g_mmi%nummat)

  ngauss = get_numqpoints(g_nsdiscr)
  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax

  dx = coord(ie+1)-coord(ie)

  do ig = 1,ngauss

    dx2 = weight(ig)/2.0 * dx

    !--- reconstruct high-order solution
    xc = 0.5*(coord(ie+1)+coord(ie))
    xg = carea(ig) * 0.5*dx + xc
    call get_basisfns(xg, xc, dx, basis)
    call linc_reconstruction(ucons(:,:,ie), uprim(:,:,ie), basis, u, pp, dx, xc)

    rho  = sum(u(g_mmi%irmin:g_mmi%irmax))
    p = 0.0
    do imat = 1,nummat
      rhom(imat) = u(g_mmi%irmin+imat-1)
      pm(imat)   = pp(apr_idx(nummat, imat))
      p = p + pm(imat)
    end do !imat

    rhorat = 1.0 !0.2 * minval(rhom)/maxval(rhom)

    ! relaxed pressure calculations
    rel_time = 0.0
    nume = 0.0
    deno = 0.0
    do imat = 1,nummat
      aimat = eos3_ss(g_gam(imat), g_pc(imat), rhom(imat), u(imat), pm(imat))
      s_max(imat) = u(imat) /dt !* (aimat+dabs(pp(vel_idx(nummat, 0))))/dx
      km(imat) = rhom(imat) * aimat*aimat
      rel_time = max( rel_time, rhorat * g_prelct * dx/aimat )
      nume = nume + pm(imat)*u(imat)/km(imat)
      deno = deno +  u(imat)*u(imat)/km(imat)
    end do !imat
    p_star = nume/deno

    s_lim = 1.0d14
    do imat = 1,nummat
      s_alp(imat) = 1.0/rel_time * (pm(imat)-p_star*u(imat))*(u(imat)/km(imat))
      s_lim = min(s_lim, s_max(imat)/(dabs(s_alp(imat))+1.0d-12))
    end do !imat

    !if (s_lim .lt. 1.0) then
    !  s_alp = s_lim * s_alp
    !end if

    do imat = 1,nummat
      rhsel(1,imat,ie) = rhsel(1,imat,ie) + dx2 * s_alp(imat)
      rhsel(1,g_mmi%iemin+imat-1,ie) = rhsel(1,g_mmi%iemin+imat-1,ie) - dx2 * p*s_alp(imat)

      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ie) >= 2)) then
        rhsel(2,imat,ie) = rhsel(2,imat,ie) + dx2 * carea(ig) * s_alp(imat)
        rhsel(2,g_mmi%iemin+imat-1,ie) = rhsel(2,g_mmi%iemin+imat-1,ie) - dx2 * carea(ig) * p*s_alp(imat)
      end if
      if ((g_nsdiscr .ge. 22) .and. (ndof_el(1,ie) >= 3)) then
        rhsel(3,imat,ie) = rhsel(3,imat,ie) + dx2 * basis(3) * s_alp(imat)
        rhsel(3,g_mmi%iemin+imat-1,ie) = rhsel(3,g_mmi%iemin+imat-1,ie) - dx2 * basis(3) * p*s_alp(imat)
      end if
    end do !imat

  end do !ig
  end do !ie

end associate

end subroutine relaxpressure_rdg

!----------------------------------------------------------------------------------------------
!----- Tiny phase treatment for mm6eq:
!----- ignore it
!----------------------------------------------------------------------------------------------

subroutine ignore_tinyphase_mm6eq(ucons, uprim)

integer :: ie, i, mmax
real*8  :: dx, xc, basis(g_tdof), &
           almat(g_mmi%nummat), pmax, tmax, uccc(g_neqns), upcc(g_nprim), &
           rhomat, rhoemat, &
           al_eps, rho, u, alsum, d_al, d_are, are_new, &
           apk, ak, kmat(g_mmi%nummat), &
           p_target, ratio, &
           ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  if (g_nmatint == 1 .or. g_intreco > 0) then

  al_eps = 1d-2

  do ie = 1,imax

    !--- get cell averaged unknowns
    dx = coord(ie+1) - coord(ie)
    xc = coord(ie) + 0.5 * dx
    call get_basisfns(xc, xc, dx, basis)
    call linc_reconstruction(ucons(:,:,ie), uprim(:,:,ie), basis, &
      uccc, upcc, dx, xc)

    !--- find material in largest quantity
    almat = uccc(g_mmi%iamin:g_mmi%iamax)
    mmax = maxloc(almat, 1)

    u = upcc(vel_idx(nummat, 0))
    rhomat  = uccc(g_mmi%irmin+mmax-1)/almat(mmax)
    rhoemat = uccc(g_mmi%iemin+mmax-1)/almat(mmax)
    rho     = sum(uccc(g_mmi%irmin:g_mmi%irmax))
    pmax = upcc(apr_idx(nummat,mmax))/almat(mmax)
    tmax = eos3_t(g_gam(mmax), g_cp(mmax), g_pc(mmax), rhomat, rhoemat, u)

    !--- get equilibrium pressure
    !p_target = 0.0
    !ratio = 0.0
    !do i = 1,nummat
    !  rhomat = uccc(g_mmi%irmin+i-1)
    !  apk = upcc(apr_idx(nummat, i))
    !  ak = eos3_ss(g_gam(i), g_pc(i), rhomat, almat(i), apk)
    !  kmat(i) = rhomat * ak * ak / almat(i)

    !  p_target = p_target + (almat(i) * apk / kmat(i))
    !  ratio = ratio + (almat(i) * almat(i) / kmat(i))
    !end do !i
    !p_target = p_target / ratio
    !p_target = max(p_target, 1d-14)
    p_target = max(pmax, 1d-14)

    !--- correct minority materials and store volume/energy changes
    d_al = 0.0
    d_are = 0.0
    do i = 1,nummat
      apk = upcc(apr_idx(nummat,i))/almat(i)
      ! positive volfrac
      if (almat(i) > 0.0) then
        ! pressure relaxation
        if ((almat(i) <= al_eps) &! .and. dabs(apk-pmax)/pmax > 1e-4) &
          .or. (apk+g_pc(i) < 0.0)) then
          rhomat = uccc(g_mmi%irmin+i-1) / almat(i)
          are_new = almat(i) * eos3_rhoe(g_gam(i), g_pc(i), p_target, rhomat, u)
          d_are = d_are + uccc(g_mmi%iemin+i-1) - are_new

          ucons(1,g_mmi%iemin+i-1,ie) = are_new
          if (g_pureco == 1) then
            if (g_pvarreco == 0) then
              uprim(1,apr_idx(nummat, i),ie) = almat(i) * p_target
            else if (g_pvarreco == 1) then
              uprim(1,apr_idx(nummat, i),ie) = p_target
            else if (g_pvarreco == 2) then
              uprim(1,apr_idx(nummat, i),ie) = p_target
            else if (g_pvarreco == 3) then
              uprim(1,apr_idx(nummat, i),ie) = p_target
              uprim(1,rhote_idx(nummat, i),ie) = are_new - uprim(1,rho_idx(nummat,i),ie)
            else if (g_pvarreco == 4) then
              uprim(1,apr_idx(nummat, i),ie) = p_target
              uprim(1,rho_idx(nummat, i),ie) = rhomat
              uprim(1,rhote_idx(nummat, i),ie) = are_new/almat(i)
            end if
          end if
        end if
      ! negative volfrac
      else if (almat(i) < 0.0) then
        rhomat = eos3_density(g_gam(i), g_cp(i), g_pc(i), p_target, tmax)
        d_al = d_al + (almat(i) - 1d-14)

        ucons(1,g_mmi%iamin+i-1,ie) = 1d-14
        ucons(1,g_mmi%irmin+i-1,ie) = 1d-14 * rhomat
        ucons(1,g_mmi%iemin+i-1,ie) = 1d-14 * eos3_rhoe(g_gam(i), g_pc(i), &
          p_target, rhomat, u)
        if (g_pureco == 1) then
          if (g_pvarreco == 0) then
            uprim(1,apr_idx(nummat, i),ie) = 1d-14 * p_target
          else if (g_pvarreco == 1) then
            uprim(1,apr_idx(nummat, i),ie) = p_target
          else if (g_pvarreco == 2) then
            uprim(1,apr_idx(nummat, i),ie) = p_target
          else if (g_pvarreco == 3) then
            uprim(1,apr_idx(nummat, i),ie) = p_target
            uprim(1,rho_idx(nummat, i),ie) = 0.5*ucons(1,g_mmi%irmin+i-1,ie)*u*u
            uprim(1,rhote_idx(nummat, i),ie) = ucons(1,g_mmi%iemin+i-1,ie) - uprim(1,rho_idx(nummat, i),ie)
          else if (g_pvarreco == 4) then
            uprim(1,apr_idx(nummat, i),ie) = p_target
            uprim(1,rho_idx(nummat, i),ie) = rhomat
            uprim(1,rhote_idx(nummat, i),ie) = ucons(1,g_mmi%iemin+i-1,ie)/1d-14
          end if
        end if
      end if
    end do !i

    !--- update state of majority material
    ucons(1,g_mmi%iamin+mmax-1,ie) = ucons(1,g_mmi%iamin+mmax-1,ie) + d_al
    almat(mmax) = ucons(1,g_mmi%iamin+mmax-1,ie)
    ucons(1,g_mmi%iemin+mmax-1,ie) = ucons(1,g_mmi%iemin+mmax-1,ie) + d_are
    if (g_pureco == 1) then
      if (g_pvarreco == 0) then
        uprim(1,apr_idx(nummat, mmax),ie) = eos3_alphapr(g_gam(mmax), g_pc(mmax), &
          almat(mmax), ucons(1,g_mmi%irmin+mmax-1,ie), &
          ucons(1,g_mmi%iemin+mmax-1,ie), u)
      else if (g_pvarreco == 1) then
        uprim(1,apr_idx(nummat, mmax),ie) = eos3_alphapr(g_gam(mmax), g_pc(mmax), &
          almat(mmax), ucons(1,g_mmi%irmin+mmax-1,ie), &
          ucons(1,g_mmi%iemin+mmax-1,ie), u) / almat(mmax)
      else if (g_pvarreco == 2) then
        uprim(1,apr_idx(nummat, mmax),ie) = eos3_alphapr(g_gam(mmax), g_pc(mmax), &
          almat(mmax), ucons(1,g_mmi%irmin+mmax-1,ie), &
          ucons(1,g_mmi%iemin+mmax-1,ie), u) / almat(mmax)
      else if (g_pvarreco == 3) then
        uprim(1,apr_idx(nummat, mmax),ie) = eos3_alphapr(g_gam(mmax), g_pc(mmax), &
          almat(mmax), ucons(1,g_mmi%irmin+mmax-1,ie), &
          ucons(1,g_mmi%iemin+mmax-1,ie), u) / almat(mmax)
        uprim(1,rhote_idx(nummat, mmax),ie) = ucons(1,g_mmi%iemin+mmax-1,ie) - uprim(1,rho_idx(nummat, mmax),ie)
      else if (g_pvarreco == 4) then
        uprim(1,apr_idx(nummat, mmax),ie) = eos3_alphapr(g_gam(mmax), g_pc(mmax), &
          almat(mmax), ucons(1,g_mmi%irmin+mmax-1,ie), &
          ucons(1,g_mmi%iemin+mmax-1,ie), u) / almat(mmax)
        uprim(1,rho_idx(nummat, mmax),ie) = ucons(1,g_mmi%irmin+mmax-1,ie)/almat(mmax)
        uprim(1,rhote_idx(nummat, mmax),ie) = ucons(1,g_mmi%iemin+mmax-1,ie)/almat(mmax)
      end if
    end if

    !--- enforce unit sum
    alsum = 0.0
    do i = 1,nummat
      alsum = alsum + ucons(1,g_mmi%iamin+i-1,ie)
    end do !i

    do i = 1,nummat
      ucons(1,g_mmi%iamin+i-1,ie) = ucons(1,g_mmi%iamin+i-1,ie) / alsum
      ucons(1,g_mmi%irmin+i-1,ie) = ucons(1,g_mmi%irmin+i-1,ie) / alsum
      ucons(1,g_mmi%iemin+i-1,ie) = ucons(1,g_mmi%iemin+i-1,ie) / alsum
      if (g_pureco == 1 .and. g_pvarreco == 0) &
        uprim(1,apr_idx(nummat, i),ie) = uprim(1,apr_idx(nummat, i),ie) / alsum
    end do !i

  end do !ie

  end if

end associate

end subroutine ignore_tinyphase_mm6eq

!-------------------------------------------------------------------------------
!----- Source term integration for MMS advection of equilibrium interface:
!-------------------------------------------------------------------------------

subroutine src_mms_tanh_rdg(ndof_el, rhsel)

integer, intent(in) :: ndof_el(2,0:imax+1)

integer :: ig, ie, ngauss, i
real*8  :: dx, dx2, xc, xg, &
           carea(5), weight(5)
real*8  :: u(g_neqns), rhsel(g_gdof,g_neqns,imax), src(g_neqns)

associate (nummat=>g_mmi%nummat)

  ngauss = get_numqpoints(g_nsdiscr)
  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax

  dx = coord(ie+1)-coord(ie)

  do ig = 1,ngauss

    dx2 = weight(ig)/2.0 * dx

    xc = 0.5*(coord(ie+1)+coord(ie))
    xg = carea(ig) * 0.5*dx + xc
    u = mms_tanh(xg, (g_time*a_nd + alpha_dt*dt))

    !--- source term calculations
    src = 0.0
    do i = 1,nummat
      ! k1 = k2 = 1.0
      src(g_mmi%irmin+i-1) = u_fs*1.0*u(i)
      src(g_mmi%imome) = src(g_mmi%imome) + u_fs*src(g_mmi%irmin+i-1)
      src(g_mmi%iemin+i-1) = 0.5*u_fs*u_fs*src(g_mmi%irmin+i-1)
    end do !i

    !--- contribute to rhs
    do i = 1,g_neqns
      rhsel(1,i,ie) = rhsel(1,i,ie) + dx2 * src(i)

      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ie) >= 2)) then
        rhsel(2,i,ie) = rhsel(2,i,ie) + dx2 * carea(ig) * src(i)
      end if
      if ((g_nsdiscr .ge. 22) .and. (ndof_el(1,ie) >= 3)) then
        rhsel(3,i,ie) = rhsel(3,i,ie) + dx2 * p2basis(xg, xc, dx) * src(i)
      end if
    end do !i

  end do !ig
  end do !ie

end associate

end subroutine src_mms_tanh_rdg

!-------------------------------------------------------------------------------
!----- Source term integration for MMS two-mat nonlinear energy growth:
!-------------------------------------------------------------------------------

subroutine src_mms_nleg_rdg(ndof_el, rhsel)

integer, intent(in) :: ndof_el(2,0:imax+1)

integer :: ig, ie, ngauss, i
real*8  :: dx, dx2, xc, xg, &
           carea(5), weight(5)
real*8  :: t, c_1, c_2, c_3, k, beta, rho0, e0, e_10, h, d_a
real*8  :: u(g_neqns), rhsel(g_gdof,g_neqns,imax), src(g_neqns)

associate (nummat=>g_mmi%nummat)

  ngauss = get_numqpoints(g_nsdiscr)
  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax

  dx = coord(ie+1)-coord(ie)

  do ig = 1,ngauss

    dx2 = weight(ig)/2.0 * dx

    xc = 0.5*(coord(ie+1)+coord(ie))
    xg = carea(ig) * 0.5*dx + xc
    t = g_time*a_nd + alpha_dt*dt
    u = mms_nleg(xg, t)
    call get_nleg_params(c_1, c_2, c_3, k, beta, rho0, e0)

    h = dcos(pi*xg)
    d_a = c_1*dexp(-t)
    e_10 = (-3.0*c_3 - 3.0*k*h*h*t)**(-1.0/3.0)

    !--- source term calculations
    src = 0.0
    src(1) = -d_a
    src(2) = d_a
    src(g_mmi%irmin) = u(g_mmi%irmin)/u(1) * src(1)
    src(g_mmi%irmin+1) = u(g_mmi%irmin+1)/u(2) * src(2)
    src(g_mmi%iemin) = u(g_mmi%iemin)/u(g_mmi%irmin)*src(g_mmi%irmin) &
      + u(g_mmi%irmin)*(k*h*h*(e_10**4.0) + beta*d_a)
    src(g_mmi%iemin+1) = u(g_mmi%iemin+1)/u(g_mmi%irmin+1)*src(g_mmi%irmin+1) &
      - u(g_mmi%irmin+1)*beta*d_a
    src(g_mmi%imome) = - (g_gam(1) - 1.0) * (u(g_mmi%irmin)*(2.0*pi*k*h*t* &
      (e_10**4.0)*dsin(pi*xg)) &
      + u(1)*u(g_mmi%iemin)*(2.0*c_2*xg)/u(g_mmi%irmin) &
      + u(g_mmi%iemin)*(pi*(dcos(pi*xg/2.0)*dsin(pi*xg/2.0))/2.0)/u(1)) &
      + (g_gam(2) - 1) * u(g_mmi%iemin+1)/u(2) &
      * (pi*(dcos(pi*xg/2.0)*dsin(pi*xg/2.0))/2.0)

    !--- contribute to rhs
    do i = 1,g_neqns
      rhsel(1,i,ie) = rhsel(1,i,ie) + dx2 * src(i)

      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ie) >= 2)) then
        rhsel(2,i,ie) = rhsel(2,i,ie) + dx2 * carea(ig) * src(i)
      end if
      if ((g_nsdiscr .ge. 22) .and. (ndof_el(1,ie) >= 3)) then
        rhsel(3,i,ie) = rhsel(3,i,ie) + dx2 * p2basis(xg, xc, dx) * src(i)
      end if
    end do !i

  end do !ig
  end do !ie

end associate

end subroutine src_mms_nleg_rdg

!-------------------------------------------------------------------------------
!----- 2-material Gaussian function in volume-fraction
!-------------------------------------------------------------------------------

function gaussian(x,t)

real*8, intent(in)  :: x, t
real*8  :: xc, al1, rho1, rho2, gaussian(g_neqns)

  xc  = 0.4 + u_fs*t
  al1 = (1.0-2.0*alpha_fs(1)) * dexp( -(x-xc)*(x-xc)/(2.0 * 0.01) ) + alpha_fs(1)

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
!----- 2-material density-Gaussian and tanh-volume fraction
!-------------------------------------------------------------------------------

function gaussiantanh(x,t)

real*8, intent(in)  :: x, t

integer :: i
real*8 :: c1, xc, al_loc(g_mmi%nummat), rhomat, gaussiantanh(g_neqns)

associate (nummat=>g_mmi%nummat)

  xc  = 0.4 + u_fs*t
  rhomat = 10.0*dexp( -(x-xc)*(x-xc)/(2.0 * 0.01) ) + rhomat_fs(1)

  c1 = 10.0
  xc = 0.45 + u_fs*t
  al_loc(1) = (1.0-2.0*alpha_fs(1)) * 0.5 * (1.0 - dtanh(c1*(x-xc))) + alpha_fs(1)
  al_loc(2) = 1.0-al_loc(1)

  ! material states
  do i = 1,nummat
    gaussiantanh(i) = al_loc(i)
    gaussiantanh(g_mmi%irmin+i-1) = gaussiantanh(i) * rhomat
    gaussiantanh(g_mmi%iemin+i-1) = gaussiantanh(i) * &
      eos3_rhoe(g_gam(i), g_pc(i), pr_fs, rhomat, u_fs)
  end do !i

  ! bulk momentum
  gaussiantanh(g_mmi%imome) = &
    sum(gaussiantanh(g_mmi%irmin:g_mmi%irmax)) * u_fs

end associate

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
!----- 2-material Gaussian function in volume-fraction
!-------------------------------------------------------------------------------

function shockentropywave(x, t)

real*8, intent(in) :: x, t

integer :: i
real*8 :: p_loc, u_loc, al_loc(g_mmi%nummat), rhomat, shockentropywave(g_neqns)

associate (nummat=>g_mmi%nummat)

  al_loc(:) = g_alphamin
  if (x < -4.0) then
    p_loc = pr_fs
    u_loc = u_fs
    al_loc(1) = alpha_fs(1)
    al_loc(2) = alpha_fs(2)
  else
    p_loc = 1.0 / p_nd
    u_loc = 0.0
    al_loc(1) = alpha_fs(1)
    al_loc(2) = alpha_fs(2)
  end if

  ! material states
  do i = 1,nummat
    if (x < -4.0) then
      rhomat = eos3_density(g_gam(i), g_cp(i), g_pc(i), pr_fs, t_fs)
    else
      rhomat = (1.0+0.2*dsin(5.0*x))/rho_nd
    end if

    shockentropywave(i) = al_loc(i)
    shockentropywave(g_mmi%irmin+i-1) = shockentropywave(i) * rhomat
    shockentropywave(g_mmi%iemin+i-1) = shockentropywave(i) * &
      eos3_rhoe(g_gam(i), g_pc(i), p_loc, rhomat, u_loc)
  end do !i

  ! bulk momentum
  shockentropywave(g_mmi%imome) = &
    sum(shockentropywave(g_mmi%irmin:g_mmi%irmax)) * u_loc

end associate

!  if (x < -4.0) then
!     rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
!     rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)
!     shockentropywave(1) = alpha_fs(1)
!     shockentropywave(2) = alpha_fs(2)
!     shockentropywave(3) = shockentropywave(1) * rho1
!     shockentropywave(4) = shockentropywave(2) * rho2
!     shockentropywave(5) = (shockentropywave(3)+shockentropywave(4)) * u_fs
!     shockentropywave(6) = shockentropywave(1) * eos3_rhoe(g_gam(1), g_pc(1), pr_fs, rho1, u_fs)
!     shockentropywave(7) = shockentropywave(2) * eos3_rhoe(g_gam(2), g_pc(2), pr_fs, rho2, u_fs)
!  else
!     rho1 = (1.0+0.2*dsin(5.0*x))/rho_nd
!     rho2 = (1.0+0.2*dsin(5.0*x))/rho_nd
!     shockentropywave(1) = alpha_fs(2)
!     shockentropywave(2) = alpha_fs(1)
!     shockentropywave(3) = shockentropywave(1) * rho1
!     shockentropywave(4) = shockentropywave(2) * rho2
!     shockentropywave(5) = 0.0
!     shockentropywave(6) = shockentropywave(1) * eos3_rhoe(g_gam(1), g_pc(1), 1.0/p_nd, rho1, 0.0)
!     shockentropywave(7) = shockentropywave(2) * eos3_rhoe(g_gam(2), g_pc(2), 1.0/p_nd, rho2, 0.0)
!  end if

end function

!-------------------------------------------------------------------------------
!----- MMS 1: advection of equilibrium interface
!-------------------------------------------------------------------------------

function mms_tanh(x, t)

real*8, intent(in) :: x, t

integer :: i
real*8 :: c1, xc, k(g_mmi%nummat), al_loc(g_mmi%nummat), rhomat, mms_tanh(g_neqns)

associate (nummat=>g_mmi%nummat)

  c1 = 10.0
  k(1) = 1.0
  k(2) = 1.0

  xc = 0.45 + u_fs*t
  al_loc(1) = (1.0-2.0*g_alphamin) * 0.5 * (1.0 - dtanh(c1*(x-xc))) + g_alphamin
  al_loc(2) = 1.0-al_loc(1)

  ! material states
  do i = 1,nummat
    rhomat = k(i)*x + rhomat_fs(i)

    mms_tanh(i) = al_loc(i)
    mms_tanh(g_mmi%irmin+i-1) = mms_tanh(i) * rhomat
    mms_tanh(g_mmi%iemin+i-1) = mms_tanh(i) * &
      eos3_rhoe(g_gam(i), g_pc(i), pr_fs, rhomat, u_fs)
  end do !i

  ! bulk momentum
  mms_tanh(g_mmi%imome) = &
    sum(mms_tanh(g_mmi%irmin:g_mmi%irmax)) * u_fs

end associate

end function

!-------------------------------------------------------------------------------
!----- MMS 2: two-material nonlinear energy growth
!-------------------------------------------------------------------------------

function mms_nleg(x, t)

real*8, intent(in) :: x, t

real*8 :: c_1, c_2, c_3, k, beta, d_a, h, rhomat, emat, rho0, e0, &
  mms_nleg(g_neqns)

associate (nummat=>g_mmi%nummat)

  call get_nleg_params(c_1, c_2, c_3, k, beta, rho0, e0)

  h = dcos(pi*x)
  d_a = c_1*dexp(-t)

  ! material volume fractions
  mms_nleg(1) = (1.0-2.0*g_alphamin) * 0.5 * (dcos(pi*x/2.0))**2.0 + g_alphamin &
    + d_a
  mms_nleg(2) = 1.0-mms_nleg(1)

  ! material states
  rhomat = rhomat_fs(1) - c_2*x*x
  emat = (-3.0*c_3 - 3.0*k*h*h*t)**(-1.0/3.0) - beta*d_a
  mms_nleg(g_mmi%irmin) = mms_nleg(1) * rhomat
  mms_nleg(g_mmi%iemin) = mms_nleg(1) * rhomat * emat

  rhomat = rhomat_fs(2)
  emat = 5.0 + beta*d_a
  mms_nleg(g_mmi%irmin+1) = mms_nleg(2) * rhomat
  mms_nleg(g_mmi%iemin+1) = mms_nleg(2) * rhomat * emat

  ! bulk momentum
  mms_nleg(g_mmi%imome) = 0.0

end associate

end function

!-------------------------------------------------------------------------------

subroutine get_nleg_params(c_1, c_2, c_3, k, beta, rho0, e0)
  real*8 :: c_1, c_2, c_3, k, beta, rho0, e0

  c_1 = 0.5
  c_2 = 1.0
  c_3 = -1.0
  k = 0.8
  beta = 1.0
  rho0 = 30.0
  e0 = 5.0
end subroutine get_nleg_params

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
!----- Evaluate degrees of freedom for element based on whether the cell
!----- contains a material interface, and initialize higher degrees of freedom
!----- for cells refined from FV2 to DGP1
!-------------------------------------------------------------------------------

subroutine fill_ndofel(ucons, ndof_el)

integer :: ndof_el(2,0:imax+1), matint(g_mmi%nummat), ie, new_ndof_el

real*8 :: ucons(g_tdof,g_neqns,0:imax+1)

  ! only applicable to DG(P1)/rDG(P1P2)/DG(P2)
  if (g_nsdiscr > 1) then

    do ie = 0, imax+1

      ! use indicator for algebraic reconstruction to determine interface cells
      if (intrecons_cell(ucons(1,g_mmi%iamin:g_mmi%iamax,ie), matint)) then
        new_ndof_el = 1
      else
        new_ndof_el = g_gdof
      end if

      ! if cell is refined from FV2 to DGP1 initialize the high-order DOFs
      if (new_ndof_el > ndof_el(1,ie)) then
        ndof_el(2,ie) = 1
        ucons(2:,:,ie) = 0.0
      else if (new_ndof_el == ndof_el(1,ie)) then
        ndof_el(2,ie) = 0
      else if (new_ndof_el < ndof_el(1,ie)) then
        ndof_el(2,ie) = -1
        !ucons(2:,:,ie) = 0.0
      end if

      ndof_el(1,ie) = new_ndof_el
    end do !ie

  end if

end subroutine fill_ndofel

!-------------------------------------------------------------------------------

END MODULE rhs_flux_mm6eq
