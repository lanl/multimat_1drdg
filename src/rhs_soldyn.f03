!!------------------------------------------------------------------------------
!!----- Hyperelastic solid mechanics (Flux and source computation module)
!!----- by
!!----- Aditya K Pandare
!!------------------------------------------------------------------------------

MODULE rhs_soldyn

USE reconstruction_soldyn

implicit none

CONTAINS

!-------------------------------------------------------------------------------
!----- rDG RHS:
!-------------------------------------------------------------------------------

subroutine rhs_rdg_soldyn(ucons, uprim, rhsel, matint_el, ndof_el)

real*8  :: rhsel(g_gdof,g_neqns,imax)

integer, intent(in) :: matint_el(0:imax+1), ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

  call flux_rdg_soldyn(ucons, uprim, ndof_el, rhsel)

end subroutine rhs_rdg_soldyn

!-------------------------------------------------------------------------------
!----- rDG Advective-flux contribution to RHS:
!-------------------------------------------------------------------------------

subroutine flux_rdg_soldyn(ucons, uprim, ndof_el, rhsel)

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

real*8  :: rhsel(g_gdof,g_neqns,imax)

  !--- surface integration
  call surfaceint_soldyn_dg(ucons, uprim, ndof_el, rhsel)
  !!--- volume integration
  !call volumeint_dg(ucons, uprim, ndof_el, riemanngrad, vriemann, rhsel)

end subroutine flux_rdg_soldyn

!-------------------------------------------------------------------------------
!----- DG surface contribution to RHS:
!-------------------------------------------------------------------------------

subroutine surfaceint_soldyn_dg(ucons, uprim, ndof_el, rhsel)

integer :: ifc, iel, ier, ieqn, imat
real*8  :: ul(g_neqns), ur(g_neqns), &
           xc, dx, basis(g_tdof), &
           intflux(g_neqns), rhsel(g_gdof,g_neqns,imax), &
           ac_l, ac_r

integer, intent(in) :: ndof_el(2,0:imax+1)
real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                      uprim(g_tdof,g_nprim,0:imax+1)

  do ifc = 1,imax+1

    intflux = 0.0

    iel = ifc - 1
    ier = ifc

    !--- sound-speeds
    ul = ucons(1,:,iel)
    ur = ucons(1,:,ier)
    call get_soldynsoundspeed(ul, ac_l)
    call get_soldynsoundspeed(ur, ac_r)

    !!--- reconstructed and limited values of conserved variables and

    !!--- left element
    !xc = 0.5*(coord(iel+1)+coord(iel))
    !dx = coord(iel+1)-coord(iel)
    !call get_basisfns(coord(ifc), xc, dx, basis)
    !call linc_reconstruction(ucons(:,:,iel), uprim(:,:,iel), basis, ul, pp_l, dx, xc)

    !!--- right element
    !xc = 0.5*(coord(ier+1)+coord(ier))
    !dx = coord(ier+1)-coord(ier)
    !call get_basisfns(coord(ifc), xc, dx, basis)
    !call linc_reconstruction(ucons(:,:,ier), uprim(:,:,ier), basis, ur, pp_r, dx, xc)

    !--- fluxes

    if (i_flux .eq. 1) then
      call llf_soldyn(ul, ur, ac_l, ac_r, intflux)
    !else if (i_flux .eq. 2) then
    !  call ausmplus_soldyn(ul, ur, ac_l, ac_r, intflux)
    !else if (i_flux .eq. 3) then
    !  call hll_soldyn(ul, ur, ac_l, ac_r, intflux)
    else
       write(*,*) "Invalid flux scheme."
       stop
    endif

    if (iel .gt. 0) then
      do ieqn = 1,g_neqns
        rhsel(1,ieqn,iel) = rhsel(1,ieqn,iel) - intflux(ieqn)
        ! high-order contributions
        if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,iel) >= 2)) &
          rhsel(2,ieqn,iel) = rhsel(2,ieqn,iel) - intflux(ieqn)
      end do !ieqn

    end if

    if (ier .lt. (imax+1)) then
      do ieqn = 1,g_neqns
        rhsel(1,ieqn,ier) = rhsel(1,ieqn,ier) + intflux(ieqn)
        ! high-order contributions
        if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ier) >= 2)) &
          rhsel(2,ieqn,ier) = rhsel(2,ieqn,ier) - intflux(ieqn)
      end do !ieqn

    end if

  end do !ifc

end subroutine surfaceint_soldyn_dg

!!-------------------------------------------------------------------------------
!!----- DG volume contribution to RHS:
!!-------------------------------------------------------------------------------
!
!subroutine volumeint_dg(ucons, uprim, ndof_el, rgrad, vriem, rhsel)
!
!integer :: ig, ie, ieqn, ngauss, imat
!
!real*8  :: dx2, b3, p, hmat, viriem, &
!           xg, xc, dx, basis(g_tdof), &
!           u(g_neqns), pp(g_nprim), &
!           rhob, y(g_mmi%nummat), dapdx, &
!           carea(2), weight(2), &
!           cflux(g_neqns), &
!           nflux(g_gdof,g_neqns), &
!           rhsel(g_gdof,g_neqns,imax)
!
!integer, intent(in) :: ndof_el(2,0:imax+1)
!real*8, intent(in) :: rgrad(g_mmi%nummat+1,imax), vriem(g_mmi%nummat+1,imax+1), &
!                      ucons(g_tdof,g_neqns,0:imax+1), &
!                      uprim(g_tdof,g_nprim,0:imax+1)
!
!associate (nummat=>g_mmi%nummat)
!
!  ngauss = get_numqpoints(g_nsdiscr)
!  call rutope(1, ngauss, carea, weight)
!
!  do ie = 1,imax
!  cflux = 0.0
!  nflux = 0.0
!  do ig = 1,ngauss
!
!    dx2 = weight(ig) ! <-- 2.0/dx * weight(ig)/2.0 * dx
!
!    !--- reconstruct high-order solution
!    xc = 0.5*(coord(ie+1)+coord(ie))
!    dx = coord(ie+1)-coord(ie)
!    xg = carea(ig) * 0.5*dx + xc
!    call get_basisfns(xg, xc, dx, basis)
!    call linc_reconstruction(ucons(:,:,ie), uprim(:,:,ie), basis, u, pp, dx, xc)
!
!    call check_volfrac(ie, u)
!
!    viriem = 0.5* (vriem(nummat+1,ie) + vriem(nummat+1,ie+1)) + carea(ig) * rgrad(nummat+1,ie)/2.0
!
!    p = 0.0
!    dapdx = 0.0
!    rhob = sum(u(g_mmi%irmin:g_mmi%irmax))
!    do imat = 1,nummat
!      !p = p + u(imat)*up(g_mmi%irmin+imat-1)
!      p = p + pp(apr_idx(nummat, imat))
!      dapdx = dapdx + rgrad(imat,ie)
!      y(imat) = u(g_mmi%irmin+imat-1) / rhob
!    end do !imat
!
!    !--- flux terms
!    ! momentum flux
!    !cflux(g_mmi%imome) = pp(vel_idx(nummat, 0)) &
!    !  * sum(pp(mmom_idx(nummat,1):mmom_idx(nummat,nummat))) + p
!    cflux(g_mmi%imome) = pp(vel_idx(nummat, 0)) * u(g_mmi%imome) + p
!    nflux(1,g_mmi%imome) = 0.0
!    if (g_nsdiscr .ge. 11) nflux(2,g_mmi%imome) = 0.0
!    do imat = 1,nummat
!      hmat = u(g_mmi%iemin+imat-1) + pp(apr_idx(nummat, imat))
!      !hmat = u(g_mmi%iemin+imat-1) + u(imat)*up(g_mmi%irmin+imat-1)
!      ! other conservative fluxes
!      cflux(imat) = 0.0
!      cflux(g_mmi%irmin+imat-1) = pp(vel_idx(nummat, 0)) * u(g_mmi%irmin+imat-1)
!      cflux(g_mmi%iemin+imat-1) = pp(vel_idx(nummat, 0)) * hmat
!
!      ! non-conservative fluxes
!      nflux(1,imat) = u(imat) * rgrad(nummat+1,ie)
!      nflux(1,g_mmi%irmin+imat-1) = 0.0
!      nflux(1,g_mmi%iemin+imat-1) = - pp(vel_idx(nummat, 0)) * ( y(imat) * dapdx &
!                                                         - rgrad(imat,ie) )
!      if (g_nsdiscr .ge. 11) then
!        nflux(2,imat) = nflux(1,imat) * carea(ig) + u(imat) * viriem * 2.0 !pp(vel_idx(nummat, 0)) * 2.0
!        nflux(2,g_mmi%irmin+imat-1) = 0.0
!        nflux(2,g_mmi%iemin+imat-1) = nflux(1,g_mmi%iemin+imat-1) * carea(ig)
!      end if
!    end do !imat
!
!    do ieqn = 1,g_neqns
!      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ie) >= 2)) &
!        rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + dx2 * cflux(ieqn)
!    end do !ieqn
!
!    do ieqn = 1,g_neqns
!      rhsel(1,ieqn,ie) = rhsel(1,ieqn,ie) + 0.5 * dx2 * nflux(1,ieqn)
!      if ((g_nsdiscr .ge. 11) .and. (ndof_el(1,ie) >= 2)) &
!        rhsel(2,ieqn,ie) = rhsel(2,ieqn,ie) + 0.5 * dx2 * nflux(2,ieqn)
!    end do !ieqn
!
!  end do !ig
!  end do !ie
!
!end associate
!
!end subroutine volumeint_dg

!-------------------------------------------------------------------------------
!----- solid sound-speed:
!-------------------------------------------------------------------------------

subroutine get_soldynsoundspeed(u, ac)

real*8, intent(in) :: u(g_neqns)

real*8, intent(out) :: ac

integer :: imat
real*8 :: pr, ucc, vcc, wcc, vmagcc

  imat = 1
  ucc = u(2)/u(1)
  vcc = u(3)/u(1)
  wcc = u(4)/u(1)
  vmagcc = dsqrt(dot_product([ucc,vcc,wcc],[ucc,vcc,wcc]))
  pr = eos3_pr(g_gam(imat), g_pc(imat), u(1), u(5), vmagcc)

  ! numerical speed of sound choice:
  ac = eos3_ss(g_gam(imat), g_pc(imat), u(1), 1.0, pr)

end subroutine get_soldynsoundspeed

!-------------------------------------------------------------------------------
!----- 2fluid Lax-Friedrichs flux:
!-------------------------------------------------------------------------------

subroutine llf_soldyn(ul, ur, ac_l, ac_r, flux)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), ac_l, ac_r
real*8, intent(out) :: flux(g_neqns)

integer :: imat
real*8 :: ffunc_l(g_neqns), ffunc_r(g_neqns), sigl(3), sigr(3), &
          rho_l, rhou_l(3), rhoe_l, u_l(3), p_l, vmag_l, &
          rho_r, rhou_r(3), rhoe_r, u_r(3), p_r, vmag_r
real*8 :: ac_12, lambda

  flux(:) = 0.0

  imat = 1

  !--- left-state
  sigl(:) = 0.0
  rho_l = ul(1)
  rhou_l(1) = ul(2)
  rhou_l(2) = ul(3)
  rhou_l(3) = ul(4)
  rhoe_l = ul(5)

  u_l(1) = rhou_l(1)/rho_l
  u_l(2) = rhou_l(2)/rho_l
  u_l(3) = rhou_l(3)/rho_l
  vmag_l = dsqrt(dot_product(u_l,u_l))
  p_l = eos3_pr(g_gam(imat), g_pc(imat), rho_l, ul(5), vmag_l)
  sigl(1) = -p_l

  !--- right-state
  sigr(:) = 0.0
  rho_r = ur(1)
  rhou_r(1) = ur(2)
  rhou_r(2) = ur(3)
  rhou_r(3) = ur(4)
  rhoe_r = ur(5)

  u_r(1) = rhou_r(1)/rho_r
  u_r(2) = rhou_r(2)/rho_r
  u_r(3) = rhou_r(3)/rho_r
  vmag_r = dsqrt(dot_product(u_r,u_r))
  p_r = eos3_pr(g_gam(imat), g_pc(imat), rho_r, ur(5), vmag_r)
  sigr(1) = -p_r

  ! numerical speed of sound:
  ac_12 = 0.5 * (ac_l+ac_r)

  lambda = ac_12 + max(dabs(vmag_l),dabs(vmag_r))

  ! flux functions
  ffunc_l = 0.0
  ffunc_r = 0.0

  ffunc_l(1) = u_l(1) * rho_l
  ffunc_l(2) = u_l(1) * rhou_l(1) - sigl(1)
  ffunc_l(3) = u_l(1) * rhou_l(2) - sigl(2)
  ffunc_l(4) = u_l(1) * rhou_l(3) - sigl(3)
  ffunc_l(5) = u_l(1) * (rhoe_l - sigl(1)) - u_l(2)*sigl(2) - u_l(3)*sigl(3)

  ffunc_r(1) = u_r(1) * rho_r
  ffunc_r(2) = u_r(1) * rhou_r(1) - sigr(1)
  ffunc_r(3) = u_r(1) * rhou_r(2) - sigr(2)
  ffunc_r(4) = u_r(1) * rhou_r(3) - sigr(3)
  ffunc_r(5) = u_r(1) * (rhoe_r - sigr(1)) - u_r(2)*sigr(2) - u_r(3)*sigr(3)

  flux = 0.5 * ( ffunc_l+ffunc_r - lambda*(ur-ul) )

end subroutine llf_soldyn

!-------------------------------------------------------------------------------
!----- Boundary conditions:
!-------------------------------------------------------------------------------

subroutine get_bc_soldyn(ucons)

integer :: imat
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  !----- left boundary

  if (g_lbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,0) = ucons(1,:,1)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,0) = ucons(1,:,1) - ucons(2,:,1)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,0) = ucons(1,:,1) - ucons(2,:,1) !+ ucons(3,:,1)/3.0
     end if

  !else if (g_lbflag .eq. 1) then
  !   !--- supersonic inflow
  !   do imat = 1,nummat
  !      ucons(1,imat,0) = alpha_fs(imat)
  !      ucons(1,g_mmi%irmin+imat-1,0) = alpha_fs(imat) * rhomat_fs(imat)
  !      ucons(1,g_mmi%iemin+imat-1,0) = alpha_fs(imat) * &
  !        eos3_rhoe(g_gam(imat), g_pc(imat), pr_fs, rhomat_fs(imat), u_fs)
  !   end do !imat
  !   ucons(1,g_mmi%imome,0) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,0)) * u_fs

  !else if (g_lbflag .eq. 2) then
  !   !--- subsonic inflow
  !   u_conv = ucons(1,g_mmi%imome,1)/sum(ucons(1,g_mmi%irmin:g_mmi%irmax,1))
  !   do imat = 1,nummat
  !      ucons(1,imat,0) = alpha_fs(imat)
  !      ucons(1,g_mmi%irmin+imat-1,0) = alpha_fs(imat) * rhomat_fs(imat)
  !      pmat = eos3_pr(g_gam(imat), g_pc(imat), &
  !                     ucons(1,g_mmi%irmin+imat-1,1)/ucons(1,imat,1), &
  !                     ucons(1,g_mmi%iemin+imat-1,1)/ucons(1,imat,1), u_conv)
  !      ucons(1,g_mmi%iemin+imat-1,0) = alpha_fs(imat) * &
  !        eos3_rhoe(g_gam(imat), g_pc(imat), pmat, rhomat_fs(imat), u_fs)
  !   end do !imat
  !   ucons(1,g_mmi%imome,0) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,0)) * u_fs

  !else if (g_lbflag .eq. 3) then
  !   !--- periodic
  !   if (g_nsdiscr .eq. 0) then
  !     ucons(1,:,0) = ucons(1,:,imax)
  !   else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
  !     ucons(1,:,0) = ucons(1,:,imax) + ucons(2,:,imax)
  !   else if (g_nsdiscr .eq. 12) then
  !     ucons(1,:,0) = ucons(1,:,imax) + ucons(2,:,imax) !+ ucons(3,:,imax)/3.0
  !   end if

  else
     write(*,*) "BC-type not set for flag ", g_lbflag

  end if

  !----- right boundary

  if (g_rbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     if (g_nsdiscr .eq. 0) then
       ucons(1,:,imax+1) = ucons(1,:,imax)
     else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
       ucons(1,:,imax+1) = ucons(1,:,imax) + ucons(2,:,imax)
     else if (g_nsdiscr .eq. 12) then
       ucons(1,:,imax+1) = ucons(1,:,imax) + ucons(2,:,imax) !+ ucons(3,:,imax)/3.0
     end if

  !else if (g_rbflag .eq. 1) then
  !   !--- supersonic inflow
  !   do imat = 1,nummat
  !      ucons(1,imat,imax+1) = alpha_fs(imat)
  !      ucons(1,g_mmi%irmin+imat-1,imax+1) = alpha_fs(imat) * rhomat_fs(imat)
  !      ucons(1,g_mmi%iemin+imat-1,imax+1) = alpha_fs(imat) * &
  !        eos3_rhoe(g_gam(imat), g_pc(imat), pr_fs, rhomat_fs(imat), u_fs)
  !   end do !imat
  !   ucons(1,g_mmi%imome,imax+1) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax+1)) * u_fs

  !else if (g_rbflag .eq. 2) then
  !   !--- subsonic inflow
  !   u_conv = ucons(1,g_mmi%imome,imax)/sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax))
  !   do imat = 1,nummat
  !      ucons(1,imat,imax+1) = alpha_fs(imat)
  !      ucons(1,g_mmi%irmin+imat-1,imax+1) = alpha_fs(imat) * rhomat_fs(imat)
  !      pmat = eos3_pr(g_gam(imat), g_pc(imat), &
  !                     ucons(1,g_mmi%irmin+imat-1,imax)/ucons(1,imat,imax), &
  !                     ucons(1,g_mmi%iemin+imat-1,imax)/ucons(1,imat,imax), u_conv)
  !      ucons(1,g_mmi%iemin+imat-1,imax+1) = alpha_fs(imat) * &
  !        eos3_rhoe(g_gam(imat), g_pc(imat), pmat, rhomat_fs(imat), u_fs)
  !   end do !imat
  !   ucons(1,g_mmi%imome,imax+1) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,imax+1)) * u_fs

  !else if (g_rbflag .eq. 3) then
  !   !--- periodic
  !   if (g_nsdiscr .eq. 0) then
  !     ucons(1,:,imax+1) = ucons(1,:,1)
  !   else if ((g_nsdiscr .eq. 1) .or. (g_nsdiscr .eq. 11)) then
  !     ucons(1,:,imax+1) = ucons(1,:,1) - ucons(2,:,1)
  !   else if (g_nsdiscr .eq. 12) then
  !     ucons(1,:,imax+1) = ucons(1,:,1) - ucons(2,:,1) !+ ucons(3,:,1)/3.0
  !   end if

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

end subroutine get_bc_soldyn

!-------------------------------------------------------------------------------

END MODULE rhs_soldyn
