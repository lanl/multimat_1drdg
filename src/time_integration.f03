!!---------------------------------------------------------------------------------------
!!----- Multi-material and multi-phase (Time-stepping module)
!!----- by
!!----- Aditya K Pandare
!!---------------------------------------------------------------------------------------

MODULE time_integration

USE preprocess

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Explicit TVD-RK3 time-stepping:
!----------------------------------------------------------------------------------------------

subroutine ExplicitRK3_mm6eq(reconst_mm6eq, ucons, uprim, matint_el, ndof_el)

procedure(), pointer :: reconst_mm6eq
integer  :: itstep, ielem, idof, ieqn, istage
integer  :: matint_el(0:imax+1), ndof_el(2,0:imax+1)
real*8   :: mm(g_tdof), err_log(2,g_neqns), linfty
real*8   :: ucons(g_tdof,g_neqns,0:imax+1),uconsn(g_tdof,g_neqns,0:imax+1), &
            uprim(g_tdof,g_nprim,0:imax+1), &
            k1(3),k2(3),dt_cfl
real*8   :: rhsel(g_gdof,g_neqns,imax), cons_err(2)

  g_time = 0.d0

  k1(1) = 0.0     !0.0
  k1(2) = 3.0/4.0 !0.5
  k1(3) = 1.0/3.0

  k2(1) = 1.0     !1.0
  k2(2) = 1.0/4.0 !0.5
  k2(3) = 2.0/3.0

  call meconservation_mm6eq(0,ucons,cons_err)

  do itstep = 1,ntstep

     !--- solution update
     uconsn(:,:,:) = ucons(:,:,:)

     if (cfl_u > 0) then
       call dtcfl_mm6eq(ucons, uprim, dt_cfl)
       dt = dt_cfl
     end if

     !---------------------------------------------------------
     !--- RK stages
     do istage = 1,3 !2

        if (istage == 1) then
          alpha_dt = 0.0
        else if (istage == 2) then
          alpha_dt = 1.0
        else
          alpha_dt = 0.5
        end if

        rhsel(:,:,:) = 0.d0
        !g_fluxch(:,:,:) = 0.d0

        call rhs_rdg_mm6eq(ucons, uprim, rhsel, matint_el, ndof_el)

        do ielem = 1,imax

        mm(1) = coord(ielem+1)-coord(ielem)
        if (g_nsdiscr .gt. 0) mm(2) = mm(1) / 3.0
        if (g_nsdiscr .gt. 12) mm(3) = mm(1) / 45.0

        do ieqn = 1,g_neqns
        do idof = 1,ndof_el(1,ielem) !g_gdof

          ucons(idof, ieqn,ielem) = &
              k1(istage) *   uconsn(idof,ieqn,ielem) &
            + k2(istage) * ( ucons(idof,ieqn,ielem) &
                           + dt * rhsel(idof,ieqn,ielem) &
                             / mm(idof) )

        end do !idof
        end do !ieqn
        end do !ielem

        call get_bc_mm6eq(ucons)
        call fill_matintel(ucons, matint_el)
        call weak_recons_primitives(ucons, uprim, ndof_el)
        call ignore_tinyphase_mm6eq(ucons, uprim)
        call reconst_mm6eq(ucons, uprim, ndof_el)
        call weak_recons_initconsr(ucons, uprim, ndof_el)

     end do !istage
     !---------------------------------------------------------

     !----- p-refinement for interface reconstruction
     if (g_intreco > 0) call fill_ndofel(ucons, ndof_el)

     g_time = g_time + (dt/a_nd)

     !----- Diagnostics:
     call meconservation_mm6eq(itstep,ucons,cons_err)

     !----- Screen-output:
     if ((itstep.eq.1) .or. (mod(itstep,n_screen).eq.0)) then
     write(*,*) "--------------------------------------------"
     write(*,*) "  itstep: ", itstep, "   Time: ", g_time
     write(*,*) "  Time step: ", dt/a_nd
     write(*,*) "  Conservation: "
     write(*,*) "  Mass:         ", cons_err(1)
     write(*,*) "  Total-energy: ", cons_err(2)
     write(*,*) "--------------------------------------------"
     write(*,*) " "
     end if

     !----- File-output:
     if ((mod(itstep,n_opfile).eq.0).or.(itstep.eq.1)) then
     call gnuplot_flow_mm6eq(ucons, uprim, matint_el, itstep)
     call gnuplot_flow_p1_mm6eq(ucons, uprim, matint_el, itstep)
     end if

     if ((mod(itstep,10).eq.0) .or. (itstep.eq.1)) then
     call gnuplot_diagnostics_mm6eq(ucons, uprim, cons_err, g_time)
     end if

  end do !itstep

  itstep = itstep - 1

  call gnuplot_flow_mm6eq(ucons, uprim, matint_el, itstep)
  call gnuplot_flow_p1_mm6eq(ucons, uprim, matint_el, itstep)
  call gnuplot_diagnostics_mm6eq(ucons, uprim, cons_err, g_time)

  !----- compute L2-error-norm
  call errorcalc_p1(ucons, uprim, g_time*a_nd, err_log, linfty)

  !----- Screen-output:
  write(*,*) "-----------------------------------------------"
  write(*,*) "-------- FINAL OUTPUT: ----- TVD-RK3: ---------"
  write(*,*) "  itstep: ", itstep, "   Time: ", g_time
  write(*,*) "  Time step: ", dt/a_nd
  write(*,*) "  Conservation: "
  write(*,*) "  Mass:         ", cons_err(1)
  write(*,*) "  Total-energy: ", cons_err(2)
  write(*,*) "  log(|e|):   ", err_log(1,1), 10.0**err_log(1,1)
  write(*,*) "  log(||e||): ", err_log(2,1), 10.0**err_log(2,1)
  write(*,*) "  |e|_inf:    ", linfty
  write(*,*) "-----------------------------------------------"
  write(*,*) "-----------------------------------------------"
  write(*,*) " "

end subroutine ExplicitRK3_mm6eq

!----------------------------------------------------------------------------------------------
!----- Determine time step size based on CFL
!----------------------------------------------------------------------------------------------

subroutine dtcfl_mm6eq(ucons, uprim, dt_cfl)

real*8, intent(in)  :: ucons(g_tdof,g_neqns,0:imax+1), &
                       uprim(g_tdof,g_nprim,0:imax+1)
real*8, intent(out) :: dt_cfl

integer :: ie, imat
real*8  :: dx, am, al, rhom, apr

associate (nummat=>g_mmi%nummat)

  dt_cfl = 1e10
  do ie = 1,imax

    dx = coord(ie+1)-coord(ie)

    am = 0.0
    do imat = 1,nummat
      al = ucons(1,imat,ie)
      if (al > 1e-4) then
        rhom = ucons(1,g_mmi%irmin+imat-1,ie)
        apr = uprim(1,apr_idx(nummat,imat),ie)
        am = max(am, eos3_ss(g_gam(imat), g_pc(imat), rhom, al, apr))
      end if
    end do !imat

    am = am + uprim(1,vel_idx(nummat, 0),ie)

    dt_cfl = min(dt_cfl, cfl_u*dx/am);

  end do !ie

end associate

end subroutine dtcfl_mm6eq

!----------------------------------------------------------------------------------------------
!----- Determine what materials are present in cell
!----------------------------------------------------------------------------------------------

subroutine checkmat_mm6eq(ie, ui, mat_present, mat_index)

integer, intent(in) :: ie

integer :: i, mat_present, mat_index(g_mmi%nummat)
real*8  :: rhomat, al_eps, rho, &
           ui(g_tdof,g_neqns)

associate (nummat=>g_mmi%nummat)

  al_eps = 1.d-14

  rho = sum(ui(1,g_mmi%irmin:g_mmi%irmax))

  mat_present = 0
  mat_index = -1

  do i = 1,g_mmi%nummat
    if (ui(1,i) .gt. al_eps) then

      rhomat  = ui(1,g_mmi%irmin+i-1)/ui(1,i)

      mat_present = mat_present + 1
      mat_index(mat_present) = i
    end if
  end do !i

  if (mat_present .eq. 0) then
    write(*,*) "Error: No materials present in cell: ", ie
    call exit
  end if

end associate

end subroutine checkmat_mm6eq

!----------------------------------------------------------------------------------------------
!----- Compute mass and energy conservation for the mm6eq system:
!----------------------------------------------------------------------------------------------

subroutine meconservation_mm6eq(its, ucons, err)

integer, intent(in) :: its
real*8,  intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: ie, imat
real*8  :: mass_m, tenergy_m, massi_m, tenergyi_m
real*8  :: err(2)

associate (nummat=>g_mmi%nummat)

  !--- initialize
  if (its .eq. 0) then
    g_mass0_m = 0.0
    g_tenergy0_m = 0.0

  else
    massi_m = 0.0
    tenergyi_m = 0.0

  end if

  do ie = 1,imax

    mass_m = 0.0
    tenergy_m = 0.0
    do imat = 1,nummat
      !--- mass
      mass_m = mass_m + ucons(1,g_mmi%irmin+imat-1,ie)

      !--- total energy
      tenergy_m = tenergy_m + ucons(1,g_mmi%iemin+imat-1,ie)
    end do !imat

    if (its .eq. 0) then
      g_mass0_m = g_mass0_m + mass_m
      g_tenergy0_m = g_tenergy0_m + tenergy_m

    else
      massi_m = massi_m + mass_m
      tenergyi_m = tenergyi_m + tenergy_m

    end if

  end do !ie

  if (its .ne. 0) then
    err(1) = (massi_m - g_mass0_m) !/ g_mass0_m
    err(2) = (tenergyi_m - g_tenergy0_m) !/ g_tenergy0_m
  else
    err(:) = 0.0;
  end if

end associate

end subroutine meconservation_mm6eq

!----------------------------------------------------------------------------------------------

END MODULE time_integration
