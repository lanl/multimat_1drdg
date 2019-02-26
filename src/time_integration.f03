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

subroutine ExplicitRK3_mm6eq(rhs_mm6eq, reconst_mm6eq, ucons, uconsn)

procedure(), pointer :: rhs_mm6eq, reconst_mm6eq
integer  :: itstep, ielem, idof, ieqn, istage
real*8   :: mm(g_tdof), err_log(g_neqns)
real*8   :: ucons(g_tdof,g_neqns,0:imax+1),uconsn(g_tdof,g_neqns,0:imax+1), &
            uconsi(g_tdof,g_neqns,0:imax+1), &
            ulim(g_tdof,g_neqns,0:imax+1), &
            k1(3),k2(3)
real*8   :: rhsel(g_gdof,g_neqns,imax), cons_err(6)

  g_time = 0.d0

  k1(1) = 0.0     !0.0
  k1(2) = 3.0/4.0 !0.5
  k1(3) = 1.0/3.0

  k2(1) = 1.0     !1.0
  k2(2) = 1.0/4.0 !0.5
  k2(3) = 2.0/3.0

  call meconservation_mm6eq(0,ucons,cons_err)

  do itstep = 1,ntstep

     uconsi = uconsn

     !---------------------------------------------------------
     !--- RK stages
     do istage = 1,3 !2

        rhsel(:,:,:) = 0.d0

        call rhs_mm6eq(uconsi, ulim, rhsel)

        do ielem = 1,imax

        mm(1) = coord(ielem+1)-coord(ielem)
        if (g_nsdiscr .gt. 0) mm(2) = mm(1) / 3.0

        do ieqn = 1,g_neqns
        do idof = 1,g_gdof

          ucons(idof, ieqn,ielem) = &
              k1(istage) *   uconsn(idof,ieqn,ielem) &
            + k2(istage) * ( uconsi(idof,ieqn,ielem) &
                           + dt * rhsel(idof,ieqn,ielem) &
                             / mm(idof) )

        end do !idof
        end do !ieqn
        end do !ielem

        call get_bc_mm6eq(ucons)
        call ignore_tinyphase_mm6eq(ucons)

        uconsi = ucons

     end do !istage
     !---------------------------------------------------------

     g_time = g_time + (dt/a_nd)

     !----- Diagnostics:
     call meconservation_mm6eq(itstep,ucons,cons_err)

     !----- Screen-output:
     if ((itstep.eq.1) .or. (mod(itstep,n_screen).eq.0)) then
     write(*,*) "--------------------------------------------"
     write(*,*) "  itstep: ", itstep, "   Time: ", g_time
     write(*,*) "  Time step: ", dt/a_nd
     write(*,*) "  Conservation: "
     write(*,*) "  Mass:         ", cons_err(3)
     write(*,*) "  Total-energy: ", cons_err(6)
     write(*,*) "--------------------------------------------"
     write(*,*) " "
     end if

     !--- solution update
     uconsn(:,:,:) = ucons(:,:,:)

     !----- File-output:
     if ((mod(itstep,n_opfile).eq.0).or.(itstep.eq.1)) then
     call gnuplot_flow_mm6eq(ulim, itstep)
     call gnuplot_flow_p1_mm6eq(ulim, itstep)
     end if

     if ((mod(itstep,10).eq.0) .or. (itstep.eq.1)) then
     call gnuplot_diagnostics_mm6eq(cons_err, itstep)
     end if

  end do !itstep

  itstep = itstep - 1

  call gnuplot_flow_mm6eq(ulim, itstep)
  call gnuplot_flow_p1_mm6eq(ulim, itstep)
  call gnuplot_diagnostics_mm6eq(cons_err, itstep)

  !----- Update ucons with reconstructed solutions, if any:
  call reconst_mm6eq(ucons)
  call errorcalc_p1(ucons, g_time*a_nd, err_log)

  !----- Screen-output:
  write(*,*) "-----------------------------------------------"
  write(*,*) "-------- FINAL OUTPUT: ----- TVD-RK3: ---------"
  write(*,*) "  itstep: ", itstep, "   Time: ", g_time
  write(*,*) "  Time step: ", dt/a_nd
  write(*,*) "  Conservation: "
  write(*,*) "  Mass:         ", cons_err(3)
  write(*,*) "  Total-energy: ", cons_err(6)
  write(*,*) "  log(||e||): ", err_log(1), 10.0**err_log(1)
  write(*,*) "-----------------------------------------------"
  write(*,*) "-----------------------------------------------"
  write(*,*) " "

end subroutine ExplicitRK3_mm6eq

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
!----- Tiny phase treatment for mm6eq:
!----- ignore it
!----------------------------------------------------------------------------------------------

subroutine ignore_tinyphase_mm6eq(ucons)

integer :: ie, i, iamax
real*8  :: almat(g_mmi%nummat), pmax, tmax, &
           rhomat, rhoemat, &
           al_eps, rho, u, alsum, &
           ucons(g_tdof,g_neqns,0:imax+1)

  al_eps = 0.01*g_alphamin

  do ie = 1,imax

    almat = ucons(1,g_mmi%iamin:g_mmi%iamax,ie)
    iamax = maxloc(almat, 1)

    rhomat  = ucons(1,g_mmi%irmin+iamax-1,ie)/almat(iamax)
    rhoemat = ucons(1,g_mmi%iemin+iamax-1,ie)/almat(iamax)
    rho     = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ie))
    u       = ucons(1,g_mmi%imome,ie)/rho
    pmax = eos3_pr(g_gam(iamax), g_pc(iamax), rhomat, rhoemat, u)
    tmax = eos3_t(g_gam(iamax), g_cp(iamax), g_pc(iamax), rhomat, rhoemat, u)

    alsum = 0.0
    do i = 1,g_mmi%nummat
      rhomat = ucons(1,g_mmi%irmin+i-1,ie)
      !--- phase-i disappearing
      if ( (almat(i) .le. al_eps) .or. (rhomat/rho .lt. 1.d-14) ) then

        ! consistently update derived quantities
        rhomat = eos3_density(g_gam(i), g_cp(i), g_pc(i), pmax, tmax)
        rhoemat = eos3_rhoe(g_gam(i), g_pc(i), pmax, rhomat, u)

        almat(i) = al_eps

        ! update conserved variables
        ucons(1,i,ie) = almat(i)
        ucons(1,g_mmi%irmin+i-1,ie) = almat(i)*rhomat
        ucons(1,g_mmi%iemin+i-1,ie) = almat(i)*rhoemat

        if (g_nsdiscr .ge. 11) then
          ucons(2,i,ie) = 0.0
          ucons(2,g_mmi%irmin+i-1,ie) = 0.0
          ucons(2,g_mmi%iemin+i-1,ie) = 0.0
          if (g_nsdiscr .ge. 12) then
            ucons(3,i,ie) = 0.0
            ucons(3,g_mmi%irmin+i-1,ie) = 0.0
            ucons(3,g_mmi%iemin+i-1,ie) = 0.0
          end if
        end if

      end if

      alsum = alsum+almat(i)
    end do !i
    almat = ucons(1,1:g_mmi%nummat,ie)/alsum

    !if (dabs(alsum-1.0) .gt. al_eps) write(*,*) &
    !  "WARNING: volumefraction sum =", alsum, " element:", ie

  end do !ie

end subroutine ignore_tinyphase_mm6eq

!----------------------------------------------------------------------------------------------
!----- Compute mass and energy conservation for the mm6eq system:
!----------------------------------------------------------------------------------------------

subroutine meconservation_mm6eq(its, ucons, err)

integer, intent(in) :: its
real*8,  intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: ie
real*8  :: mass_1, tenergy_1, massi_1, tenergyi_1, &
           mass_2, tenergy_2, massi_2, tenergyi_2, &
           mass_m, tenergy_m, massi_m, tenergyi_m
real*8  :: err(6)

  !--- initialize
  if (its .eq. 0) then
    g_mass0_1 = 0.0
    g_mass0_2 = 0.0
    g_mass0_m = 0.0
    g_tenergy0_1 = 0.0
    g_tenergy0_2 = 0.0
    g_tenergy0_m = 0.0

  else
    massi_1 = 0.0
    massi_2 = 0.0
    massi_m = 0.0
    tenergyi_1 = 0.0
    tenergyi_2 = 0.0
    tenergyi_m = 0.0

  end if

  do ie = 1,imax

    !--- mass
    mass_1 = ucons(1,g_mmi%irmin,ie)
    mass_2 = ucons(1,g_mmi%irmin+1,ie)
    mass_m = mass_1 + mass_2

    !--- total energy
    tenergy_1 = ucons(1,g_mmi%iemin,ie)
    tenergy_2 = ucons(1,g_mmi%iemin+1,ie)
    tenergy_m = tenergy_1 + tenergy_2

    if (its .eq. 0) then
      g_mass0_1 = g_mass0_1 + mass_1
      g_mass0_2 = g_mass0_2 + mass_2
      g_mass0_m = g_mass0_m + mass_m
      g_tenergy0_1 = g_tenergy0_1 + tenergy_1
      g_tenergy0_2 = g_tenergy0_2 + tenergy_2
      g_tenergy0_m = g_tenergy0_m + tenergy_m

    else
      massi_1 = massi_1 + mass_1
      massi_2 = massi_2 + mass_2
      massi_m = massi_m + mass_m
      tenergyi_1 = tenergyi_1 + tenergy_1
      tenergyi_2 = tenergyi_2 + tenergy_2
      tenergyi_m = tenergyi_m + tenergy_m

    end if

  end do !ie

  if (its .ne. 0) then
    err(1) = (massi_1 - g_mass0_1) / g_mass0_1
    err(2) = (massi_2 - g_mass0_2) / g_mass0_2
    err(3) = (massi_m - g_mass0_m) / g_mass0_m
    err(4) = (tenergyi_1 - g_tenergy0_1) / g_tenergy0_1
    err(5) = (tenergyi_2 - g_tenergy0_2) / g_tenergy0_2
    err(6) = (tenergyi_m - g_tenergy0_m) / g_tenergy0_m
  else
    err(:) = 0.0;
  end if

end subroutine meconservation_mm6eq

!----------------------------------------------------------------------------------------------

END MODULE time_integration
