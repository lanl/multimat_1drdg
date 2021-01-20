!!---------------------------------------------------------------------------------------
!!----- Non-equilibrium velocity Isothermal Multiphase and
!!----- non-equilibrium pressure multi-material flow code
!!----- by
!!----- Aditya K Pandare
!!
!!
!!----- Start date : 27 Sept 2017
!!---------------------------------------------------------------------------------------

PROGRAM twofluid_main

!----- Modules:

USE preprocess
USE time_integration

implicit none

!----- Local variable definitions:
integer :: imat
integer, allocatable :: matint_el(:), ndof_el(:,:)
real*8, allocatable :: ucons(:,:,:), uprim(:,:,:), err_log(:)
real*8 :: t_start, t_end, linfty

procedure(reconstruction_p0p1), pointer :: reconst_mm6eq => NULL()
procedure(limiting_p1), pointer :: tvdlimiting_mm6eq => NULL()

!----- Read control file:
call read_cntl()

!--- system dimensionality:

write(*,*) " "
if (i_system .eq. -1) then
   ! k-exactness check
   g_neqns = 3
else if (i_system .eq. 1) then
   ! pressure non-equilibrium velocity equilibrium multi-material system
   g_mmi%iamin = 1
   g_mmi%iamax = g_mmi%nummat
   g_mmi%irmin = g_mmi%iamax+1
   g_mmi%irmax = g_mmi%irmin+g_mmi%nummat-1
   g_mmi%imome = g_mmi%irmax+1
   g_mmi%iemin = g_mmi%imome+1
   g_mmi%iemax = g_mmi%iemin+g_mmi%nummat-1

   g_neqns = g_mmi%iemax
   ! primitive variable vector. See glob_var.f03 for indexing
   g_nprim = 2*g_mmi%nummat+1
else if (i_system .eq. 2) then
   ! hyperelastic solid dynamics
   g_neqns = 8
   g_nprim = 0
end if

if (g_nsdiscr .eq. 0) then
   g_tdof = 1
   g_gdof = 1
else if (g_nsdiscr .eq. 1) then
   g_tdof = 2
   g_gdof = 1
else if (g_nsdiscr .eq. 11) then
   g_tdof = 2
   g_gdof = 2
else if (g_nsdiscr .eq. 12) then
   g_tdof = 3
   g_gdof = 2
else
   write(*,*) "Error: Incorrect discretization scheme selected:", g_nsdiscr
end if

!----- Interface reconstruction option:
if (g_nlim == 6 .or. g_nlim == 7) then
  g_intreco = 1
else
  g_intreco = 0
end if

!----- Allocation:
allocate(ucons(g_tdof,g_neqns,0:imax+1), &
         uprim(g_tdof,g_nprim,0:imax+1), &
         matint_el(0:imax+1), &
         ndof_el(2,0:imax+1), &
         err_log(g_neqns))

allocate(coord(0:imax+2))

!----- Mesh generation:
call gen_mesh()

!----- Screen output:
call screen_output()

call cpu_time(t_start)

!----- Diagnostics file:
open(33,file='diag.dat',status='unknown')

write(33,'(3A12)',advance='no') "# time,", & !1
  "massmix," , &                              !2
  "tenergymix"                                !3
do imat = 1,g_mmi%nummat
  write(33,'(A12)',advance='no') "pkavg,"
end do !imat
write(33,*) "pbavg"

!--- Explicit TVD-RK3:
if (i_system .eq. -1) then

  select case(g_nsdiscr)

  case(1)
    reconst_mm6eq => reconstruction_p0p1
  case(12)
    reconst_mm6eq => reconstruction_p1p2
  case default
    write(*,*) "FATAL ERROR: Main3d: Incorrect spatial discretization:", &
               g_nsdiscr
    call exit

  end select

  ! Initialization
  call init_soln_kex(ucons, ndof_el)
  dt = dt_u
  call reconst_mm6eq(ucons, uprim, ndof_el)
  call errorcalc_p1(ucons, 0.0, err_log, linfty)
  write(*,*) "  quadratic: log(||e||): ", err_log(1), 10.0**err_log(1)
  write(*,*) "      cubic: log(||e||): ", err_log(2), 10.0**err_log(2)
  write(*,*) "   gaussian: log(||e||): ", err_log(3), 10.0**err_log(3)

else if (i_system .eq. 1) then

  select case(g_nsdiscr)

  case(0)
    reconst_mm6eq => reconstruction_p0
    tvdlimiting_mm6eq => limiting_p0
  case(1)
    reconst_mm6eq => reconstruction_p0p1
    tvdlimiting_mm6eq => limiting_p1
  case(11)
    reconst_mm6eq => reconstruction_p1
    tvdlimiting_mm6eq => limiting_p1
  case(12)
    reconst_mm6eq => reconstruction_p1p2
    tvdlimiting_mm6eq => limiting_p2
  case default
    write(*,*) "FATAL ERROR: Main3d: Incorrect spatial discretization:", &
               g_nsdiscr
    call exit

  end select

  ! Initialization
  call init_soln_mm6eq(reconst_mm6eq, ucons, uprim, matint_el, ndof_el)
  dt = dt_u
  call ExplicitRK3_mm6eq(reconst_mm6eq, ucons, uprim, matint_el, ndof_el)

else if (i_system .eq. 2) then

  select case(g_nsdiscr)

  case(0)
    reconst_mm6eq => reconstruction_p0
  !case(1)
  !  reconst_mm6eq => reconstruction_p0p1
  !case(11)
  !  reconst_mm6eq => reconstruction_p1
  !case(12)
  !  reconst_mm6eq => reconstruction_p1p2
  case default
    write(*,*) "FATAL ERROR: Main3d: Incorrect spatial discretization:", &
               g_nsdiscr
    call exit

  end select

  ! Initialization
  call init_soln_soldyn(reconst_mm6eq, ucons, uprim, matint_el, ndof_el)
  dt = dt_u
  call ExplicitRK3_soldyn(reconst_mm6eq, ucons, uprim, matint_el, ndof_el)

end if

close(33)

call cpu_time(t_end)

write(*,*)
write(*,*) "------------------------------------------------------------------"
write(*,*) "Run-time:", t_end-t_start
write(*,*) "Time per time-step:", (t_end-t_start) / ntstep
write(*,*) "------------------------------------------------------------------"
write(*,*)

!----- Cleanup:
deallocate(ucons, &
           uprim, &
           matint_el, &
           ndof_el, &
           err_log)

deallocate(g_gam, g_cp, g_pc, &
           alpha_fs, rhomat_fs)

deallocate(coord)

!!---------------------------------------------------------------------------------------

END PROGRAM twofluid_main
