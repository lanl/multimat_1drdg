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
real*8, allocatable :: ucons(:,:,:), uconsn(:,:,:), err_log(:)

procedure(rhs_p0p1_mm6eq), pointer :: rhs_mm6eq => NULL()
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

!----- Allocation:
allocate(ucons(g_tdof,g_neqns,0:imax+1), uconsn(g_tdof,g_neqns,0:imax+1), &
         err_log(g_neqns))

allocate(coord(0:imax+2))

!----- Mesh generation:
call gen_mesh()

!----- Screen output:
call screen_output()

!----- Initialization:
if (i_system .eq. -1) then
   call init_soln_kex(ucons)
else if (i_system .eq. 1) then
   call init_soln_mm6eq(ucons, uconsn)
end if

!----- Time-stepping:
dt = dt_u

!----- Diagnostics file:
open(33,file='diag.dat',status='unknown')

write(33,'(7A12)') "# tstep,", &   !1
                   "mass1,", &     !2
                   "mass2,", &     !3
                   "massmix," , &  !4
                   "tenergy1,", &  !5
                   "tenergy2,", &  !6
                   "tenergymix"    !7

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

  call reconst_mm6eq(ucons)
  call errorcalc_p1(ucons, 0.0, err_log)
  write(*,*) "  quadratic: log(||e||): ", err_log(1), 10.0**err_log(1)
  write(*,*) "      cubic: log(||e||): ", err_log(2), 10.0**err_log(2)
  write(*,*) "   gaussian: log(||e||): ", err_log(3), 10.0**err_log(3)

else if (i_system .eq. 1) then

  select case(g_nsdiscr)

  case(0)
    rhs_mm6eq => rhs_p0_mm6eq
    reconst_mm6eq => reconstruction_p0
    tvdlimiting_mm6eq => limiting_p0
  case(1)
    rhs_mm6eq => rhs_p0p1_mm6eq
    reconst_mm6eq => reconstruction_p0p1
    tvdlimiting_mm6eq => limiting_p1
  case(11)
    rhs_mm6eq => rhs_p1_mm6eq
    reconst_mm6eq => reconstruction_p1
    tvdlimiting_mm6eq => limiting_p1
  case(12)
    rhs_mm6eq => rhs_p1p2_mm6eq
    reconst_mm6eq => reconstruction_p1p2
    tvdlimiting_mm6eq => limiting_p2
  case default
    write(*,*) "FATAL ERROR: Main3d: Incorrect spatial discretization:", &
               g_nsdiscr
    call exit

  end select

  call reconst_mm6eq(ucons)
  call ExplicitRK3_mm6eq(rhs_mm6eq, reconst_mm6eq, ucons, uconsn)

end if

close(33)

!----- Cleanup:
deallocate(ucons, uconsn, &
           err_log)

deallocate(g_gam, g_cp, g_pc, &
           alpha_fs, rhomat_fs)

deallocate(coord)

!!---------------------------------------------------------------------------------------

END PROGRAM twofluid_main
