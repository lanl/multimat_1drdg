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
real*8, allocatable :: uprim(:,:,:), uprimn(:,:,:), &
                       ucons(:,:,:), uconsn(:,:,:), err_log(:)

procedure(rhs_p0p1_mm6eq), pointer :: rhs_mm6eq => NULL()
procedure(reconstruction_p0p1), pointer :: reconst_mm6eq => NULL()

!----- Read control file:
call read_cntl()

!--- system dimensionality:

write(*,*) " "
if (i_system .eq. -1) then
   ! k-exactness check
   g_neqns = 3
else if (i_system .eq. 0) then
   ! isothermal single-pressure two-fluid system
   g_neqns = 4
else if (i_system .eq. 1) then
   ! pressure non-equilibrium velocity equilibrium two-fluid system
   g_neqns = 6
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
allocate(uprim(g_tdof,g_neqns,0:imax+1), uprimn(g_tdof,g_neqns,0:imax+1), &
         ucons(g_tdof,g_neqns,0:imax+1), uconsn(g_tdof,g_neqns,0:imax+1), &
         err_log(g_neqns))

allocate(coord(0:imax+2))

!----- Mesh generation:
call gen_mesh()

!----- Screen output:
call screen_output()

!----- Initialization:
if (i_system .eq. -1) then
   call init_soln_kex(ucons)
else if (i_system .eq. 0) then
   call init_soln_4eq(uprim, uprimn, ucons, uconsn)
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

else if (i_system .eq. 0) then
  call ExplicitRK3_4eq(ucons, uconsn, uprim, uprimn)

else if (i_system .eq. 1) then

  select case(g_nsdiscr)

  case(0)
    rhs_mm6eq => rhs_p0_mm6eq
    reconst_mm6eq => reconstruction_p0
  case(1)
    rhs_mm6eq => rhs_p0p1_mm6eq
    reconst_mm6eq => reconstruction_p0p1
  case(11)
    rhs_mm6eq => rhs_p1_mm6eq
    reconst_mm6eq => reconstruction_p1
  case(12)
    rhs_mm6eq => rhs_p1p2_mm6eq
    reconst_mm6eq => reconstruction_p1p2
  case default
    write(*,*) "FATAL ERROR: Main3d: Incorrect spatial discretization:", &
               g_nsdiscr
    call exit

  end select

  call ExplicitRK3_mm6eq(rhs_mm6eq, reconst_mm6eq, ucons, uconsn)

end if

close(33)

!----- Cleanup:
deallocate(uprim, uprimn, &
           ucons, uconsn, &
           err_log)

deallocate(coord)

!!---------------------------------------------------------------------------------------

END PROGRAM twofluid_main
