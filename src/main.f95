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
                       ucons(:,:,:), uconsn(:,:,:)

!----- Read control file:
call read_cntl()

!--- system dimensionality:

write(*,*) " "
if (i_system .eq. 0) then
   ! isothermal single-pressure two-fluid system
   write(*,*) " Isothermal single-pressure two-fluid system: "
   write(*,*) " Pressure equilibrium " 
   write(*,*) " Velocity non-equilibrium "
   g_neqns = 4
else if (i_system .eq. 1) then
   ! pressure non-equilibrium velocity equilibrium two-fluid system
   write(*,*) " Multi-material system: "
   write(*,*) " Pressure non-equilibrium " 
   write(*,*) " Velocity equilibrium "
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
else
   write(*,*) "Error: Incorrect discretization scheme selected:", g_nsdiscr
end if

!----- Allocation:
allocate(uprim(g_tdof,g_neqns,0:imax+1), uprimn(g_tdof,g_neqns,0:imax+1), &
         ucons(g_tdof,g_neqns,0:imax+1), uconsn(g_tdof,g_neqns,0:imax+1))

allocate(coord(0:imax+2))

!----- Mesh generation:
call gen_mesh()

!----- Screen output:
write(*,*) " "
write(*,*) "PREPROCESSING FINISHED."
write(*,*) " nelem = ", imax
write(*,*) " dt =    ", dt_u
write(*,*) " "
write(*,*) " Spatial discretization: "
if (g_nsdiscr .eq. 0) then
  write(*,*) "   P0"
else if (g_nsdiscr .eq. 1) then
  write(*,*) "   P0P1"
else if (g_nsdiscr .eq. 11) then
  write(*,*) "   P1"
end if
write(*,*) " "
write(*,*) " Using appr. Riemann solver: "
if (i_flux .eq. 1) then
  write(*,*) "   Lax-Friedrichs flux."
else if (i_flux .eq. 2) then
  write(*,*) "   AUSM+ flux."
else
  write(*,*) "Invalid flux scheme."
  stop
end if
write(*,*) " "

!----- Initialization:
if (i_system .eq. 0) then
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
if (i_system .eq. 0) then
  call ExplicitRK3_4eq(ucons, uconsn, uprim, uprimn)

else if (i_system .eq. 1) then

  select case(g_nsdiscr)

  case(0)
    call ExplicitRK3_mm6eq(flux_p0_mm6eq, ucons, uconsn)

  case(1)
    call ExplicitRK3_mm6eq(flux_p0p1_mm6eq, ucons, uconsn)

  case(11)
    call ExplicitRK3_mm6eq(flux_p1_mm6eq, ucons, uconsn)

  case default
    write(*,*) "FATAL ERROR: Main3d: Incorrect spatial discretization:", &
               g_nsdiscr
    call exit

  end select

end if

close(33)

!----- Cleanup:
deallocate(uprim, uprimn, &
           ucons, uconsn)

deallocate(coord)

!!---------------------------------------------------------------------------------------

END PROGRAM twofluid_main
