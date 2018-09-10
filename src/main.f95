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
                       ucons(:,:), uconsn(:,:)

!----- Read control file:
call read_cntl()

!--- system dimensionality:

write(*,*) " "
if (i_system .eq. 0) then
   ! isothermal single-pressure two-fluid system
   write(*,*) " Isothermal single-pressure two-fluid system: "
   write(*,*) " Pressure equilibrium " 
   write(*,*) " Velocity non-equilibrium "
   neqns = 4
else if (i_system .eq. 1) then
   ! pressure non-equilibrium velocity equilibrium two-fluid system
   write(*,*) " Multi-material system: "
   write(*,*) " Pressure non-equilibrium " 
   write(*,*) " Velocity equilibrium "
   neqns = 6
end if

if (nsdiscr .eq. 0) then
   ndof = 1
else
   ndof = 2
end if

!----- Allocation:
allocate(uprim(ndof,neqns,0:imax+1), uprimn(ndof,neqns,0:imax+1), &
         ucons(neqns,0:imax+1), uconsn(neqns,0:imax+1))

allocate(coord(0:imax+2))

!----- Mesh generation:
call gen_mesh()

write(*,*) " "
write(*,*) "PREPROCESSING FINISHED."
write(*,*) " nelem = ", imax
write(*,*) " dt =    ", dt_u
write(*,*) " "

!----- Initialization:
if (i_system .eq. 0) then
   call init_soln_4eq(uprim, uprimn, ucons, uconsn)
else if (i_system .eq. 1) then
   call init_soln_mm6eq(uprim, uprimn, ucons, uconsn)
end if

!----- Time-stepping:
dt = dt_u

!--- Explicit TVD-RK3:
if (i_system .eq. 0) then
   call ExplicitRK3_4eq(ucons, uconsn, uprim, uprimn)
else if (i_system .eq. 1) then
   call ExplicitRK3_mm6eq(ucons, uconsn, uprim, uprimn)
end if

!----- Cleanup:
deallocate(uprim, uprimn, &
           ucons, uconsn)

deallocate(coord)

!!---------------------------------------------------------------------------------------

END PROGRAM twofluid_main
