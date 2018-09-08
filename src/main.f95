!!---------------------------------------------------------------------------------------
!!----- Preliminary exam for Dr Edwards -----
!!----- Non-equilibrium Isothermal Multiphase
!!----- by
!!----- Aditya K Pandare
!!----- Department of Mechanical and Aerospace Engineering,
!!----- North Carolina State University.
!!
!!
!!----- Start date : 27 Sept 2017
!!---------------------------------------------------------------------------------------

PROGRAM twofluid_main

!----- Modules:

USE glob_var
USE preprocess
USE time_integration

implicit none

!----- Local variable definitions:
real*8, allocatable :: uprim(:,:,:), uprimn(:,:,:), &
                       ucons(:,:), uconsn(:,:)

!----- Read control file:
call read_cntl()

!--- system dimensionality:
neqns = 4
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
call init_soln(uprim, uprimn, ucons, uconsn)

!----- Time-stepping:
dt = dt_u

!--- Backward Euler:
call EulerExplicit(ucons, uconsn, uprim, uprimn)

!----- Cleanup:
deallocate(uprim, uprimn, &
           ucons, uconsn)

deallocate(coord)

!!---------------------------------------------------------------------------------------

END PROGRAM twofluid_main
