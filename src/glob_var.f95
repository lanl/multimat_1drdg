MODULE glob_var

implicit none

integer :: imax, neqns, ntstep, &
           nsdiscr, ndof, &
           iprob, i_restart, i_flux, n_opfile, n_screen
real*8  :: dt_u, dt, dt_s

!--- flow properties
real*8  :: u1_fs, rho1_fs, alpha1_fs, alphamin
real*8  :: u2_fs, rho2_fs, pr_fs

!--- isothermal gas
real*8, parameter :: rgas = 287.15
real*8, parameter :: tinf = 300.0
real*8  :: agas = dsqrt(rgas*tinf)

!--- Tait equation constants
real*8  :: rho0 = 1000.0
real*8  :: p0   = 101325.0
real*8  :: k0   = 2.15d9
real*8  :: n0   = 7.15

!--- interface pressure correction coeff
real*8  :: ipsig = 2.0

!--- geometry
real*8,  allocatable :: coord(:)

real*8,  parameter :: Pi = 4.d0*datan(1.d0)

CONTAINS

!-----------------------------------------------------------------------
!----- MC limiter
!-----------------------------------------------------------------------

function mclim(a,b,c)
real*8,intent(in) :: a,b,c
real*8 :: mclim

        if ((a*b>0.0) .and. (a*c>0.0)) then
           mclim = min(2.0*dabs(a), dabs(b), 2.0*dabs(c)) * dsign(1.0,b)
        else
           mclim = 0.0
        end if

end function

!-----------------------------------------------------------------------

END MODULE glob_var
