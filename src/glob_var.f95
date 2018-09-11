MODULE glob_var

implicit none

integer :: i_system, &
           imax, g_neqns, ntstep, &
           nsdiscr, ndof, &
           g_lbflag, g_rbflag, &
           iprob, i_restart, i_flux, n_opfile, n_screen
real*8  :: dt_u, dt, dt_s

!--- flow properties
real*8  :: u1_fs, rho1_fs, pr1_fs, t1_fs, alpha1_fs, alphamin
real*8  :: u2_fs, rho2_fs, pr2_fs, t2_fs
real*8  :: u_fs, pr_fs

!--- reference quantities for nondimensionalization
real*8  :: a_nd, t_nd, p_nd, rho_nd

!--- isothermal gas
real*8, parameter :: rgas = 287.15
real*8, parameter :: tinf = 300.0
real*8  :: agas = dsqrt(rgas*tinf)

!--- Tait equation constants
real*8  :: rho0 = 1000.0
real*8  :: p0   = 101325.0
real*8  :: k0   = 2.15d9
real*8  :: n0   = 7.15

!--- stiffened gas constants
real*8  :: g_gam1 = 1.4
real*8  :: g_gam2 = 2.8
real*8  :: g_cp1 = 1004.5
real*8  :: g_cp2 = 4186.0
real*8  :: g_pc1 = 0.0
real*8  :: g_pc2 = 8.5d8

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
