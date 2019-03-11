MODULE glob_var

implicit none

integer :: i_system, &
           imax, g_neqns, ntstep, &
           g_nsdiscr, g_gdof, g_tdof, g_nlim, g_nmatint, &
           g_nprelx, &
           g_lbflag, g_rbflag, &
           iprob, i_restart, i_flux, n_opfile, n_screen
real*8  :: dt_u, dt, dt_s, g_time

!--- indexing into the multimaterial system
type mmindex
  integer :: nummat, iamin, iamax, irmin, irmax, imome, iemin, iemax
end type mmindex

type(mmindex) :: g_mmi

!--- flow properties
real*8  :: rho1_fs, alpha1_fs, g_alphamin
real*8  :: rho2_fs
real*8  :: u_fs, pr_fs, t_fs
real*8, allocatable  :: alpha_fs(:), rhomat_fs(:)

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
real*8, allocatable :: g_gam(:), g_cp(:), g_pc(:)

!--- interface pressure correction coeff
real*8  :: g_prelct

!--- conservation
real*8  :: g_mass0_m, g_tenergy0_m

!--- geometry
real*8,  allocatable :: coord(:)

real*8,  parameter :: Pi = 4.d0*datan(1.d0)

CONTAINS

!-----------------------------------------------------------------------
!----- MC limiter
!-----------------------------------------------------------------------

real*8 function mclim(a,b,c)
real*8,intent(in) :: a,b,c

  if ((a*b>0.0) .and. (a*c>0.0)) then
    mclim = min(2.0*dabs(a), dabs(b), 2.0*dabs(c)) * dsign(1.0,b)
  else
    mclim = 0.0
  end if

end function mclim

!-----------------------------------------------------------------------

real*8 function p1basis(x,xc,dx)
real*8, intent(in) :: x, xc, dx

  p1basis = (x-xc)/(0.5*dx)

end function p1basis

!-----------------------------------------------------------------------

real*8 function p2basis(x,xc,dx)
real*8, intent(in) :: x, xc, dx
real*8 :: hdx

  hdx = 0.5*dx
  p2basis = (x-xc)*(x-xc)/(2.0*hdx*hdx) - 1.0/6.0

end function p2basis

!-----------------------------------------------------------------------

END MODULE glob_var
