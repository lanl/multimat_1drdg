!!---------------------------------------------------------------------------------------
!!----- Multi-material and multi-phase (EOS module)
!!----- by
!!----- Aditya K Pandare
!!---------------------------------------------------------------------------------------

MODULE eos
      
USE glob_var

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Density from pressure using isothermal 'ideal-gas' eos
!----------------------------------------------------------------------------------------------

function eos1_density(pres)

real*8, intent(in) :: pres

real*8  :: eos1_density

    eos1_density = pres / (rgas*tinf)

end function

!----------------------------------------------------------------------------------------------
!----- Density from pressure using Tait eos
!----------------------------------------------------------------------------------------------

function eos2_density(pres)

real*8, intent(in) :: pres

real*8  :: eos2_density

    eos2_density = rho0 * (1.0 + n0/k0 * (pres - p0))**(1.0/n0)

end function

!----------------------------------------------------------------------------------------------
!----- speed of sound using Tait eos
!----------------------------------------------------------------------------------------------

function eos2_ssound(rho)

real*8, intent(in) :: rho

real*8  :: eos2_ssound

    eos2_ssound = dsqrt(k0/(rho0**n0) * (rho**(n0-1.0)))

end function

!----------------------------------------------------------------------------------------------
!----- Density from pressure and temperature using stiffened gas eos
!----------------------------------------------------------------------------------------------

function eos3_density(gam, cp, p_c, pres, temp)

real*8, intent(in) :: gam, cp, p_c, pres, temp

real*8  :: eos3_density

    eos3_density = gam * (pres + p_c) / ((gam - 1.0) * cp * temp)

end function

!----------------------------------------------------------------------------------------------
!----- Total energy from pressure and density using stiffened gas eos
!----------------------------------------------------------------------------------------------

function eos3_rhoe(gam, p_c, pres, rho, u)

real*8, intent(in) :: gam, p_c, pres, rho, u

real*8  :: eos3_rhoe

    eos3_rhoe = (pres + p_c)/(gam - 1.0) + (0.5*rho*u*u) + p_c;

end function

!----------------------------------------------------------------------------------------------
!----- Total energy from pressure and density using stiffened gas eos
!----------------------------------------------------------------------------------------------

function eos3_arhoe(gam, p_c, apres, arho, u, alp)

real*8, intent(in) :: gam, p_c, apres, arho, u, alp

real*8  :: eos3_arhoe

    eos3_arhoe = (apres + alp*p_c)/(gam - 1.0) + (0.5*arho*u*u) + alp*p_c;

end function

!----------------------------------------------------------------------------------------------
!----- Pressure from total energy and density using stiffened gas eos
!----------------------------------------------------------------------------------------------

function eos3_pr(gam, p_c, rho, rhoe, u)

real*8, intent(in) :: gam, p_c, rho, rhoe, u

real*8  :: eos3_pr

    eos3_pr = (gam - 1.0) * (rhoe - 0.5*rho*u*u - p_c) - p_c;

end function

!----------------------------------------------------------------------------------------------
!----- Vol-frac times pressure from total energy and density using sgeos
!----------------------------------------------------------------------------------------------

function eos3_alphapr(gam, p_c, alpha, arho, arhoe, u)

real*8, intent(in) :: gam, p_c, alpha, arho, arhoe, u

real*8  :: eos3_alphapr

    eos3_alphapr = (gam - 1.0) * (arhoe - 0.5*arho*u*u - alpha*p_c) - alpha*p_c;

end function

!----------------------------------------------------------------------------------------------
!----- Temperature from total energy and density using stiffened gas eos
!----------------------------------------------------------------------------------------------

function eos3_t(gam, cp, p_c, rho, rhoe, u)

real*8, intent(in) :: gam, cp, p_c, rho, rhoe, u

real*8  :: eos3_t

    eos3_t = (rhoe - (0.5*rho*u*u) - p_c) * (gam/(rho * cp));

end function

!----------------------------------------------------------------------------------------------
!----- speed of sound from density and pressure using stiffened gas eos
!----------------------------------------------------------------------------------------------

function eos3_ss(gam, p_c, arho, al, apr)

real*8, intent(in) :: gam, p_c, arho, al, apr

real*8  :: pre, eos3_ss

    pre = max(apr+al*p_c, 1.0d-14);
    if (arho .lt. 1.0d-16) then
      if (al .le. 100.0*g_alphamin) then
        eos3_ss = 0.0
      else
        write(*,*) "Error: zero/negative density encountered in speed of sound calculation: ", arho
        write(*,*) "  Volume-fraction: ", al
        call exit
      end if
    else
      eos3_ss = dsqrt(gam *pre/arho);
    end if

end function

!----------------------------------------------------------------------------------------------
!----- Internal energy from total energy by subtracting KE; not an EOS call
!----------------------------------------------------------------------------------------------

real*8 function intenergy(ret, r, ru)
  real*8, intent(in) :: ret, r, ru

  intenergy = ret - 0.5*ru*ru/r

end function intenergy

!----------------------------------------------------------------------------------------------
!----- Get primitive variables for a mesh element from conserved
!----------------------------------------------------------------------------------------------

subroutine get_uprim_mm6eq(ucons, uprim)

real*8, intent(in) :: ucons(g_neqns)

integer :: i
real*8  :: u, uprim(g_nprim)

associate (nummat=>g_mmi%nummat)

  ! bulk state
  u = ucons(g_mmi%imome) / sum( ucons(g_mmi%irmin:g_mmi%irmax) )
  uprim(vel_idx(nummat,0)) = u

  ! material states
  do i = 1,nummat
    uprim(apr_idx(nummat,i)) = eos3_alphapr(g_gam(i), g_pc(i), ucons(i), &
      ucons(g_mmi%irmin+i-1), ucons(g_mmi%iemin+i-1), u)
  end do !i

end associate

end subroutine get_uprim_mm6eq

!----------------------------------------------------------------------------------------------

END MODULE eos
