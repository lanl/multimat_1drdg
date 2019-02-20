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
!----- Pressure from total energy and density using stiffened gas eos
!----------------------------------------------------------------------------------------------

function eos3_pr(gam, p_c, rho, rhoe, u)

real*8, intent(in) :: gam, p_c, rho, rhoe, u

real*8  :: eos3_pr

    eos3_pr = (gam - 1.0) * (rhoe - 0.5*rho*u*u - p_c) - p_c;

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

function eos3_ss(gam, p_c, rho, pr)

real*8, intent(in) :: gam, p_c, rho, pr

real*8  :: pre, eos3_ss

    pre = max(pr+p_c, 1.0d-14);
    if (rho .lt. 1.0d-16) &
      write(*,*) "Warning: negative density encountered in speed of sound calculation: ", rho
    eos3_ss = dsqrt(gam *pre/rho);

end function

!----------------------------------------------------------------------------------------------
!----- Decode primitive variables from conserved 
!----------------------------------------------------------------------------------------------

subroutine decode_uprim(ucons, uprim)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: it, itNewt, ie
real*8  :: pres, pres0, coefa, coefb, taitrel, funcp, dfuncp, rho1
real*8  :: uprim(g_tdof,g_neqns,0:imax+1)

        do ie = 0,imax+1

        itNewt = 100
        pres0 = uprim(1,2,ie)
        pres  = pres0

        !--- Newton loop to find correct pressure
        do it = 1,itNewt

           coefa = ucons(1,1,ie) * rgas * tinf
           coefb = ucons(1,3,ie) / rho0

           taitrel = 1.0 + n0/k0*(pres - p0)

           funcp = (coefa/pres) + (coefb/(taitrel**(1.0/n0))) - 1.0

           ! convergence check
           if (dabs(funcp).lt.1.0e-10) exit
           !print*, ie, ": ", it, pres, funcp

           dfuncp = - ( coefa/(pres*pres) + coefb/k0 * (taitrel**(-(n0+1.0)/n0)) )

           pres = pres - (funcp/dfuncp)

           if (it .eq. itNewt) then
              write(*,*) "Pressure decoding did not converge!"
              write(*,*) " p = ", pres, " p_init = ", pres0
              write(*,*) " Res = ", funcp
              stop
           end if

        end do !it

        if (pres .lt. 0.0) then
           write(*,*) "Negative pressure detected!"
           write(*,*) " cell = ", ie, pres
        end if

        rho1 = eos1_density(pres)

        uprim(1,1,ie) = ucons(1,1,ie)/rho1
        uprim(1,2,ie) = pres
        uprim(1,3,ie) = ucons(1,2,ie)/ucons(1,1,ie)
        uprim(1,4,ie) = ucons(1,4,ie)/ucons(1,3,ie)

        end do !ie

end subroutine decode_uprim

!----------------------------------------------------------------------------------------------
!----- Get primitive variables for a mesh element from conserved
!----------------------------------------------------------------------------------------------

subroutine get_uprim_mm6eq(ucons, uprim)

real*8, intent(in) :: ucons(g_neqns)

integer :: i
real*8  :: almat(g_mmi%nummat), pmat(g_mmi%nummat), tmat(g_mmi%nummat), &
           rhomat(g_mmi%nummat), rhoemat(g_mmi%nummat), u, rho
real*8  :: uprim(g_neqns)

associate (nummat=>g_mmi%nummat)

  ! bulk state
  rho = sum( ucons(g_mmi%irmin:g_mmi%irmax) )
  u = ucons(g_mmi%imome) / rho
  uprim(g_mmi%imome) = u

  ! material states
  do i = 1,nummat
    ! conserved variables
    almat(i) = ucons(i)
    rhomat(i) = ucons(g_mmi%irmin+i-1) / almat(i)
    rhoemat(i) = ucons(g_mmi%iemin+i-1) / almat(i)
    ! primitive variables
    pmat(i) = eos3_pr(g_gam(i), g_pc(i), rhomat(i), rhoemat(i), u)
    tmat(i) = eos3_t(g_gam(i), g_cp(i), g_pc(i), rhomat(i), rhoemat(i), u)

    uprim(i) = almat(i)
    uprim(g_mmi%irmin+i-1) = pmat(i)
    uprim(g_mmi%iemin+i-1) = tmat(i)
  end do !i

end associate

end subroutine get_uprim_mm6eq

!!----- uprimitive calculations with disappearing phase treatment
!subroutine get_uprim_mm6eq(ucons, uprim)
!
!real*8, intent(in) :: ucons(g_neqns)
!
!integer :: i, iamax
!real*8  :: al_eps, rhomax, rhoemax, pmax, tmax, u, rho
!real*8  :: uprim(g_neqns)
!real*8,dimension(g_mmi%nummat) :: almat, pmat, tmat, rhomat, rhoemat
!
!associate (nummat=>g_mmi%nummat)
!
!  al_eps = g_alphamin
!
!  ! bulk state
!  rho = sum( ucons(g_mmi%irmin:g_mmi%irmax) )
!  u = ucons(g_mmi%imome) / rho
!  uprim(g_mmi%imome) = u
!
!  ! majority material
!  almat = ucons(g_mmi%iamin:g_mmi%iamax)
!  iamax = maxloc(almat, 1)
!
!  rhomax  = ucons(g_mmi%irmin+iamax-1)/almat(iamax)
!  rhoemax = ucons(g_mmi%iemin+iamax-1)/almat(iamax)
!  pmax = eos3_pr(g_gam(iamax), g_pc(iamax), rhomax, rhoemax, u)
!  tmax = eos3_t(g_gam(iamax), g_cp(iamax), g_pc(iamax), rhomax, rhoemax, u)
!
!  ! material states
!  do i = 1,nummat
!    if ( (dabs(almat(i)-al_eps) .le. 0.1*al_eps) &
!         .or. (almat(i) .le. al_eps) ) then
!    ! disappearing phase
!      !almat(i) = al_eps
!      pmat(i) = pmax
!      tmat(i) = tmax
!
!    else
!    ! substantial phase
!      rhomat(i) = ucons(g_mmi%irmin+i-1) / almat(i)
!      rhoemat(i) = ucons(g_mmi%iemin+i-1) / almat(i)
!      pmat(i) = eos3_pr(g_gam(i), g_pc(i), rhomat(i), rhoemat(i), u)
!      tmat(i) = eos3_t(g_gam(i), g_cp(i), g_pc(i), rhomat(i), rhoemat(i), u)
!
!    end if
!
!    ! primitive variables
!    uprim(i) = almat(i)
!    uprim(g_mmi%irmin+i-1) = pmat(i)
!    uprim(g_mmi%iemin+i-1) = tmat(i)
!  end do !i
!
!end associate
!
!end subroutine get_uprim_mm6eq

!----------------------------------------------------------------------------------------------

END MODULE eos
