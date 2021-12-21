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
!----- Tiny phase treatment for mm6eq:
!----- ignore it
!----------------------------------------------------------------------------------------------

subroutine ignore_tinyphase_mm6eq(ucons, uprim)

integer :: ie, i, mmax
real*8  :: almat(g_mmi%nummat), apmat(g_mmi%nummat), pmax, tmax, &
           rhomat, rhoemat, &
           al_eps, rho, u, alsum, d_al, d_are, are_new, &
           apk, ak, kmat(g_mmi%nummat), &
           p_target, ratio, &
           ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

associate (nummat=>g_mmi%nummat)

  al_eps = 1d-2

  do ie = 1,imax

    !--- find material in largest quantity
    almat = ucons(1,g_mmi%iamin:g_mmi%iamax,ie)
    mmax = maxloc(almat, 1)

    !--- get material pressures
    if (g_pureco == 0) then
      u = ucons(1,g_mmi%imome,ie)/sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ie))
      do i = 1,nummat
        apmat(i) = eos3_alphapr(g_gam(i), g_pc(i), almat(i), &
        ucons(1,g_mmi%irmin+i-1,ie), ucons(1,g_mmi%iemin+i-1,ie), u)
      end do !i
    else
      u = uprim(1,vel_idx(nummat, 0),ie)
      apmat = uprim(1,apr_idx(nummat,1):apr_idx(nummat,nummat),ie)
    end if

    rhomat  = ucons(1,g_mmi%irmin+mmax-1,ie)/almat(mmax)
    rhoemat = ucons(1,g_mmi%iemin+mmax-1,ie)/almat(mmax)
    rho     = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ie))
    pmax = apmat(mmax)/almat(mmax)
    tmax = eos3_t(g_gam(mmax), g_cp(mmax), g_pc(mmax), rhomat, rhoemat, u)

    !--- get equilibrium pressure
    !p_target = 0.0
    !ratio = 0.0
    !do i = 1,nummat
    !  rhomat = ucons(1,g_mmi%irmin+i-1,ie)
    !  apk = uprim(1,apr_idx(nummat, i),ie)
    !  ak = eos3_ss(g_gam(i), g_pc(i), rhomat, almat(i), apk)
    !  kmat(i) = rhomat * ak * ak / almat(i)

    !  p_target = p_target + (almat(i) * apk / kmat(i))
    !  ratio = ratio + (almat(i) * almat(i) / kmat(i))
    !end do !i
    !p_target = p_target / ratio
    !p_target = max(p_target, 1d-14)
    p_target = max(pmax, 1d-14)

    !--- correct minority materials and store volume/energy changes
    d_al = 0.0
    d_are = 0.0
    do i = 1,nummat
      apk = apmat(i)/almat(i)
      ! positive volfrac
      if (almat(i) > 0.0) then
        ! pressure relaxation
        if ((almat(i) <= al_eps) &! .and. dabs(apk-pmax)/pmax > 1e-4) &
          .or. (apk+g_pc(i) < 0.0)) then
          rhomat = ucons(1,g_mmi%irmin+i-1,ie) / almat(i)
          are_new = almat(i) * eos3_rhoe(g_gam(i), g_pc(i), p_target, rhomat, u)
          d_are = d_are + ucons(1,g_mmi%iemin+i-1,ie) - are_new

          ucons(1,g_mmi%iemin+i-1,ie) = are_new
          if (g_pureco == 1) &
            uprim(1,apr_idx(nummat, i),ie) = almat(i) * p_target
        end if
      ! negative volfrac
      else if (almat(i) < 0.0) then
        rhomat = eos3_density(g_gam(i), g_cp(i), g_pc(i), p_target, tmax)
        d_al = d_al + (almat(i) - 1d-14)

        ucons(1,g_mmi%iamin+i-1,ie) = 1d-14
        ucons(1,g_mmi%irmin+i-1,ie) = 1d-14 * rhomat
        ucons(1,g_mmi%iemin+i-1,ie) = 1d-14 * eos3_rhoe(g_gam(i), g_pc(i), &
          p_target, rhomat, u)
        if (g_pureco == 1) uprim(1,apr_idx(nummat, i),ie) = 1d-14 * p_target
      end if
    end do !i

    !--- update state of majority material
    ucons(1,g_mmi%iamin+mmax-1,ie) = ucons(1,g_mmi%iamin+mmax-1,ie) + d_al
    almat(mmax) = ucons(1,g_mmi%iamin+mmax-1,ie)
    ucons(1,g_mmi%iemin+mmax-1,ie) = ucons(1,g_mmi%iemin+mmax-1,ie) + d_are
    if (g_pureco == 1) then
      uprim(1,apr_idx(nummat, mmax),ie) = eos3_alphapr(g_gam(mmax), g_pc(mmax), &
        almat(mmax), ucons(1,g_mmi%irmin+mmax-1,ie), &
        ucons(1,g_mmi%iemin+mmax-1,ie), u)
    end if

    !--- enforce unit sum
    alsum = 0.0
    do i = 1,nummat
      alsum = alsum + ucons(1,g_mmi%iamin+i-1,ie)
    end do !i

    do i = 1,nummat
      ucons(1,g_mmi%iamin+i-1,ie) = ucons(1,g_mmi%iamin+i-1,ie) / alsum
      ucons(1,g_mmi%irmin+i-1,ie) = ucons(1,g_mmi%irmin+i-1,ie) / alsum
      ucons(1,g_mmi%iemin+i-1,ie) = ucons(1,g_mmi%iemin+i-1,ie) / alsum
      if (g_pureco == 1) &
        uprim(1,apr_idx(nummat, i),ie) = uprim(1,apr_idx(nummat, i),ie) / alsum
    end do !i

  end do !ie

end associate

end subroutine ignore_tinyphase_mm6eq

!----------------------------------------------------------------------------------------------

END MODULE eos
