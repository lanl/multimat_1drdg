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

!----------------------------------------------------------------------------------------------
!----- Tiny phase treatment for mm6eq:
!----- ignore it
!----------------------------------------------------------------------------------------------

subroutine ignore_tinyphase_mm6eq(ucons, uprim)

integer :: ie, i, iamax
real*8  :: almat(g_mmi%nummat), pmax, tmax, &
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
    iamax = maxloc(almat, 1)

    rhomat  = ucons(1,g_mmi%irmin+iamax-1,ie)/almat(iamax)
    rhoemat = ucons(1,g_mmi%iemin+iamax-1,ie)/almat(iamax)
    rho     = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ie))
    u       = uprim(1,vel_idx(nummat, 0),ie)
    pmax = uprim(1,apr_idx(nummat, iamax),ie)/almat(iamax)
    tmax = eos3_t(g_gam(iamax), g_cp(iamax), g_pc(iamax), rhomat, rhoemat, u)

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
    p_target = pmax

    !--- correct minority materials and store volume/energy changes
    d_al = 0.0
    d_are = 0.0
    alsum = 0.0
    do i = 1,nummat
      apk = uprim(1,apr_idx(nummat, i),ie)/almat(i)
      ! positive volfrac
      if (almat(i) > 0.0) then
        ! pressure relaxation
        if ((almat(i) <= al_eps) &! .and. dabs(apk-pmax)/pmax > 1e-4) &
          .or. (apk+g_pc(i) < 0.0)) then
          rhomat = ucons(1,g_mmi%irmin+i-1,ie) / almat(i)
          are_new = almat(i) * eos3_rhoe(g_gam(i), g_pc(i), p_target, rhomat, u)
          d_are = d_are + ucons(1,g_mmi%iemin+i-1,ie) - are_new

          ucons(1,g_mmi%iemin+i-1,ie) = are_new
          uprim(1,apr_idx(nummat, i),ie) = almat(i) * p_target
        end if
      ! negative volfrac
      else if (almat(i) < 0.0) then
        rhomat = eos3_density(g_gam(i), g_cp(i), g_pc(i), p_target, tmax)
        d_al = d_al + (almat(i) - 1d-12)

        ucons(1,g_mmi%iamin+i-1,ie) = 1d-12
        ucons(1,g_mmi%irmin+i-1,ie) = 1d-12 * rhomat
        ucons(1,g_mmi%iemin+i-1,ie) = 1d-12 * eos3_rhoe(g_gam(i), g_pc(i), &
          p_target, rhomat, u)
        uprim(1,apr_idx(nummat, i),ie) = 1d-12 * p_target
      end if
    end do !i

    !--- update state of majority material
    ucons(1,g_mmi%iamin+iamax-1,ie) = ucons(1,g_mmi%iamin+iamax-1,ie) + d_al
    almat(iamax) = ucons(1,g_mmi%iamin+iamax-1,ie)
    ucons(1,g_mmi%iemin+iamax-1,ie) = ucons(1,g_mmi%iemin+iamax-1,ie) + d_are
    uprim(1,apr_idx(nummat, iamax),ie) = almat(iamax) * eos3_pr(g_gam(iamax), &
      g_pc(iamax), ucons(1,g_mmi%irmin+iamax-1,ie)/almat(iamax), &
      ucons(1,g_mmi%iemin+iamax-1,ie)/almat(iamax), u)

  end do !ie

end associate

end subroutine ignore_tinyphase_mm6eq

!----------------------------------------------------------------------------------------------
!----- Check if volume fraction satisfies [0,1] bounds
!----------------------------------------------------------------------------------------------

subroutine check_volfrac(ie,uc)

integer, intent(in) :: ie
real*8,  intent(in) :: uc(g_neqns)

integer :: imat, iamax
real*8  :: eps

associate (nummat=>g_mmi%nummat)

  eps = 0.0001*g_alphamin

  !do imat = 1,nummat

  !  if ( (uc(imat) .lt. eps) .or. (uc(imat) .gt. 1.0-eps) ) then
  !    write(*,*) "------------------------------------------------------------------"
  !    write(*,*) "Bound-violating volume fraction: ", uc(imat)
  !    write(*,*) "  in cell ", ie
  !    write(*,*) "  material-id:  ", imat
  !    write(*,*) "  density:      ", uc(g_mmi%irmin+imat-1)
  !    write(*,*) "  total-energy: ", uc(g_mmi%iemin+imat-1)
  !    write(*,*) "------------------------------------------------------------------"
  !    write(*,*) " "
  !  end if

  !end do !imat

end associate

end subroutine check_volfrac

!----------------------------------------------------------------------------------------------

END MODULE eos
