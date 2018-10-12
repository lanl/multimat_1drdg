!!---------------------------------------------------------------------------------------
!!----- Pressure non-equilibrium multi-material (Flux computation module)
!!----- by
!!----- Aditya K Pandare
!!---------------------------------------------------------------------------------------

MODULE rhs_flux_mm6eq

USE glob_var
USE eos

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Advective-flux contribution to RHS:
!----------------------------------------------------------------------------------------------

subroutine flux_p0_mm6eq(uprim, ucons, rhsel)

integer :: ifc, iel, ier, ieqn
real*8  :: ul(g_neqns), ur(g_neqns), &
           intflux(g_neqns), rhsel(g_neqns,imax), &
           u_conv_l, u_conv_r, uf, pr, p1, p2, &
           alp1f, alp2f, ncnfl, ncnfr, lplus, lminu, lmag

real*8, intent(in) :: uprim(ndof,g_neqns,0:imax+1), ucons(g_neqns,0:imax+1)

  do ifc = 1,imax+1

  iel = ifc - 1
  ier = ifc

  ul(:) = ucons(:,iel)
  ur(:) = ucons(:,ier)

  !--- Conservative fluxes

  if (i_flux .eq. 1) then
     call llf_mm6eq(ul, ur, intflux, lplus, lminu)
  else if (i_flux .eq. 2) then
     call ausmplus_mm6eq(ul, ur, intflux, lplus, lminu, lmag)
  else
     print*, "Invalid flux scheme."
     stop
  endif

  !--- Non-conservative terms
  alp1f = lplus*ucons(1,iel) + lminu*ucons(1,ier)
  alp2f = lplus*(1.0-ucons(1,iel)) + lminu*(1.0-ucons(1,ier))

  ! left element
  u_conv_l = ucons(4,iel)/(ucons(2,iel)+ucons(3,iel))
  p1 = eos3_pr(g_gam1, g_pc1, ucons(2,iel)/ucons(1,iel), &
               ucons(5,iel)/ucons(1,iel), u_conv_l)
  p2 = eos3_pr(g_gam2, g_pc2, ucons(3,iel)/(1.0-ucons(1,iel)),&
               ucons(6,iel)/(1.0-ucons(1,iel)), u_conv_l)
  pr = ucons(1,iel)*p1 + (1.0-ucons(1,iel))*p2
  ncnfl = pr * u_conv_l

  ! right element
  u_conv_r = ucons(4,ier)/(ucons(2,ier)+ucons(3,ier))
  p1 = eos3_pr(g_gam1, g_pc1, ucons(2,ier)/ucons(1,ier), &
               ucons(5,ier)/ucons(1,ier), u_conv_r)
  p2 = eos3_pr(g_gam2, g_pc2, ucons(3,ier)/(1.0-ucons(1,ier)),&
               ucons(6,ier)/(1.0-ucons(1,ier)), u_conv_r)
  pr = ucons(1,ier)*p1 + (1.0-ucons(1,ier))*p2
  ncnfr = pr * u_conv_r

  uf = 0.5 * (u_conv_l+u_conv_r) !lplus*u_conv_l + lminu*u_conv_r

  !print*, ifc, intflux(1), intflux(2), intflux(3), intflux(4), intflux(5), intflux(6)
  !print*, ifc, ncnfl, ncnfr, alp1f, alp2f

  if (iel .gt. 0) then
    rhsel(1,iel) = rhsel(1,iel) - intflux(1) * u_conv_l/lmag
    do ieqn = 2,g_neqns
          rhsel(ieqn,iel) = rhsel(ieqn,iel) - intflux(ieqn)
    end do !ieqn
    !rhsel(1,iel) = rhsel(1,iel) + ucons(1,iel) * uf
    rhsel(5,iel) = rhsel(5,iel) + alp1f * ncnfl
    rhsel(6,iel) = rhsel(6,iel) + alp2f * ncnfl
  end if

  if (ier .lt. (imax+1)) then
    rhsel(1,ier) = rhsel(1,ier) + intflux(1) * u_conv_r/lmag
    do ieqn = 2,g_neqns
          rhsel(ieqn,ier) = rhsel(ieqn,ier) + intflux(ieqn)
    end do !ieqn
    !rhsel(1,ier) = rhsel(1,ier) - ucons(1,ier) * uf
    rhsel(5,ier) = rhsel(5,ier) - alp1f * ncnfr
    rhsel(6,ier) = rhsel(6,ier) - alp2f * ncnfr
  end if

  end do !ifc

  !do iel = 1,imax
  !print*, iel, rhsel(1,iel), rhsel(2,iel), rhsel(3,iel), rhsel(4,iel), rhsel(5,iel), rhsel(6,iel)
  !end do !iel
  !print*, " "

end subroutine flux_p0_mm6eq

!----------------------------------------------------------------------------------------------
!----- Interface pressure correction (hyperbolization)
!----------------------------------------------------------------------------------------------

subroutine get_pugradalpha(ucons, rhsel)

real*8,  intent(in) :: ucons(g_neqns,0:imax+1)

integer :: ie
real*8  :: dx, al1, al2, p1, p2, u_conv, pr, &
           fwddiff,cntdiff,bwddiff, &
           gradal1, rhsel(g_neqns,imax)

  do ie = 1,imax

    dx = 2.0 * (coord(ie+1)-coord(ie))

    u_conv = ucons(4,ie)/(ucons(2,ie)+ucons(3,ie))
    al1 = ucons(1,ie)
    al2 = 1.0 - al1
    p1  = eos3_pr(g_gam1, g_pc1, ucons(2,ie)/al1, ucons(5,ie)/al1, u_conv)
    p2  = eos3_pr(g_gam2, g_pc2, ucons(3,ie)/al2, ucons(6,ie)/al2, u_conv)
    pr  = al1*p1+al2*p2

    fwddiff = ( ucons(1,ie+1) - ucons(1,ie) ) / (0.5*dx)
    bwddiff = ( ucons(1,ie) - ucons(1,ie-1) ) / (0.5*dx)
    cntdiff = (ucons(1,ie+1) - ucons(1,ie-1)) / dx

    gradal1 = mclim(fwddiff,cntdiff,bwddiff)

    rhsel(5,ie) = rhsel(5,ie) + pr*u_conv*gradal1*dx*0.5
    rhsel(6,ie) = rhsel(6,ie) - pr*u_conv*gradal1*dx*0.5

  end do !ie

end subroutine get_pugradalpha

!-------------------------------------------------------------------------------------
!----- 2fluid Lax-Friedrichs flux:
!-------------------------------------------------------------------------------------

subroutine llf_mm6eq(ul, ur, flux, lplus, lminu)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns)

real*8 :: flux(g_neqns), lplus, lminu
real*8 :: al1_l, al1_r, al2_l, al2_r
real*8 :: ffunc_l(g_neqns), ffunc_r(g_neqns), &
          arho1_l,rho1_l,e1_l,a1_l,h1_l,p1_l, &
          arho2_l,rho2_l,e2_l,a2_l,h2_l,p2_l, &
          arho1_r,rho1_r,e1_r,a1_r,h1_r,p1_r, &
          arho2_r,rho2_r,e2_r,a2_r,h2_r,p2_r, &
          rhou_l, u_l, rho_l, p_l, &
          rhou_r, u_r, rho_r, p_r
real*8 :: a1_12,rho1_12,al1_12, &
          a2_12,rho2_12,al2_12, &
          rho_12,ac_12
  
real*8 :: lambda

  flux(:) = 0.0
  
  ! ul
  al1_l   = ul(1)
  arho1_l = ul(2)
  arho2_l = ul(3)
  rhou_l  = ul(4)
  e1_l    = ul(5)
  e2_l    = ul(6)
  al2_l   = 1.0 - al1_l
  
  rho_l  = arho1_l + arho2_l
  rho1_l = arho1_l / al1_l
  rho2_l = arho2_l / al2_l
  u_l    = rhou_l / rho_l
  p1_l   = eos3_pr(g_gam1, g_pc1, rho1_l, (e1_l/al1_l), u_l)
  p2_l   = eos3_pr(g_gam2, g_pc2, rho2_l, (e2_l/al2_l), u_l)
  p_l    = al1_l*p1_l+al2_l*p2_l
  a1_l   = eos3_ss(g_gam1, g_pc1, rho1_l, p1_l)
  h1_l   = e1_l + al1_l*p1_l
  a2_l   = eos3_ss(g_gam2, g_pc2, rho2_l, p2_l)
  h2_l   = e2_l + al2_l*p2_l
  
  ! ur
  al1_r   = ur(1)
  arho1_r = ur(2)
  arho2_r = ur(3)
  rhou_r  = ur(4)
  e1_r    = ur(5)
  e2_r    = ur(6)
  al2_r   = 1.0 - al1_r
  
  rho_r  = arho1_r + arho2_r
  rho1_r = arho1_r / al1_r
  rho2_r = arho2_r / al2_r
  u_r    = rhou_r / rho_r
  p1_r   = eos3_pr(g_gam1, g_pc1, rho1_r, (e1_r/al1_r), u_r)
  p2_r   = eos3_pr(g_gam2, g_pc2, rho2_r, (e2_r/al2_r), u_r)
  p_r    = al1_r*p1_r+al2_r*p2_r
  a1_r   = eos3_ss(g_gam1, g_pc1, rho1_r, p1_r)
  h1_r   = e1_r + al1_r*p1_r
  a2_r   = eos3_ss(g_gam2, g_pc2, rho2_r, p2_r)
  h2_r   = e2_r + al2_r*p2_r
  
  rho_12  = 0.5*(rho_l + rho_r)
  rho1_12 = 0.5*(rho1_l + rho1_r)
  a1_12   = 0.5*(a1_l + a1_r)
  rho2_12 = 0.5*(rho2_l + rho2_r)
  a2_12   = 0.5*(a2_l + a2_r)
  al1_12  = 0.5*(al1_l + al1_r)
  al2_12  = 0.5*(al2_l + al2_r)
  
  ! numerical speed of sound choice:
  ac_12 = dsqrt( (  al1_12*rho1_12*a1_12*a1_12 &
                  + al2_12*rho2_12*a2_12*a2_12 ) / rho_12 )
  
  lambda = ac_12 + max(dabs(u_l),dabs(u_r));

  ffunc_l(1) = u_l * al1_l
  ffunc_l(2) = u_l * arho1_l
  ffunc_l(3) = u_l * arho2_l
  ffunc_l(4) = u_l * rhou_l + p_l
  ffunc_l(5) = u_l * h1_l
  ffunc_l(6) = u_l * h2_l

  ffunc_r(1) = u_r * al1_r
  ffunc_r(2) = u_r * arho1_r
  ffunc_r(3) = u_r * arho2_r
  ffunc_r(4) = u_r * rhou_r + p_r
  ffunc_r(5) = u_r * h1_r
  ffunc_r(6) = u_r * h2_r
  
  flux(1) = 0.5 * ( ffunc_l(1)+ffunc_r(1) - lambda*(ur(1)-ul(1)) )
  flux(2) = 0.5 * ( ffunc_l(2)+ffunc_r(2) - lambda*(ur(2)-ul(2)) )
  flux(3) = 0.5 * ( ffunc_l(3)+ffunc_r(3) - lambda*(ur(3)-ul(3)) )
  flux(4) = 0.5 * ( ffunc_l(4)+ffunc_r(4) - lambda*(ur(4)-ul(4)) )
  flux(5) = 0.5 * ( ffunc_l(5)+ffunc_r(5) - lambda*(ur(5)-ul(5)) )
  flux(6) = 0.5 * ( ffunc_l(6)+ffunc_r(6) - lambda*(ur(6)-ul(6)) )

  lplus = 0.5
  lminu = 0.5

end subroutine llf_mm6eq

!-------------------------------------------------------------------------------------
!----- 2fluid AUSM+UP:
!-------------------------------------------------------------------------------------

subroutine ausmplus_mm6eq(ul, ur, flux, lambda_plus, lambda_minu, lambda_mag)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns)

real*8 :: flux(g_neqns)
real*8 :: al1_l, al1_r, al2_l, al2_r
real*8 :: arho1_l,rho1_l,e1_l,a1_l,h1_l,p1_l, &
          arho2_l,rho2_l,e2_l,a2_l,h2_l,p2_l, &
          arho1_r,rho1_r,e1_r,a1_r,h1_r,p1_r, &
          arho2_r,rho2_r,e2_r,a2_r,h2_r,p2_r, &
          rhou_l, u_l, m_l, rho_l, p_l, &
          rhou_r, u_r, m_r, rho_r, p_r
real*8 :: a1_12,rho1_12,al1_12,mdot1_12, &
          a2_12,rho2_12,al2_12,mdot2_12, &
          f_a,rho_12,m_12,p_12,ac_12!,m_p,p_u
real*8 :: msplus_l(3),msplus_r(3),msminu_l(3),msminu_r(3)
real*8 :: psplus_l,psplus_r,psminu_l,psminu_r
!real*8 :: temp,temp1,temp2,num,den
  
real*8 :: lambda,lambda_plus, lambda_minu, lambda_mag
  
!real*8 :: k_p, k_u

  flux(:) = 0.0
  
  ! ul
  al1_l   = ul(1)
  arho1_l = ul(2)
  arho2_l = ul(3)
  rhou_l  = ul(4)
  e1_l    = ul(5)
  e2_l    = ul(6)
  al2_l   = 1.0 - al1_l
  
  rho_l  = arho1_l + arho2_l
  rho1_l = arho1_l / al1_l
  rho2_l = arho2_l / al2_l
  u_l    = rhou_l / rho_l
  p1_l   = eos3_pr(g_gam1, g_pc1, rho1_l, (e1_l/al1_l), u_l)
  p2_l   = eos3_pr(g_gam2, g_pc2, rho2_l, (e2_l/al2_l), u_l)
  p_l    = al1_l*p1_l+al2_l*p2_l
  a1_l   = eos3_ss(g_gam1, g_pc1, rho1_l, p1_l)
  h1_l   = e1_l + al1_l*p1_l
  a2_l   = eos3_ss(g_gam2, g_pc2, rho2_l, p2_l)
  h2_l   = e2_l + al2_l*p2_l
  
  ! ur
  al1_r   = ur(1)
  arho1_r = ur(2)
  arho2_r = ur(3)
  rhou_r  = ur(4)
  e1_r    = ur(5)
  e2_r    = ur(6)
  al2_r   = 1.0 - al1_r
  
  rho_r  = arho1_r + arho2_r
  rho1_r = arho1_r / al1_r
  rho2_r = arho2_r / al2_r
  u_r    = rhou_r / rho_r
  p1_r   = eos3_pr(g_gam1, g_pc1, rho1_r, (e1_r/al1_r), u_r)
  p2_r   = eos3_pr(g_gam2, g_pc2, rho2_r, (e2_r/al2_r), u_r)
  p_r    = al1_r*p1_r+al2_r*p2_r
  a1_r   = eos3_ss(g_gam1, g_pc1, rho1_r, p1_r)
  h1_r   = e1_r + al1_r*p1_r
  a2_r   = eos3_ss(g_gam2, g_pc2, rho2_r, p2_r)
  h2_r   = e2_r + al2_r*p2_r
  
  rho_12  = 0.5*(rho_l + rho_r)
  rho1_12 = 0.5*(rho1_l + rho1_r)
  a1_12   = 0.5*(a1_l + a1_r)
  rho2_12 = 0.5*(rho2_l + rho2_r)
  a2_12   = 0.5*(a2_l + a2_r)
  al1_12  = 0.5*(al1_l + al1_r)
  al2_12  = 0.5*(al2_l + al2_r)
  
  ! numerical speed of sound choice:
  ac_12 = dsqrt( (  al1_12*rho1_12*a1_12*a1_12 &
                  + al2_12*rho2_12*a2_12*a2_12 ) / rho_12 )
  
  m_l = u_l/ac_12
  m_r = u_r/ac_12
  
  ! all-speed scaling:
  !mbar2 = c05 * (vn_l*vn_l + vn_r*vn_r)/(ac_12*ac_12)
  !umag_fs = sqrt(g_freestream.u*g_freestream.u + g_freestream.v*g_freestream.v)
  !m_0   = sqrt(min(1.0, max(mbar2, umag_fs/ac_12)))
  f_a = 1.0 !m_0 * (2.0 - m_0)
  
  ! split mach number functions
  call splitmach_as(f_a,m_l,msplus_l,msminu_l,psplus_l,psminu_l)
  call splitmach_as(f_a,m_r,msplus_r,msminu_r,psplus_r,psminu_r)
  
  ! "p"
  !temp  = c10 - (c05*(vn_l*vn_l + vn_r*vn_r)/(ac_12*ac_12))
  !m_p   = -k_p* (max(temp,c00))* (p_r-p_l) / (f_a* rho_12*ac_12*ac_12)
  m_12 = msplus_l(3) + msminu_r(3)! + m_p
  
  ! "u"
  !p_u   = -k_u* (c10-psplus_l* psminu_r)* f_a* rho_12* vnrel* (vn_r-vn_l)
  p_12 = psplus_l*p_l + psminu_r*p_r! + p_u
  
  lambda = ac_12 * m_12
  
  lambda_plus = 0.5 * (lambda + dabs(lambda))
  lambda_minu = 0.5 * (lambda - dabs(lambda))
  
  flux(1) = lambda_plus*al1_l   + lambda_minu*al1_r
  flux(2) = lambda_plus*arho1_l + lambda_minu*arho1_r
  flux(3) = lambda_plus*arho2_l + lambda_minu*arho2_r
  flux(4) = lambda_plus*rhou_l  + lambda_minu*rhou_r + p_12
  flux(5) = lambda_plus*(h1_l)  + lambda_minu*(h1_r)
  flux(6) = lambda_plus*(h2_l)  + lambda_minu*(h2_r)

  lambda_mag = dabs(lambda) + 1.d-16

  lambda_plus = lambda_plus/(lambda_mag) 
  lambda_minu = lambda_minu/(lambda_mag)

end subroutine ausmplus_mm6eq

!-------------------------------------------------------------------------------------
!----- Split Mach polynomials for AUSM+UP:
!-------------------------------------------------------------------------------------
subroutine splitmach_as(fa, mach, msplus, msminu, psplus, psminu)

real*8, intent(in) :: fa, mach

real*8             :: alph_fa, msplus(3), msminu(3), psplus, psminu

    msplus(1) = 0.5d0*(mach + dabs(mach))
    msminu(1) = 0.5d0*(mach - dabs(mach))

    msplus(2) = +0.25d0*(mach + 1.0d0)*(mach + 1.0d0)
    msminu(2) = -0.25d0*(mach - 1.0d0)*(mach - 1.0d0)

    alph_fa = (3.d0/16.d0) * (-4.d0 + 5.d0*fa*fa)

    if (dabs(mach) .ge. 1.d0) then
    
        msplus(3) = msplus(1)
        msminu(3) = msminu(1)
        psplus    = msplus(1)/mach
        psminu    = msminu(1)/mach
   
    else
    
        msplus(3) = msplus(2)* (1.d0 - 2.d0*msminu(2))
        msminu(3) = msminu(2)* (1.d0 + 2.d0*msplus(2))
        psplus    = msplus(2)* ((+2.d0 - mach) - (16.d0 * alph_fa)*mach*msminu(2))
        psminu    = msminu(2)* ((-2.d0 - mach) + (16.d0 * alph_fa)*mach*msplus(2))

    end if

end subroutine splitmach_as

!----------------------------------------------------------------------------------------------
!----- Boundary conditions:
!----------------------------------------------------------------------------------------------

subroutine get_bc_mm6eq(ucons)

real*8  :: ucons(g_neqns,0:imax+1)
real*8  :: p1, p2, rho1, rho2, u_conv

  !----- left boundary

  if (g_lbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     ucons(:,0) = ucons(:,1)

  else if (g_lbflag .eq. 1) then
     !--- supersonic inflow
     ucons(1,0) = alpha1_fs
     ucons(2,0) = alpha1_fs       * rho1_fs
     ucons(3,0) = (1.0-alpha1_fs) * rho2_fs
     ucons(4,0) = (ucons(3,0)+ucons(2,0)) * u_fs
     ucons(5,0) = alpha1_fs       * eos3_rhoe(g_gam1, g_pc1, pr1_fs, rho1_fs, u_fs)
     ucons(6,0) = (1.0-alpha1_fs) * eos3_rhoe(g_gam2, g_pc2, pr2_fs, rho2_fs, u_fs)

  else if (g_lbflag .eq. 2) then
     !--- subsonic inflow
     u_conv = ucons(4,1)/(ucons(2,1)+ucons(3,1))
     p1  = eos3_pr(g_gam1, g_pc1, ucons(2,1)/ucons(1,1), &
                   ucons(5,1)/ucons(1,1), u_conv)
     p2  = eos3_pr(g_gam2, g_pc2, ucons(3,1)/(1.0-ucons(1,1)),&
                   ucons(6,1)/(1.0-ucons(1,1)), u_conv)
     ucons(1,0) = alpha1_fs
     ucons(2,0) = alpha1_fs       * rho1_fs
     ucons(3,0) = (1.0-alpha1_fs) * rho2_fs
     ucons(4,0) = (ucons(3,0)+ucons(2,0)) * u_fs
     ucons(5,0) = alpha1_fs       * eos3_rhoe(g_gam1, g_pc1, p1, rho1_fs, u_fs)
     ucons(6,0) = (1.0-alpha1_fs) * eos3_rhoe(g_gam2, g_pc2, p2, rho2_fs, u_fs)

  else
     write(*,*) "BC-type not set for flag ", g_lbflag

  end if

  !----- right boundary

  if (g_rbflag .eq. 0) then
     !--- extrapolation / supersonic outflow
     ucons(:,imax+1) = ucons(:,imax)

  else if (g_rbflag .eq. 1) then
     !--- supersonic inflow
     ucons(1,imax+1) = alpha1_fs
     ucons(2,imax+1) = alpha1_fs       * rho1_fs
     ucons(3,imax+1) = (1.0-alpha1_fs) * rho2_fs
     ucons(4,imax+1) = (ucons(3,imax+1)+ucons(2,imax+1)) * u_fs
     ucons(5,imax+1) = alpha1_fs       * eos3_rhoe(g_gam1, g_pc1, pr1_fs, rho1_fs, u_fs)
     ucons(6,imax+1) = (1.0-alpha1_fs) * eos3_rhoe(g_gam2, g_pc2, pr2_fs, rho2_fs, u_fs)

  else if (g_rbflag .eq. 2) then
     !--- subsonic inflow
     u_conv = ucons(4,imax)/(ucons(2,imax)+ucons(3,imax))
     p1  = eos3_pr(g_gam1, g_pc1, ucons(2,imax)/ucons(1,imax), &
                   ucons(5,imax)/ucons(1,imax), u_conv)
     p2  = eos3_pr(g_gam2, g_pc2, ucons(3,imax)/(1.0-ucons(1,imax)),&
                   ucons(6,imax)/(1.0-ucons(1,imax)), u_conv)
     ucons(1,imax+1) = alpha1_fs
     ucons(2,imax+1) = alpha1_fs       * rho1_fs
     ucons(3,imax+1) = (1.0-alpha1_fs) * rho2_fs
     ucons(4,imax+1) = (ucons(3,imax+1)+ucons(2,imax+1)) * u_fs
     ucons(5,imax+1) = alpha1_fs       * eos3_rhoe(g_gam1, g_pc1, p1, rho1_fs, u_fs)
     ucons(6,imax+1) = (1.0-alpha1_fs) * eos3_rhoe(g_gam2, g_pc2, p2, rho2_fs, u_fs)

  else if (g_rbflag .eq. 3) then
     !--- subsonic outflow
     u_conv = ucons(4,imax)/(ucons(2,imax)+ucons(3,imax))
     rho1 = ucons(2,imax) / ucons(1,imax)
     rho2 = ucons(3,imax) / (1.0-ucons(1,imax))
     ucons(1,imax+1) = ucons(1,imax)
     ucons(2,imax+1) = ucons(2,imax)
     ucons(3,imax+1) = ucons(3,imax)
     ucons(4,imax+1) = ucons(4,imax)
     ucons(5,imax+1) = ucons(1,imax)       * eos3_rhoe(g_gam1, g_pc1, pr1_fs, rho1, u_conv)
     ucons(6,imax+1) = (1.0-ucons(1,imax)) * eos3_rhoe(g_gam2, g_pc2, pr2_fs, rho2, u_conv)

  else
     write(*,*) "BC-type not set for flag ", g_rbflag

  end if

end subroutine get_bc_mm6eq

!----------------------------------------------------------------------------------------------

END MODULE rhs_flux_mm6eq
