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

integer :: ifc, iel, ier
real*8  :: ul(neqns), ur(neqns), &
           intflux(neqns), intpflux(2), rhsel(neqns,imax), &
           alpha1_l, alpha2_l, alpha1_r, alpha2_r, &
           pflux1_l, pflux2_l, pflux1_r, pflux2_r

real*8, intent(in) :: uprim(ndof,neqns,0:imax+1), ucons(neqns,0:imax+1)

  do ifc = 1,imax+1

  iel = ifc - 1
  ier = ifc

  ul(:) = ucons(:,iel)
  ur(:) = ucons(:,ier)

  !--- Conservative fluxes

  call ausmplus_mm6eq(ul, ur, intflux)

  if (iel .gt. 0) then
          rhsel(1,iel) = rhsel(1,iel) - intflux(1)
          rhsel(2,iel) = rhsel(2,iel) - intflux(2)
          rhsel(3,iel) = rhsel(3,iel) - intflux(3)
          rhsel(4,iel) = rhsel(4,iel) - intflux(4)
          rhsel(5,iel) = rhsel(5,iel) - intflux(5)
          rhsel(6,iel) = rhsel(6,iel) - intflux(6)
  end if

  if (ier .lt. (imax+1)) then
          rhsel(1,ier) = rhsel(1,ier) + intflux(1)
          rhsel(2,ier) = rhsel(2,ier) + intflux(2)
          rhsel(3,ier) = rhsel(3,ier) + intflux(3)
          rhsel(4,ier) = rhsel(4,ier) + intflux(4)
          rhsel(5,ier) = rhsel(5,ier) + intflux(5)
          rhsel(6,ier) = rhsel(6,ier) + intflux(6)
  end if

  !--- Non-conservative fluxes
  !--- incomplete

  end do !ifc

end subroutine flux_p0_mm6eq

!-------------------------------------------------------------------------------------
!----- 2fluid AUSM+UP:
!-------------------------------------------------------------------------------------

subroutine ausmplus_mm6eq(ul, ur, flux)

real*8, intent(in) :: ul(neqns), ur(neqns)

real*8 :: flux(neqns)
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
  
real*8 :: lambda,lambda_plus, lambda_minu
  
!real*8 :: k_p, k_u
  
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
  
  !p_12 = max(p_12, 1.0e-6*min(p_l,p_r))
  
  lambda = ac_12 * m_12
  
  lambda_plus = 0.5 * (lambda + dabs(lambda))
  lambda_minu = 0.5 * (lambda - dabs(lambda))
  
  flux(1) = lambda_plus*al1_l   + lambda_minu*al1_r
  flux(2) = lambda_plus*arho1_l + lambda_minu*arho1_r
  flux(3) = lambda_plus*arho2_l + lambda_minu*arho2_r
  flux(4) = lambda_plus*rhou_l  + lambda_minu*rhou_r + p_12
  flux(5) = lambda_plus*(h1_l)  + lambda_minu*(h1_r)
  flux(6) = lambda_plus*(h2_l)  + lambda_minu*(h2_r)

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

real*8  :: ucons(neqns,0:imax+1)

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
     ucons(4,imax+1) = (ucons(3,0)+ucons(2,0)) * u_fs
     ucons(5,imax+1) = alpha1_fs       * eos3_rhoe(g_gam1, g_pc1, pr1_fs, rho1_fs, u_fs)
     ucons(6,imax+1) = (1.0-alpha1_fs) * eos3_rhoe(g_gam2, g_pc2, pr2_fs, rho2_fs, u_fs)

  else
     write(*,*) "BC-type not set for flag ", g_rbflag

  end if

end subroutine get_bc_mm6eq

!----------------------------------------------------------------------------------------------

END MODULE rhs_flux_mm6eq
