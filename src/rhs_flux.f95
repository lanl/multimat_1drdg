!!---------------------------------------------------------------------------------------
!!----- Velocity non-equilibrium isothermal Multiphase (Flux computation module)
!!----- by
!!----- Aditya K Pandare
!!---------------------------------------------------------------------------------------

MODULE rhs_flux

USE glob_var
USE eos

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Interface pressure correction (hyperbolization)
!----------------------------------------------------------------------------------------------

subroutine get_pgradalpha(uprim, rhsel)

real*8,  intent(in) :: uprim(ndof,g_neqns,0:imax+1)

integer :: ielem
real*8  :: dx, alpha1, alpha2, u1, u2, pr, rho1, rho2, &
           fwddiff,cntdiff,bwddiff, &
           gradalpha1, urel, temp, temp1, pint, deltap, rhsel(g_neqns,imax)

        do ielem = 1,imax

           dx = 2.0 * (coord(ielem+1)-coord(ielem))

           alpha1 = uprim(1,1,ielem)
           pr = uprim(1,2,ielem)
           u1 = uprim(1,3,ielem)
           u2 = uprim(1,4,ielem)
           alpha2 = 1.0 - alpha1

           rho1 = eos1_density(pr)
           rho2 = eos2_density(pr)

           fwddiff = ( uprim(1,1,ielem+1) - uprim(1,1,ielem) ) / (dx/2.0)
           bwddiff = ( uprim(1,1,ielem) - uprim(1,1,ielem-1) ) / (dx/2.0)
           cntdiff = (uprim(1,1,ielem+1) - uprim(1,1,ielem-1)) / dx

           gradalpha1 = mclim(fwddiff,cntdiff,bwddiff)

           urel = u1-u2
           temp  = alpha1*alpha2*rho1*rho2
           temp1 = alpha1*rho2 + alpha2*rho1
           deltap = ipsig * temp/temp1 * (urel*urel)
           pint = pr - min(deltap,0.01*pr)

           rhsel(2,ielem) = rhsel(2,ielem) + pint*gradalpha1*dx/2.0
           rhsel(4,ielem) = rhsel(4,ielem) - pint*gradalpha1*dx/2.0

        end do !ielem

end subroutine get_pgradalpha

!----------------------------------------------------------------------------------------------
!----- Interface pressure correction (hyperbolization) P0P1
!----------------------------------------------------------------------------------------------

subroutine get_pgradalpha_p0p1(uprim, rhsel)

real*8,  intent(in) :: uprim(ndof,g_neqns,0:imax+1)

integer :: ielem
real*8  :: dx, alpha1, alpha2, u1, u2, pr, rho1, rho2, &
           gradalpha1, urel, temp, temp1, pint, deltap, rhsel(g_neqns,imax)

        do ielem = 1,imax

           dx = 2.0 * (coord(ielem+1)-coord(ielem))

           alpha1 = uprim(1,1,ielem)
           pr = uprim(1,2,ielem)
           u1 = uprim(1,3,ielem)
           u2 = uprim(1,4,ielem)
           alpha2 = 1.0 - alpha1

           rho1 = eos1_density(pr)
           rho2 = eos2_density(pr)

           gradalpha1 = uprim(2,1,ielem)

           urel = u1-u2
           temp  = alpha1*alpha2*rho1*rho2
           temp1 = alpha1*rho2 + alpha2*rho1
           deltap = ipsig * temp/temp1 * (urel*urel)
           pint = pr - min(deltap,0.01*pr)

           rhsel(2,ielem) = rhsel(2,ielem) + pint*gradalpha1*dx/2.0
           rhsel(4,ielem) = rhsel(4,ielem) - pint*gradalpha1*dx/2.0

        end do !ielem

end subroutine get_pgradalpha_p0p1

!----------------------------------------------------------------------------------------------
!----- Advective-flux contribution to RHS:
!----------------------------------------------------------------------------------------------

subroutine flux_p0(uprim, ucons, rhsel)

integer :: ifc, iel, ier
real*8  :: ul(g_neqns), ur(g_neqns), &
           vl(g_neqns), vr(g_neqns), &
           intflux(g_neqns), intpflux(2), rhsel(g_neqns,imax), &
           alpha1_l, alpha2_l, alpha1_r, alpha2_r, &
           pflux1_l, pflux2_l, pflux1_r, pflux2_r

real*8, intent(in) :: uprim(ndof,g_neqns,0:imax+1), ucons(g_neqns,0:imax+1)

        do ifc = 1,imax+1

        iel = ifc - 1
        ier = ifc

        ul(:) = ucons(:,iel)
        ur(:) = ucons(:,ier)

        vl(:) = uprim(1,:,iel)
        vr(:) = uprim(1,:,ier)

        if (i_flux .eq. 1) then
           call LDFSS_real(ul, ur, vl, vr, intflux, intpflux)
        elseif (i_flux .eq. 2) then
           call LDFSS_psreal(ul, ur, vl, vr, intflux, intpflux)
        else
           print*, "Invalid flux scheme."
           stop
        endif

        alpha1_l = vl(1)
        alpha2_l = 1.0 - alpha1_l
        alpha1_r = vr(1)
        alpha2_r = 1.0 - alpha1_r

        ! nonconservative treatment of pressure-flux
        pflux1_l = alpha1_l * intpflux(1) 
        pflux2_l = alpha2_l * intpflux(2) 
        pflux1_r = alpha1_r * intpflux(1) 
        pflux2_r = alpha2_r * intpflux(2) 

        if (iel .gt. 0) then
                rhsel(1,iel) = rhsel(1,iel) - (intflux(1)         )
                rhsel(2,iel) = rhsel(2,iel) - (intflux(2)+pflux1_l)
                rhsel(3,iel) = rhsel(3,iel) - (intflux(3)         )
                rhsel(4,iel) = rhsel(4,iel) - (intflux(4)+pflux2_l)
        end if

        if (ier .lt. (imax+1)) then
                rhsel(1,ier) = rhsel(1,ier) + (intflux(1)         )
                rhsel(2,ier) = rhsel(2,ier) + (intflux(2)+pflux1_r)
                rhsel(3,ier) = rhsel(3,ier) + (intflux(3)         )
                rhsel(4,ier) = rhsel(4,ier) + (intflux(4)+pflux2_r)
        end if

        end do !ifc

end subroutine flux_p0

!----------------------------------------------------------------------------------------------
!----- Advective-flux contribution to RHS second order:
!----------------------------------------------------------------------------------------------

subroutine flux_p0p1(uprim, ucons, rhsel)

integer :: ifc, iel, ier
real*8  :: ul(g_neqns), ur(g_neqns), &
           vl(g_neqns), vr(g_neqns), &
           dxl, dxr, intflux(g_neqns), intpflux(2), rhsel(g_neqns,imax), &
           alpha1_l, alpha2_l, alpha1_r, alpha2_r, &
           pflux1_l, pflux2_l, pflux1_r, pflux2_r

real*8, intent(in) :: uprim(ndof,g_neqns,0:imax+1), ucons(g_neqns,0:imax+1)

        call reconstruct_uprim(uprim)

        do ifc = 1,imax+1

        iel = ifc - 1
        ier = ifc

        dxl = 0.5 * (coord(ifc)-coord(ifc-1))
        dxr = 0.5 * (coord(ifc)-coord(ifc+1))

        !--- primitve variables to be reconstructed
        vl(1) = uprim(1,1,iel) + dxl * uprim(2,1,iel)
        vl(2) = uprim(1,2,iel) + dxl * uprim(2,2,iel)
        vl(3) = uprim(1,3,iel) + dxl * uprim(2,3,iel)
        vl(4) = uprim(1,4,iel) + dxl * uprim(2,4,iel)

        vr(1) = uprim(1,1,ier) + dxr * uprim(2,1,ier)
        vr(2) = uprim(1,2,ier) + dxr * uprim(2,2,ier)
        vr(3) = uprim(1,3,ier) + dxr * uprim(2,3,ier)
        vr(4) = uprim(1,4,ier) + dxr * uprim(2,4,ier)

        ul(1) = vl(1)*eos1_density(vl(2))
        ul(2) = ul(1)*vl(3)
        ul(3) = (1.0-vl(1))*eos2_density(vl(2))
        ul(4) = ul(3)*vl(4)

        ur(1) = vr(1)*eos1_density(vr(2))
        ur(2) = ur(1)*vr(3)
        ur(3) = (1.0-vr(1))*eos2_density(vr(2))
        ur(4) = ur(3)*vr(4)

        if (i_flux .eq. 1) then
           call LDFSS_real(ul, ur, vl, vr, intflux, intpflux)
        elseif (i_flux .eq. 2) then
           call LDFSS_psreal(ul, ur, vl, vr, intflux, intpflux)
        else
           print*, "Invalid flux scheme."
           stop
        endif

        alpha1_l = vl(1)
        alpha2_l = 1.0 - alpha1_l
        alpha1_r = vr(1)
        alpha2_r = 1.0 - alpha1_r

        ! nonconservative treatment of pressure-flux
        pflux1_l = alpha1_l * intpflux(1) 
        pflux2_l = alpha2_l * intpflux(2) 
        pflux1_r = alpha1_r * intpflux(1) 
        pflux2_r = alpha2_r * intpflux(2) 

        if (iel .gt. 0) then
                rhsel(1,iel) = rhsel(1,iel) - (intflux(1)         )
                rhsel(2,iel) = rhsel(2,iel) - (intflux(2)+pflux1_l)
                rhsel(3,iel) = rhsel(3,iel) - (intflux(3)         )
                rhsel(4,iel) = rhsel(4,iel) - (intflux(4)+pflux2_l)
        end if

        if (ier .lt. (imax+1)) then
                rhsel(1,ier) = rhsel(1,ier) + (intflux(1)         )
                rhsel(2,ier) = rhsel(2,ier) + (intflux(2)+pflux1_r)
                rhsel(3,ier) = rhsel(3,ier) + (intflux(3)         )
                rhsel(4,ier) = rhsel(4,ier) + (intflux(4)+pflux2_r)
        end if

        end do !ifc

end subroutine flux_p0p1

!----------------------------------------------------------------------------------------------
!----- All-speed + real-gas LDFSS flux by Edwards:
!----------------------------------------------------------------------------------------------

subroutine LDFSS_real(ul, ur, vl, vr, intflux, intpflux)

real*8 :: arho1_l, arhou1_l, u1_l, alpha1_l, a1_l, rho1_l, m1_l, &
          arho2_l, arhou2_l, u2_l, alpha2_l, a2_l, rho2_l, m2_l, &
          arho1_r, arhou1_r, u1_r, alpha1_r, a1_r, rho1_r, m1_r, &
          arho2_r, arhou2_r, u2_r, alpha2_r, a2_r, rho2_r, m2_r, &
          a1_12, rho1_12, &
          a2_12, rho2_12, &
          p_l, p_r, ac_12, deltap, p_12, &
          mplus, pplus, mplus_12, mdotplus, &
          mminu, pminu, mminu_12, mdotminu, &
          vrel, disp, mcorr_12, macht

real*8 :: intflux(g_neqns), intpflux(2)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), vl(g_neqns), vr(g_neqns)

        !----- Left state
        arho1_l  = ul(1)
        arhou1_l = ul(2)
        arho2_l  = ul(3)
        arhou2_l = ul(4)

        alpha1_l = vl(1)
        p_l      = vl(2)
        u1_l     = vl(3)
        u2_l     = vl(4)
        alpha2_l = 1.0 - alpha1_l

        rho1_l = arho1_l / alpha1_l
        rho2_l = arho2_l / alpha2_l

        a1_l = agas
        a2_l = eos2_ssound(rho2_l)

        !----- Right state
        arho1_r  = ur(1)
        arhou1_r = ur(2)
        arho2_r  = ur(3)
        arhou2_r = ur(4)

        alpha1_r = vr(1)
        p_r      = vr(2)
        u1_r     = vr(3)
        u2_r     = vr(4)
        alpha2_r = 1.0 - alpha1_r

        rho1_r = arho1_r / alpha1_r
        rho2_r = arho2_r / alpha2_r

        a1_r = agas
        a2_r = eos2_ssound(rho2_r)

        ! interface speed of sound
        a1_12 = 0.5 * (a1_l+a1_r)
        a2_12 = 0.5 * (a2_l+a2_r)

        deltap = p_l - p_r
        rho1_12 = 0.5 * (rho1_l+rho1_r)
        rho2_12 = 0.5 * (rho2_l+rho2_r)

        ! common speed of sound choice
        ac_12 = 0.5 * (a1_12+a2_12)
        vrel = max(dabs(u1_l-u2_l), dabs(u1_r-u2_r))

        ! mach numbers
        m1_l = u1_l/ac_12
        m1_r = u1_r/ac_12
        m2_l = u2_l/ac_12
        m2_r = u2_r/ac_12

        !--- phase 1
        !--------------------------------------------------------

        ! split mach numbers
        if (dabs(m1_l) .lt. 1.0) then
           mplus = 0.25 * (m1_l + 1.0) * (m1_l + 1.0)
           pplus = 0.5  * (1.0 + m1_l)
        else
           macht = 0.5 * (m1_l + dabs(m1_l))
           mplus = macht
           pplus = macht / m1_l
        end if

        if (dabs(m1_r) .lt. 1.0) then
           mminu = - 0.25 * (m1_r - 1.0) * (m1_r - 1.0)
           pminu = 0.5  * (1.0 - m1_r)
        else
           macht = 0.5 * (m1_r - dabs(m1_r))
           mminu = macht
           pminu = macht / m1_r
        end if

        ! correction
        if ((dabs(m1_l)<1.0) .and. (dabs(m1_r)<1.0)) then
           mcorr_12 = 0.25 * ( ( dsqrt(0.5*(m1_l*m1_l + m1_r*m1_r)) - 1.0 )**2.0 )
        else
           mcorr_12 = 0.0
        end if

        mplus_12 = mcorr_12 * (1.0 - (deltap + dabs(deltap)) / (2.0*rho1_l*ac_12*ac_12))
        mminu_12 = mcorr_12 * (1.0 + (deltap - dabs(deltap)) / (2.0*rho1_r*ac_12*ac_12))

        disp = vrel * max(alpha1_l,alpha1_r) * (arho1_r-arho1_l) * 10.0

        ! mass fluxes
        mdotplus = arho1_l * ac_12 * (mplus - mplus_12) !- disp
        mdotminu = arho1_r * ac_12 * (mminu + mminu_12) !- disp

        ! convection flux
        intflux(1) = mdotplus        + mdotminu
        intflux(2) = mdotplus * u1_l + mdotminu * u1_r

        ! interface pressure
        p_12 =  0.5 * (p_l+p_r) &
              + 0.5 * (deltap) * (pplus-pminu) &
              + rho1_12*ac_12*ac_12 * (pplus+pminu-1.0)

        ! pressure flux
        intpflux(1) = p_12

        !--- phase 2
        !--------------------------------------------------------

        ! split mach numbers
        if (dabs(m2_l) .lt. 1.0) then
           mplus = 0.25 * (m2_l + 1.0) * (m2_l + 1.0)
           pplus = 0.5  * (1.0 + m2_l)
        else
           macht = 0.5 * (m2_l + dabs(m2_l))
           mplus = macht
           pplus = macht / m2_l
        end if

        if (dabs(m2_r) .lt. 1.0) then
           mminu = - 0.25 * (m2_r - 1.0) * (m2_r - 1.0)
           pminu = 0.5  * (1.0 - m2_r)
        else
           macht = 0.5 * (m2_r - dabs(m2_r))
           mminu = macht
           pminu = macht / m2_r
        end if

        ! correction
        if ((dabs(m2_l)<1.0) .and. (dabs(m2_r)<1.0)) then
           mcorr_12 = 0.25 * ( ( dsqrt(0.5*(m2_l*m2_l + m2_r*m2_r)) - 1.0 )**2.0 )
        else
           mcorr_12 = 0.0
        end if

        mplus_12 = mcorr_12 * (1.0 - (deltap + dabs(deltap)) / (2.0*rho2_l*ac_12*ac_12))
        mminu_12 = mcorr_12 * (1.0 + (deltap - dabs(deltap)) / (2.0*rho2_r*ac_12*ac_12))

        disp = vrel * max(alpha2_l,alpha2_r) * (arho2_r-arho2_l) * 10.0

        ! mass fluxes
        mdotplus = arho2_l * ac_12 * (mplus - mplus_12) !- disp
        mdotminu = arho2_r * ac_12 * (mminu + mminu_12) !- disp

        ! convection flux
        intflux(3) = mdotplus        + mdotminu
        intflux(4) = mdotplus * u2_l + mdotminu * u2_r

        ! interface pressure
        p_12 =  0.5 * (p_l+p_r) &
              + 0.5 * (deltap) * (pplus-pminu) &
              + rho2_12*ac_12*ac_12 * (pplus+pminu-1.0)

        ! pressure flux
        intpflux(2) = p_12

end subroutine LDFSS_real

!----------------------------------------------------------------------------------------------
!----- All-speed + real-gas LDFSS flux by Edwards:
!----- a single pressure flux is calculated using mixture mach numbers and densities
!----------------------------------------------------------------------------------------------

subroutine LDFSS_psreal(ul, ur, vl, vr, intflux, intpflux)

real*8 :: arho1_l, arhou1_l, u1_l, alpha1_l, a1_l, rho1_l, m1_l, &
          arho2_l, arhou2_l, u2_l, alpha2_l, a2_l, rho2_l, m2_l, &
          arho1_r, arhou1_r, u1_r, alpha1_r, a1_r, rho1_r, m1_r, &
          arho2_r, arhou2_r, u2_r, alpha2_r, a2_r, rho2_r, m2_r, &
          a1_12, rho1_12, &
          a2_12, rho2_12, &
          p_l, p_r, ac_12, deltap, p_12, &
          mplus, pplus, mplus_12, mdotplus, &
          mminu, pminu, mminu_12, mdotminu, &
          vrel, disp, mcorr_12, macht, &
          vmix_l, vmix_r, mm_l, mm_r, rhoc_12

real*8 :: intflux(g_neqns), intpflux(2)

real*8, intent(in) :: ul(g_neqns), ur(g_neqns), vl(g_neqns), vr(g_neqns)

        !----- Left state
        arho1_l  = ul(1)
        arhou1_l = ul(2)
        arho2_l  = ul(3)
        arhou2_l = ul(4)

        alpha1_l = vl(1)
        p_l      = vl(2)
        u1_l     = vl(3)
        u2_l     = vl(4)
        alpha2_l = 1.0 - alpha1_l

        rho1_l = arho1_l / alpha1_l
        rho2_l = arho2_l / alpha2_l

        a1_l = agas
        a2_l = eos2_ssound(rho2_l)

        !----- Right state
        arho1_r  = ur(1)
        arhou1_r = ur(2)
        arho2_r  = ur(3)
        arhou2_r = ur(4)

        alpha1_r = vr(1)
        p_r      = vr(2)
        u1_r     = vr(3)
        u2_r     = vr(4)
        alpha2_r = 1.0 - alpha1_r

        rho1_r = arho1_r / alpha1_r
        rho2_r = arho2_r / alpha2_r

        a1_r = agas
        a2_r = eos2_ssound(rho2_r)

        ! interface speed of sound
        a1_12 = 0.5 * (a1_l+a1_r)
        a2_12 = 0.5 * (a2_l+a2_r)

        deltap = p_l - p_r
        rho1_12 = 0.5 * (rho1_l+rho1_r)
        rho2_12 = 0.5 * (rho2_l+rho2_r)

        ! common speed of sound choice
        ac_12 = 0.5 * (a1_12+a2_12)
        vrel = max(dabs(u1_l-u2_l), dabs(u1_r-u2_r))

        ! mixture velocities
        vmix_l = (arhou1_l + arhou2_l)/(arho1_l+arho2_l)
        vmix_r = (arhou1_r + arhou2_r)/(arho1_r+arho2_r)
        rhoc_12 = 0.5 * (arho1_r+arho2_r + arho1_l+arho2_l)

        ! mach numbers
        m1_l = u1_l/ac_12
        m1_r = u1_r/ac_12
        m2_l = u2_l/ac_12
        m2_r = u2_r/ac_12

        mm_l = vmix_l/ac_12
        mm_r = vmix_r/ac_12

        !--- phase 1
        !--------------------------------------------------------

        ! split mach numbers
        if (dabs(m1_l) .lt. 1.0) then
           mplus = 0.25 * (m1_l + 1.0) * (m1_l + 1.0)
        else
           macht = 0.5 * (m1_l + dabs(m1_l))
           mplus = macht
        end if

        if (dabs(m1_r) .lt. 1.0) then
           mminu = - 0.25 * (m1_r - 1.0) * (m1_r - 1.0)
        else
           macht = 0.5 * (m1_r - dabs(m1_r))
           mminu = macht
        end if

        ! correction
        if ((dabs(m1_l)<1.0) .and. (dabs(m1_r)<1.0)) then
           mcorr_12 = 0.25 * ( ( dsqrt(0.5*(m1_l*m1_l + m1_r*m1_r)) - 1.0 )**2.0 )
        else
           mcorr_12 = 0.0
        end if

        mplus_12 = mcorr_12 * (1.0 - (deltap + dabs(deltap)) / (2.0*rho1_l*ac_12*ac_12))
        mminu_12 = mcorr_12 * (1.0 + (deltap - dabs(deltap)) / (2.0*rho1_r*ac_12*ac_12))

        disp = vrel * max(alpha1_l,alpha1_r) * (arho1_r-arho1_l) * 10.0

        ! mass fluxes
        mdotplus = arho1_l * ac_12 * (mplus - mplus_12) !- disp
        mdotminu = arho1_r * ac_12 * (mminu + mminu_12) !- disp

        ! convection flux
        intflux(1) = mdotplus        + mdotminu
        intflux(2) = mdotplus * u1_l + mdotminu * u1_r

        !--- phase 2
        !--------------------------------------------------------

        ! split mach numbers
        if (dabs(m2_l) .lt. 1.0) then
           mplus = 0.25 * (m2_l + 1.0) * (m2_l + 1.0)
        else
           macht = 0.5 * (m2_l + dabs(m2_l))
           mplus = macht
        end if

        if (dabs(m2_r) .lt. 1.0) then
           mminu = - 0.25 * (m2_r - 1.0) * (m2_r - 1.0)
        else
           macht = 0.5 * (m2_r - dabs(m2_r))
           mminu = macht
        end if

        ! correction
        if ((dabs(m2_l)<1.0) .and. (dabs(m2_r)<1.0)) then
           mcorr_12 = 0.25 * ( ( dsqrt(0.5*(m2_l*m2_l + m2_r*m2_r)) - 1.0 )**2.0 )
        else
           mcorr_12 = 0.0
        end if

        mplus_12 = mcorr_12 * (1.0 - (deltap + dabs(deltap)) / (2.0*rho2_l*ac_12*ac_12))
        mminu_12 = mcorr_12 * (1.0 + (deltap - dabs(deltap)) / (2.0*rho2_r*ac_12*ac_12))

        disp = vrel * max(alpha2_l,alpha2_r) * (arho2_r-arho2_l) * 10.0

        ! mass fluxes
        mdotplus = arho2_l * ac_12 * (mplus - mplus_12) !- disp
        mdotminu = arho2_r * ac_12 * (mminu + mminu_12) !- disp

        ! convection flux
        intflux(3) = mdotplus        + mdotminu
        intflux(4) = mdotplus * u2_l + mdotminu * u2_r

        !--- pressure flux
        !--------------------------------------------------------

        ! split mach numbers
        if (dabs(mm_l) .lt. 1.0) then
           pplus = 0.5  * (1.0 + mm_l)
        else
           macht = 0.5 * (mm_l + dabs(mm_l))
           pplus = macht / mm_l
        end if

        if (dabs(mm_r) .lt. 1.0) then
           pminu = 0.5  * (1.0 - mm_r)
        else
           macht = 0.5 * (mm_r - dabs(mm_r))
           pminu = macht / mm_r
        end if

        ! interface pressure
        p_12 =  0.5 * (p_l+p_r) &
              + 0.5 * (deltap) * (pplus-pminu) &
              + rhoc_12*ac_12*ac_12 * (pplus+pminu-1.0)

        ! pressure flux
        intpflux(1) = p_12
        intpflux(2) = p_12

end subroutine LDFSS_psreal

!----------------------------------------------------------------------------------------------
!----- Boundary conditions:
!----------------------------------------------------------------------------------------------

subroutine get_bc_4eq(ucons)

real*8  :: rho_g, u_g, pr_g, ucons(g_neqns,0:imax+1)

        !--- extrapolation
        if ((iprob .eq. 0) .or. (iprob .eq. 2)) then
           ucons(1,0) = alpha1_fs       * rho1_fs
           ucons(2,0) = alpha1_fs       * rho1_fs * u1_fs
           ucons(3,0) = (1.0-alpha1_fs) * rho2_fs
           ucons(4,0) = (1.0-alpha1_fs) * rho2_fs * u2_fs
           ucons(:,imax+1) = ucons(:,imax)
        else if (iprob .eq. 1) then
           ucons(:,0) = ucons(:,1)
           ucons(:,imax+1) = ucons(:,imax)
        else
           write(*,*) "BC's not set for iprob = ", iprob
        end if

end subroutine get_bc_4eq

!----------------------------------------------------------------------------------------------
!----- Primitive variable reconstruction (MC):
!----------------------------------------------------------------------------------------------

subroutine reconstruct_uprim(uprim)

integer :: ie, ieqn
real*8  :: vol, fwddiff, bwddiff, cntdiff, &
           uprim(ndof,g_neqns,0:imax+1)

real*8  :: umin, umax, ui, ug, uplus, psil, psir, psi(imax)

        do ie = 1,imax

        vol = coord(ie+1)-coord(ie)

        do ieqn = 1,g_neqns

           fwddiff = ( uprim(1,ieqn,ie+1) - uprim(1,ieqn,ie) ) / vol
           bwddiff = ( uprim(1,ieqn,ie) - uprim(1,ieqn,ie-1) ) / vol
           cntdiff = (uprim(1,ieqn,ie+1) - uprim(1,ieqn,ie-1)) / (2.0*vol)

           uprim(2,ieqn,ie) = cntdiff

        end do !ieqn
        end do !ie

        !--- hardset boundary conditions
        uprim(2,:,0)      = 0.0
        uprim(2,:,imax+1) = 0.0

        !--- barth limiter

        do ieqn = 1,g_neqns

        do ie = 1,imax

           vol = coord(ie+1) - coord(ie)

           umin = min(uprim(1,ieqn,ie-1),uprim(1,ieqn,ie),uprim(1,ieqn,ie+1))
           umax = max(uprim(1,ieqn,ie-1),uprim(1,ieqn,ie),uprim(1,ieqn,ie+1))

           ui = uprim(1,ieqn,ie)

           ! left face
           ug = - 0.5*vol*uprim(2,ieqn,ie)

           if (ug > 0.0) then
              uplus = umax-ui
           else
              uplus = umin-ui
           end if

           if (dabs(ug) > 1.d-12) then
                   psil = min(1.0,uplus/ug)
           else
                   psil = 1.0
           end if

           ! right face
           ug = + 0.5*vol*uprim(2,ieqn,ie)

           if (ug > 0.0) then
              uplus = umax-ui
           else
              uplus = umin-ui
           end if

           if (dabs(ug) > 1.d-12) then
                   psir = min(1.0,uplus/ug)
           else
                   psir = 1.0
           end if

           psi(ie) = min(psil,psir)

        end do !ie

        do ie = 1,imax
           uprim(2,ieqn,ie) = psi(ie) * uprim(2,ieqn,ie)
        end do

        end do

end subroutine reconstruct_uprim

!----------------------------------------------------------------------------------------------

END MODULE rhs_flux
