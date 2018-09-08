!!---------------------------------------------------------------------------------------
!!----- Non-equilibrium Isothermal Multiphase (EOS Module)
!!----- by
!!----- Aditya K Pandare
!!----- Department of Mechanical and Aerospace Engineering,
!!----- North Carolina State University.
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
!----- Decode primitive variables from conserved 
!----------------------------------------------------------------------------------------------

subroutine decode_uprim(ucons, uprim)

real*8, intent(in) :: ucons(neqns,0:imax+1)

integer :: it, itNewt, ie
real*8  :: pres, pres0, coefa, coefb, taitrel, funcp, dfuncp, rho1
real*8  :: uprim(ndof,neqns,0:imax+1)

        do ie = 0,imax+1

        itNewt = 100
        pres0 = uprim(1,2,ie)
        pres  = pres0

        !--- Newton loop to find correct pressure
        do it = 1,itNewt

           coefa = ucons(1,ie) * rgas * tinf
           coefb = ucons(3,ie) / rho0

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

        uprim(1,1,ie) = ucons(1,ie)/rho1
        uprim(1,2,ie) = pres
        uprim(1,3,ie) = ucons(2,ie)/ucons(1,ie)
        uprim(1,4,ie) = ucons(4,ie)/ucons(3,ie)

        end do !ie

end subroutine decode_uprim

!----------------------------------------------------------------------------------------------

END MODULE eos
