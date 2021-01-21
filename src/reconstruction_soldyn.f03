!!------------------------------------------------------------------------------
!!----- Hyperelastic solid dynamics (Reconstruction-limiting module)
!!----- by
!!----- Aditya K Pandare
!!------------------------------------------------------------------------------

MODULE reconstruction_soldyn

USE glob_var
USE eos
USE reconstruction_mm6eq

implicit none

CONTAINS

!-------------------------------------------------------------------------------
!----- P0P1 reconstruction:
!-------------------------------------------------------------------------------

subroutine reconst_p0p1_soldyn(ucons, uprim, ndof_el)

integer, intent(in) :: ndof_el(2,0:imax+1)

integer :: ie
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  !--- central difference reconstruction (least-squares for uniform meshes)
  do ie = 1,imax
    ucons(2,:,ie) = 0.25 * (ucons(1,:,ie+1) - ucons(1,:,ie-1))
  end do !ie

  !--- limit second-order solution
  call limiting_p1_soldyn(ucons, uprim)

end subroutine reconst_p0p1_soldyn

!-------------------------------------------------------------------------------
!----- P1 limiting:
!-------------------------------------------------------------------------------

subroutine limiting_p1_soldyn(ucons, uprim)

real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)

  select case (g_nlim)

  case(0)

  case(2)
    call vertexbased_p1_soldyn(ucons, uprim)

  case default
    write(*,*) "Error: incorrect p1-limiter index in control file: ", g_nlim
    call exit

  end select

end subroutine limiting_p1_soldyn

!-------------------------------------------------------------------------------
!----- vertexbased limiter for P1 dofs:
!-------------------------------------------------------------------------------

subroutine vertexbased_p1_soldyn(ucons, uprim)

integer :: ie, ieqn
real*8  :: theta(g_neqns)
real*8  :: uneigh(2,g_neqns,-1:1), ucons(g_tdof,g_neqns,0:imax+1), &
           uprim(g_tdof,g_nprim,0:imax+1)

  do ie = 1,imax

    !--- 1. obtain limiter function for individual unknowns
    uneigh(1:2,:,-1) = ucons(1:2,:,ie-1)
    uneigh(1:2,:,0)  = ucons(1:2,:,ie)
    uneigh(1:2,:,1)  = ucons(1:2,:,ie+1)

    call vertexbased_fn(g_neqns, 2.0, 1.0, uneigh, theta)

    !--- 3b. Obtain limiter functions for equation system in single-material cell

    do ieqn = 1,g_neqns
      ucons(2,ieqn,ie) = theta(ieqn) * ucons(2,ieqn,ie)
    end do !ieqn

  end do !ie

end subroutine vertexbased_p1_soldyn

!-------------------------------------------------------------------------------

END MODULE reconstruction_soldyn
