!!---------------------------------------------------------------------------------------
!!----- Multi-material and multi-phase (Time-stepping module)
!!----- by
!!----- Aditya K Pandare
!!---------------------------------------------------------------------------------------

MODULE time_integration

USE preprocess

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Explicit TVD-RK3 time-stepping:
!----------------------------------------------------------------------------------------------

subroutine ExplicitRK3_4eq(ucons, uconsn, uprim, uprimn)

integer :: itstep, ielem, ieqn, istage
real*8  :: vol,time
real*8  :: uprim(ndof,g_neqns,0:imax+1),uprimn(ndof,g_neqns,0:imax+1), &
           ucons(g_neqns,0:imax+1),uconsn(g_neqns,0:imax+1),uconsi(g_neqns,0:imax+1), &
           k1(3),k2(3)
real*8  :: rhsel(g_neqns,imax)

        time = 0.d0

        k1(1) = 0.0     !0.0
        k1(2) = 3.0/4.0 !0.5
        k1(3) = 1.0/3.0

        k2(1) = 1.0     !1.0
        k2(2) = 1.0/4.0 !0.5
        k2(3) = 2.0/3.0

        do itstep = 1,ntstep

           uconsi = uconsn

           !---------------------------------------------------------
           !--- RK stages
           do istage = 1,3 !2

              rhsel(:,:) = 0.d0

              if (nsdiscr .eq. 0) then
                 call flux_p0(uprim,uconsi,rhsel)
                 call get_pgradalpha(uprim,rhsel)
              else if (nsdiscr .eq. 1) then
                 call flux_p0p1(uprim,uconsi,rhsel)
                 call get_pgradalpha_p0p1(uprim,rhsel)
              end if

              do ielem = 1,imax

              vol = coord(ielem+1)-coord(ielem)

              ! gravity
              if (iprob.eq.2) then
                 rhsel(2,ielem) = rhsel(2,ielem) + uconsn(1,ielem)*9.81*vol
                 rhsel(4,ielem) = rhsel(4,ielem) + uconsn(3,ielem)*9.81*vol
              end if

              do ieqn  = 1,g_neqns

                   ucons(ieqn,ielem) =   k1(istage) *   uconsn(ieqn,ielem) &
                                       + k2(istage) * ( uconsi(ieqn,ielem) &
                                                      + dt * rhsel(ieqn,ielem)/ vol)

              end do !ieqn
              end do !ielem

              call get_bc_4eq(ucons)
              call decode_uprim(ucons,uprim)
              !call limit_disphase(ucons, uprim, uprimn)
              call blend_disphase(ucons, uprim)

              uconsi = ucons
              uprimn = uprim

           end do !istage
           !---------------------------------------------------------

           time = time + dt

           !----- Screen-output:
           if ((itstep.eq.1) .or. (mod(itstep,n_screen).eq.0)) then
           write(*,*) "--------------------------------------------"
           write(*,*) "  itstep: ", itstep, "   Time: ", time
           write(*,*) "  Time step: ", dt
           write(*,*) "--------------------------------------------"
           write(*,*) " "
           end if

           !--- solution update
           uconsn(:,:)   = ucons(:,:)

           !----- File-output:
           if ((mod(itstep,n_opfile).eq.0).or.(itstep.eq.1)) then
           call gnuplot_flow_4eq(uprim, ucons, itstep)
           end if

        end do !itstep

        itstep = itstep - 1

        call gnuplot_flow_4eq(uprim, ucons, itstep)

        !----- Screen-output:
        write(*,*) "-----------------------------------------------"
        write(*,*) "-------- FINAL OUTPUT: ----- TVD-RK3: ---------"
        write(*,*) "  itstep: ", itstep, "   Time: ", time
        write(*,*) "  Time step: ", dt
        write(*,*) "-----------------------------------------------"
        write(*,*) "-----------------------------------------------"
        write(*,*) " "

end subroutine ExplicitRK3_4eq

!----------------------------------------------------------------------------------------------

subroutine ExplicitRK3_mm6eq(ucons, uconsn, uprim, uprimn)

integer :: itstep, ielem, ieqn, istage
real*8  :: vol,time
real*8  :: uprim(ndof,g_neqns,0:imax+1),uprimn(ndof,g_neqns,0:imax+1), &
           ucons(g_neqns,0:imax+1),uconsn(g_neqns,0:imax+1),uconsi(g_neqns,0:imax+1), &
           k1(3),k2(3)
real*8  :: rhsel(g_neqns,imax)

  time = 0.d0

  k1(1) = 0.0     !0.0
  k1(2) = 3.0/4.0 !0.5
  k1(3) = 1.0/3.0

  k2(1) = 1.0     !1.0
  k2(2) = 1.0/4.0 !0.5
  k2(3) = 2.0/3.0

  do itstep = 1,ntstep

     uconsi = uconsn

     !---------------------------------------------------------
     !--- RK stages
     do istage = 1,3 !2

        rhsel(:,:) = 0.d0

        if (nsdiscr .eq. 0) then
           call flux_p0_mm6eq(uprim,uconsi,rhsel)
           !call get_pugradalpha(uconsi,rhsel)
        else if (nsdiscr .eq. 1) then
           write(*,*) "p0p1 not set up for mm6eq!" 
           exit
        end if

        do ielem = 1,imax

        vol = coord(ielem+1)-coord(ielem)

        do ieqn  = 1,g_neqns

             ucons(ieqn,ielem) =   k1(istage) *   uconsn(ieqn,ielem) &
                                 + k2(istage) * ( uconsi(ieqn,ielem) &
                                                + dt * rhsel(ieqn,ielem)/ vol)

        end do !ieqn
        end do !ielem

        call get_bc_mm6eq(ucons)
        !call decode_uprim(ucons,uprim)
        !call blend_disphase(ucons, uprim)

        uconsi = ucons
        uprimn = uprim

     end do !istage
     !---------------------------------------------------------

     time = time + (dt/a_nd)

     !----- Screen-output:
     if ((itstep.eq.1) .or. (mod(itstep,n_screen).eq.0)) then
     write(*,*) "--------------------------------------------"
     write(*,*) "  itstep: ", itstep, "   Time: ", time
     write(*,*) "  Time step: ", dt/a_nd
     write(*,*) "--------------------------------------------"
     write(*,*) " "
     end if

     !--- solution update
     uconsn(:,:)   = ucons(:,:)

     !----- File-output:
     if ((mod(itstep,n_opfile).eq.0).or.(itstep.eq.1)) then
     call gnuplot_flow_mm6eq(ucons, itstep)
     end if

  end do !itstep

  itstep = itstep - 1

  call gnuplot_flow_mm6eq(ucons, itstep)

  !----- Screen-output:
  write(*,*) "-----------------------------------------------"
  write(*,*) "-------- FINAL OUTPUT: ----- TVD-RK3: ---------"
  write(*,*) "  itstep: ", itstep, "   Time: ", time
  write(*,*) "  Time step: ", dt/a_nd
  write(*,*) "-----------------------------------------------"
  write(*,*) "-----------------------------------------------"
  write(*,*) " "

end subroutine ExplicitRK3_mm6eq

!----------------------------------------------------------------------------------------------
!----- Disappearing phase treatment:
!----- using a hard-bound
!----------------------------------------------------------------------------------------------

subroutine limit_disphase(ucons, uprim, uprimn)

real*8, intent(in) :: uprimn(ndof,g_neqns,0:imax+1)

integer :: ie
real*8  :: alpha1, u1, u2, &
           ucons(g_neqns,0:imax+1), uprim(ndof,g_neqns,0:imax+1)

        do ie = 1,imax
           alpha1 = uprim(1,1,ie)
           u1     = uprim(1,3,ie)
           u2     = uprim(1,4,ie)
        
           !--- phase-1 disappearing
           if (alpha1 .le. 10.0 * alphamin) then
              uprim(1,3,ie) = uprimn(1,3,ie)
              ucons(2,ie)   = ucons(1,ie) * uprim(1,3,ie)

           !--- phase-2 disappearing
           elseif (alpha1 .ge. 1.0-(10.0*alphamin)) then
              uprim(1,4,ie) = uprimn(1,4,ie)
              ucons(4,ie)   = ucons(3,ie) * uprim(1,4,ie)

           end if

        end do !ie

end subroutine limit_disphase

!----------------------------------------------------------------------------------------------
!----- Disappearing phase treatment:
!----- using blending
!----------------------------------------------------------------------------------------------

subroutine blend_disphase(ucons, uprim)

integer :: ie
real*8  :: alpha1, alpha2, u1, u2, pres, xi, gxi, rho, epsmin, epsmax, &
           ucons(g_neqns,0:imax+1), uprim(ndof,g_neqns,0:imax+1)

        epsmin = 1.0  * alphamin
        epsmax = 10.0 * alphamin

        do ie = 1,imax
           alpha1 = uprim(1,1,ie)
           pres   = uprim(1,2,ie)
           u1     = uprim(1,3,ie)
           u2     = uprim(1,4,ie)
           alpha2 = 1.0 -alpha1
        
           !--- phase-1 disappearing
           if (alpha1 .lt. epsmax) then

              if ((alpha1 .ge. epsmin)) then
                    xi = (alpha1 - epsmin)/(epsmax - epsmin)
              elseif (alpha1 .lt. epsmin) then
                    alpha1 = epsmin
                    xi = 0.0
              end if

              Gxi = - xi*xi * ((2.0 * xi) - 3.0)
              u1 = (Gxi*u1) + ((1.0 - Gxi)*u2)
              rho = eos1_density(pres)

              uprim(1,3,ie) = u1

              ucons(1,ie) = alpha1 * rho
              ucons(2,ie) = ucons(1,ie) * uprim(1,3,ie)

           !--- phase-2 disappearing
           elseif (alpha2 .lt. epsmax) then

              if ((alpha2 .ge. epsmin)) then
                    xi = (alpha2 - epsmin)/(epsmax - epsmin)
              elseif (alpha2 .lt. epsmin) then
                    alpha2 = epsmin
                    xi = 0.0
              end if

              Gxi = - xi*xi * ((2.0 * xi) - 3.0)
              u2 = (Gxi*u2) + ((1.0 - Gxi)*u1)
              rho = eos2_density(pres)

              uprim(1,4,ie) = u2

              ucons(3,ie) = alpha2 * rho
              ucons(4,ie) = ucons(3,ie) * uprim(1,4,ie)

           end if

        end do !ie

end subroutine blend_disphase

!----------------------------------------------------------------------------------------------

END MODULE time_integration
