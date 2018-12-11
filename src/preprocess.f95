!!---------------------------------------------------------------------------------------
!!----- Multi-material and multi-phase (Preprocessing module)
!!----- by
!!----- Aditya K Pandare
!!---------------------------------------------------------------------------------------

MODULE preprocess

USE rhs_flux
USE rhs_flux_mm6eq

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Read the control file:
!----------------------------------------------------------------------------------------------

subroutine read_cntl()

integer :: i

        open(12, file = 'setflow.cntl')

        read(12,*) i_system
        read(12,*) ! blank line
        read(12,*) iprob
        read(12,*) imax
        read(12,*) g_lbflag, g_rbflag
        read(12,*) i_restart
        read(12,*) i_flux
        read(12,*) ! blank line
        read(12,*) g_nprelx, g_prelct
        read(12,*) g_gam1, g_pc1, g_cp1
        read(12,*) g_gam2, g_pc2, g_cp2
        read(12,*) ! blank line
        read(12,*) u_fs
        read(12,*) pr1_fs, pr2_fs
        read(12,*) ! blank line
        read(12,*) g_nsdiscr, g_nlim
        read(12,*) dt_u
        read(12,*) ntstep
        read(12,*) ! blank line
        read(12,*) n_opfile
        read(12,*) n_screen
        read(12,*) ! blank line

        close(12)

end subroutine read_cntl

!----------------------------------------------------------------------------------------------
!----- Mesh and CS-area:
!----------------------------------------------------------------------------------------------

subroutine gen_mesh()

integer :: ipoin
real*8  :: ldomn, dx

        if ((iprob.eq.2) .and. (i_system.eq.0)) then
                ldomn = 12.0
        else
                ldomn = 1.0
        end if

        dx = ldomn/imax
        coord(0)  = -dx
        coord(1)  = 0.d0

        do ipoin = 2,imax+1
           coord(ipoin) = coord(ipoin-1) + dx
        end do !ipoin

        coord(imax+2) = coord(imax+1) + dx

end subroutine gen_mesh

!----------------------------------------------------------------------------------------------
!----- Reference quantities and nondimensionalization:
!----------------------------------------------------------------------------------------------

subroutine nondimen_mm6eq()

  a_nd = dsqrt(pr1_fs/rho1_fs)
  t_nd = t1_fs
  p_nd = pr1_fs
  rho_nd = rho1_fs

  write(*,*)" Reference quantities used in non-dimensionalization:"
  write(*,*)"  Speed of sound: ", a_nd
  write(*,*)"  Temperature:    ", t_nd
  write(*,*)"  Pressure:       ", p_nd
  write(*,*)"  Density:        ", rho_nd

  !--- nondimensionalization
  rho1_fs = rho1_fs/rho_nd
  rho2_fs = rho2_fs/rho_nd
  pr1_fs = pr1_fs/p_nd
  pr2_fs = pr2_fs/p_nd
  t1_fs = t1_fs/t_nd
  t2_fs = t2_fs/t_nd
  u_fs = u_fs/a_nd
  dt_u = dt_u*a_nd

  g_pc1 = g_pc1/p_nd
  g_pc2 = g_pc2/p_nd
  g_cp1 = g_cp1/ (a_nd*a_nd/t_nd)
  g_cp2 = g_cp2/ (a_nd*a_nd/t_nd)

end subroutine nondimen_mm6eq

!----------------------------------------------------------------------------------------------
!----- Solution initialization for 4-eq isothermal 2fluid model:
!----------------------------------------------------------------------------------------------

subroutine init_soln_4eq(uprim, uprimn, ucons, uconsn)

integer :: i, ielem
real*8  :: uprim(g_tdof,g_neqns,0:imax+1), uprimn(g_tdof,g_neqns,0:imax+1), &
           ucons(g_tdof,g_neqns,0:imax+1), uconsn(g_tdof,g_neqns,0:imax+1)
real*8  :: xf, pl, pr

        !--- SCD
        if (iprob .eq. 0) then

           g_alphamin = 1.d-6

           alpha1_fs = g_alphamin
           u1_fs = 10.0
           u2_fs = 10.0
           pr_fs = 1.1d5
           rho1_fs = eos1_density(pr_fs)
           rho2_fs = eos2_density(pr_fs)

           do ielem = 0,imax+1

              xf = coord(ielem)
              pl = 1.1d5
              pr = 1.1d5

              if (xf .le. 0.5) then
                 ucons(1,1,ielem) = g_alphamin * eos1_density(pl)
                 ucons(1,2,ielem) = ucons(1,1,ielem) * u1_fs
                 ucons(1,3,ielem) = (1.0-g_alphamin) * eos2_density(pl)
                 ucons(1,4,ielem) = ucons(1,3,ielem) * u2_fs
                 uprim(1,2,ielem) = pl
              else
                 ucons(1,1,ielem) = (1.0-g_alphamin) * eos1_density(pr)
                 ucons(1,2,ielem) = ucons(1,1,ielem) * u1_fs
                 ucons(1,3,ielem) = g_alphamin * eos2_density(pr)
                 ucons(1,4,ielem) = ucons(1,3,ielem) * u2_fs
                 uprim(1,2,ielem) = pr
              end if

           end do !ielem

        !--- Saurel 2001 water-air shocktube
        elseif (iprob .eq. 1) then

           g_alphamin = 1.d-6

           do ielem = 0,imax+1

              xf = coord(ielem)
              pl = 1.d9
              pr = 101325.0

              if (xf .le. 0.7) then
                 ucons(1,1,ielem) = g_alphamin * eos1_density(pl)
                 ucons(1,2,ielem) = 0.0
                 ucons(1,3,ielem) = (1.0-g_alphamin) * eos2_density(pl)
                 ucons(1,4,ielem) = 0.0
                 uprim(1,2,ielem) = pl
              else
                 ucons(1,1,ielem) = (1.0-g_alphamin) * eos1_density(pr)
                 ucons(1,2,ielem) = 0.0
                 ucons(1,3,ielem) = g_alphamin * eos2_density(pr)
                 ucons(1,4,ielem) = 0.0
                 uprim(1,2,ielem) = pr
              end if

           end do !ielem

        !--- Ransom's faucet   
        else if (iprob .eq. 2) then

           g_alphamin = 1.d-6

           alpha1_fs = 0.2
           u1_fs = 0.0
           u2_fs = 10.0
           pr_fs = 1.1d5
           rho1_fs = eos1_density(pr_fs)
           rho2_fs = eos2_density(pr_fs)

           do ielem = 0,imax+1

              xf = coord(ielem)
              pl = 1.1d5
              pr = 1.1d5

              ucons(1,1,ielem) = alpha1_fs * eos1_density(pl)
              ucons(1,2,ielem) = ucons(1,1,ielem) * u1_fs
              ucons(1,3,ielem) = (1.0-alpha1_fs) * eos2_density(pl)
              ucons(1,4,ielem) = ucons(1,3,ielem) * u2_fs
              uprim(1,2,ielem) = pl

           end do !ielem

        else
           write(*,*) "Incorrect problem setup code!"

        end if

        ! boundary conditions:
        call get_bc_4eq(ucons)
        call decode_uprim(ucons,uprim)

        uprimn(:,:,:) = uprim(:,:,:)
        uconsn(:,:,:) = ucons(:,:,:)

        call gnuplot_flow_4eq(uprim, ucons, 0)

end subroutine init_soln_4eq

!----------------------------------------------------------------------------------------------
!----- Solution initialization for multimaterial 6-eq pressure non-equilibrium 
!----- velocity equilibrium 2fluid model:
!----------------------------------------------------------------------------------------------

subroutine init_soln_mm6eq(ucons, uconsn)

integer :: i, ielem
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uconsn(g_tdof,g_neqns,0:imax+1)
real*8  :: s(g_neqns), xf, p1l, p1r, t1l, t1r, &
           ul, ur, p2l, p2r, t2l, t2r, rho1, rho2

  !--- Gaussian
  !----------
  if (iprob .eq. -1) then

     g_alphamin = 1.d-10

     alpha1_fs = g_alphamin
     t1_fs = 300.0
     t2_fs = 300.0
     rho1_fs = eos3_density(g_gam1, g_cp1, g_pc1, pr1_fs, t1_fs)
     rho2_fs = eos3_density(g_gam2, g_cp2, g_pc2, pr2_fs, t2_fs)

     call nondimen_mm6eq()

     if (g_nsdiscr .ge. 1) then
       call weakinit_p1(ucons)

     else
       do ielem = 1,imax
         xf = 0.5*(coord(ielem) + coord(ielem+1))

         s = gaussian(xf,0.0)
         ucons(1,:,ielem) = s(:)

         if (g_nsdiscr .ge. 1) then
           ucons(2,:,ielem) = 0.0
         end if
       end do !ielem

     end if

  !--- SCD
  !----------
  else if (iprob .eq. 0) then

     g_alphamin = 1.d-10

     alpha1_fs = g_alphamin
     t1_fs = 300.0
     t2_fs = 300.0
     rho1_fs = eos3_density(g_gam1, g_cp1, g_pc1, pr1_fs, t1_fs)
     rho2_fs = eos3_density(g_gam2, g_cp2, g_pc2, pr2_fs, t2_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr1_fs
     p2l = pr2_fs
     t1l = t1_fs
     t2l = t2_fs
     ul  = u_fs
     ! right state
     p1r = pr1_fs
     p2r = pr2_fs
     t1r = t1_fs
     t2r = t2_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.5) then
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1l, t1l)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2l, t2l)
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = g_alphamin * rho1
           ucons(1,3,ielem) = (1.0-g_alphamin) * rho2
           ucons(1,4,ielem) = (ucons(1,2,ielem)+ucons(1,3,ielem)) * u_fs
           ucons(1,5,ielem) = g_alphamin * eos3_rhoe(g_gam1, g_pc1, p1l, rho1, ul)
           ucons(1,6,ielem) = (1.0-g_alphamin) * eos3_rhoe(g_gam2, g_pc2, p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1r, t1r)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2r, t2r)
           ucons(1,1,ielem) = 1.0-g_alphamin
           ucons(1,2,ielem) = (1.0-g_alphamin) * rho1
           ucons(1,3,ielem) = g_alphamin * rho2
           ucons(1,4,ielem) = (ucons(1,2,ielem)+ucons(1,3,ielem)) * u_fs
           ucons(1,5,ielem) = (1.0-g_alphamin) * eos3_rhoe(g_gam1, g_pc1, p1r, rho1, ur)
           ucons(1,6,ielem) = g_alphamin * eos3_rhoe(g_gam2, g_pc2, p2r, rho2, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
        end if

     end do !ielem

  !--- Sod Shocktube
  !----------
  else if (iprob .eq. 1) then

     g_alphamin = 1.d-10

     alpha1_fs = g_alphamin
     u_fs = 0.0
     pr1_fs = 1.0
     pr2_fs = 1.0
     t1_fs = 3.484321d-3
     t2_fs = 3.484321d-3
     rho1_fs = eos3_density(g_gam1, g_cp1, g_pc1, pr1_fs, t1_fs)
     rho2_fs = eos3_density(g_gam2, g_cp2, g_pc2, pr2_fs, t2_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr1_fs
     p2l = pr2_fs
     t1l = t1_fs
     t2l = t2_fs
     ul  = u_fs
     ! right state
     p1r = 0.1*pr1_fs
     p2r = 0.1*pr2_fs
     t1r = 0.8*t1_fs
     t2r = 0.8*t2_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.5) then
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1l, t1l)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2l, t2l)
           ucons(1,1,ielem) = 1.0-alpha1_fs
           ucons(1,2,ielem) = (1.0-alpha1_fs) * rho1
           ucons(1,3,ielem) = alpha1_fs * rho2
           ucons(1,4,ielem) = (ucons(1,2,ielem)+ucons(1,3,ielem)) * u_fs
           ucons(1,5,ielem) = (1.0-alpha1_fs) * eos3_rhoe(g_gam1, g_pc1, p1l, rho1, ul)
           ucons(1,6,ielem) = alpha1_fs * eos3_rhoe(g_gam2, g_pc2, p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1r, t1r)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2r, t2r)
           ucons(1,1,ielem) = alpha1_fs
           ucons(1,2,ielem) = alpha1_fs * rho1
           ucons(1,3,ielem) = (1.0-alpha1_fs) * rho2
           ucons(1,4,ielem) = (ucons(1,2,ielem)+ucons(1,3,ielem)) * u_fs
           ucons(1,5,ielem) = alpha1_fs * eos3_rhoe(g_gam1, g_pc1, p1r, rho1, ur)
           ucons(1,6,ielem) = (1.0-alpha1_fs) * eos3_rhoe(g_gam2, g_pc2, p2r, rho2, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
        end if

     end do !ielem

  !--- Abgrall's water-air Shocktube
  !----------
  else if (iprob .eq. 2) then

     g_alphamin = 1.d-10

     alpha1_fs = g_alphamin
     u_fs = 0.0
     pr1_fs = 1.0d9 !2.0d8
     pr2_fs = 1.0d9 !2.0d8
     t1_fs = 494.646 !247.323
     t2_fs = 494.646 !247.323
     rho1_fs = eos3_density(g_gam1, g_cp1, g_pc1, pr1_fs, t1_fs)
     rho2_fs = eos3_density(g_gam2, g_cp2, g_pc2, pr2_fs, t2_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr1_fs
     p2l = pr2_fs
     t1l = t1_fs
     t2l = t2_fs
     ul  = u_fs
     ! right state
     p1r = pr1_fs/1.0d4 !pr1_fs/2.0d3
     p2r = pr2_fs/1.0d4 !pr2_fs/2.0d3
     t1r = 0.070441*t1_fs !0.028176*t1_fs
     t2r = 0.070441*t2_fs !0.028176*t2_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.75) then
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1l, t1l)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2l, t2l)
           ucons(1,1,ielem) = 1.0-g_alphamin
           ucons(1,2,ielem) = (1.0-g_alphamin) * rho1
           ucons(1,3,ielem) = g_alphamin * rho2
           ucons(1,4,ielem) = (ucons(1,2,ielem)+ucons(1,3,ielem)) * u_fs
           ucons(1,5,ielem) = (1.0-g_alphamin) * eos3_rhoe(g_gam1, g_pc1, p1l, rho1, ul)
           ucons(1,6,ielem) = g_alphamin * eos3_rhoe(g_gam2, g_pc2, p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1r, t1r)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2r, t2r)
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = g_alphamin * rho1
           ucons(1,3,ielem) = (1.0-g_alphamin) * rho2
           ucons(1,4,ielem) = (ucons(1,2,ielem)+ucons(1,3,ielem)) * u_fs
           ucons(1,5,ielem) = g_alphamin * eos3_rhoe(g_gam1, g_pc1, p1r, rho1, ur)
           ucons(1,6,ielem) = (1.0-g_alphamin) * eos3_rhoe(g_gam2, g_pc2, p2r, rho2, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
        end if

     end do !ielem

  else
     write(*,*) "Incorrect problem setup code!"

  end if

  ! boundary conditions:
  call get_bc_mm6eq(ucons)

  uconsn(:,:,:) = ucons(:,:,:)

  call gnuplot_flow_mm6eq(ucons, 0)
  call gnuplot_flow_p1_mm6eq(ucons, 0)

end subroutine init_soln_mm6eq

!-------------------------------------------------------------------------------
!----- Weak initialization for DG(P1) for smooth problems:
!-------------------------------------------------------------------------------

subroutine weakinit_p1(ucons)

integer :: ig, ie, ieqn, ngauss
data       ngauss/2/

real*8  :: wi, vol, xc, x, s(g_neqns), rhs(g_tdof,g_neqns), carea(2), weight(2)
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax

    vol = coord(ie+1)-coord(ie)
    xc  = 0.5*(coord(ie+1)+coord(ie))

    rhs(:,:) = 0.0

    do ig = 1,ngauss

      wi = 0.5 * weight(ig) * vol
      x = carea(ig) * 0.5 * vol + xc
      s = gaussian(x,0.0)

      do ieqn = 1,g_neqns
        rhs(1,ieqn) = rhs(1,ieqn) + wi * s(ieqn)
        rhs(2,ieqn) = rhs(2,ieqn) + wi * s(ieqn) * carea(ig)
      end do !ieqn

    end do !ig

    do ieqn = 1,g_neqns
      ucons(1,ieqn,ie) = rhs(1,ieqn) / vol
      ucons(2,ieqn,ie) = rhs(2,ieqn) / (vol/3.0)
    end do !ieqn

  end do !ie

end subroutine weakinit_p1

!-------------------------------------------------------------------------------

function gaussian(x,t)

real*8, intent(in)  :: x, t
real*8  :: xc, al1, rho1, rho2, gaussian(g_neqns)

  xc  = 0.25 + u_fs*t
  al1 = (1.0-alpha1_fs) * dexp( -(x-xc)*(x-xc)/(2.0 * 0.002) ) + alpha1_fs

  rho1 = eos3_density(g_gam1, g_cp1, g_pc1, pr1_fs, t1_fs)
  rho2 = eos3_density(g_gam2, g_cp2, g_pc2, pr2_fs, t2_fs)
  gaussian(1) = al1
  gaussian(2) = al1 * rho1
  gaussian(3) = (1.0-al1) * rho2
  gaussian(4) = (gaussian(2)+gaussian(3)) * u_fs
  gaussian(5) = al1 * eos3_rhoe(g_gam1, g_pc1, pr1_fs, rho1, u_fs)
  gaussian(6) = (1.0-al1) * eos3_rhoe(g_gam2, g_pc2, pr2_fs, rho2, u_fs)

end function

!----------------------------------------------------------------------------------------------
!----- GNUPLOT outputs:
!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow_4eq(uprim, ucons, itstep)

integer, intent(in) :: itstep
real*8,  intent(in) :: uprim(g_tdof,g_neqns,0:imax+1), &
                       ucons(g_tdof,g_neqns,0:imax+1)

integer :: ielem
real*8  :: xcc, pres, rhomix, umix, &
           rho1, rho2, u1, u2, alp1, alp2

character(len=100) :: filename2,filename3

        write(filename2,'(1I50)')itstep
        filename3 = trim(adjustl(filename2)) // '.twofluid.'//'dat'
        open(23,file=trim(adjustl(filename3)),status='unknown')

        do ielem = 1,imax

           xcc = 0.5d0 * (coord(ielem) + coord(ielem+1))
           alp1 = uprim(1,1,ielem)
           pres = uprim(1,2,ielem)
           u1   = uprim(1,3,ielem)
           u2   = uprim(1,4,ielem)
           alp2 = 1.0 - alp1

           rho1 = eos1_density(pres)
           rho2 = eos2_density(pres)

           rhomix = alp1*rho1 + alp2*rho2
           umix   = alp1*  u1 + alp2*  u2

           write(23,'(11E16.6)') xcc, alp1, rhomix, pres, umix, u1, u2, &
                                ucons(1,1,ielem), ucons(1,2,ielem), ucons(1,3,ielem), ucons(1,4,ielem)

        end do !ielem

        close(23)

end subroutine gnuplot_flow_4eq

!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow_mm6eq(ucons, itstep)

integer, intent(in) :: itstep
real*8,  intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: ielem
real*8  :: xcc, pmix, tmix, rhomix, e_mix, &
           arho1, arho2, arhoe1, arhoe2, alp2, trcell
real*8  :: uconsi(g_neqns), uprimi(g_neqns)

character(len=100) :: filename2,filename3

  write(filename2,'(1I50)')itstep
  filename3 = trim(adjustl(filename2)) // '.twofluid.'//'dat'
  open(23,file=trim(adjustl(filename3)),status='unknown')

  write(23,'(12A8)') "# xcc,", &   !1
                     "alp1,", &    !2
                     "rhomix,", &  !3
                     "umix," , &   !4
                     "pmix,", &    !5
                     "tmix,", &    !6
                     "p1,", &      !7
                     "p2,", &      !8
                     "t1,", &      !9
                     "t2,", &      !10
                     "e_m", &      !11
                     "int_cell"    !12

  do ielem = 1,imax

     uconsi = ucons(1,:,ielem)
     call get_uprim_mm6eq(uconsi, uprimi)

     xcc = 0.5d0 * (coord(ielem) + coord(ielem+1))
     alp2 = 1.0 - uprimi(1)
     arho1 = ucons(1,2,ielem)
     arho2 = ucons(1,3,ielem)
     rhomix = ucons(1,2,ielem) + ucons(1,3,ielem)
     arhoe1 = ucons(1,5,ielem)
     arhoe2 = ucons(1,6,ielem)
     pmix = uprimi(1)*uprimi(2) + alp2*uprimi(3)
     tmix = uprimi(1)*uprimi(5) + alp2*uprimi(6)
     e_mix = ((arhoe1 - 0.5*arho1*uprimi(4)*uprimi(4)) + &
              (arhoe2 - 0.5*arho2*uprimi(4)*uprimi(4))) / rhomix

     if ( (uprimi(1) .gt. 10.0*g_alphamin) &
         .and. (uprimi(1) .lt. 1.0-10.0*g_alphamin) ) then
       trcell = 1.0
     else
       trcell = 0.0
     end if

     write(23,'(12E16.6)') xcc, &             !1
                           uprimi(1), &       !2
                           rhomix*rho_nd, &   !3
                           uprimi(4)*a_nd , & !4
                           pmix*p_nd, &       !5
                           tmix*t_nd, &       !6
                           uprimi(2)*p_nd, &  !7
                           uprimi(3)*p_nd, &  !8
                           uprimi(5)*t_nd, &  !9
                           uprimi(6)*t_nd, &  !10
                           e_mix, &           !11
                           trcell             !12

  end do !ielem

  close(23)

end subroutine gnuplot_flow_mm6eq

!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow_p1_mm6eq(ucons, itstep)

integer, intent(in) :: itstep
real*8,  intent(in) :: ucons(g_tdof,g_neqns,0:imax+1)

integer :: ielem
real*8  :: xp, pmix, tmix, rhomix, e_mix, &
           arho1, arho2, arhoe1, arhoe2, alp2, trcell
real*8  :: uconsi(g_neqns), uprimi(g_neqns)

character(len=100) :: filename2,filename3

  write(filename2,'(1I50)')itstep
  filename3 = trim(adjustl(filename2)) // '.dgtwofluid.'//'dat'
  open(24,file=trim(adjustl(filename3)),status='unknown')

  write(24,'(12A8)') "# xcc,", &   !1
                     "alp1,", &    !2
                     "rhomix,", &  !3
                     "umix," , &   !4
                     "pmix,", &    !5
                     "tmix,", &    !6
                     "p1,", &      !7
                     "p2,", &      !8
                     "t1,", &      !9
                     "t2,", &      !10
                     "e_m", &      !11
                     "int_cell"    !12

  do ielem = 1,imax

     ! left face
     if (g_nsdiscr .gt. 0) then
     uconsi = ucons(1,:,ielem) - ucons(2,:,ielem)
     else
     uconsi = ucons(1,:,ielem)
     end if
     call get_uprim_mm6eq(uconsi, uprimi)

     xp = coord(ielem)
     alp2 = 1.0 - uprimi(1)
     arho1 = uconsi(2)
     arho2 = uconsi(3)
     rhomix = uconsi(2) + uconsi(3)
     arhoe1 = uconsi(5)
     arhoe2 = uconsi(6)
     pmix = uprimi(1)*uprimi(2) + alp2*uprimi(3)
     tmix = uprimi(1)*uprimi(5) + alp2*uprimi(6)
     e_mix = ((arhoe1 - 0.5*arho1*uprimi(4)*uprimi(4)) + &
              (arhoe2 - 0.5*arho2*uprimi(4)*uprimi(4))) / rhomix

     if ( (uprimi(1) .gt. 10.0*g_alphamin) &
         .and. (uprimi(1) .lt. 1.0-10.0*g_alphamin) ) then
       trcell = 1.0
     else
       trcell = 0.0
     end if

     write(24,'(12E16.6)') xp, &              !1
                           uprimi(1), &       !2
                           rhomix*rho_nd, &   !3
                           uprimi(4)*a_nd , & !4
                           pmix*p_nd, &       !5
                           tmix*t_nd, &       !6
                           uprimi(2)*p_nd, &  !7
                           uprimi(3)*p_nd, &  !8
                           uprimi(5)*t_nd, &  !9
                           uprimi(6)*t_nd, &  !10
                           e_mix, &           !11
                           trcell             !12

     ! right face
     if (g_nsdiscr .gt. 0) then
     uconsi = ucons(1,:,ielem) + ucons(2,:,ielem)
     else
     uconsi = ucons(1,:,ielem)
     end if
     call get_uprim_mm6eq(uconsi, uprimi)

     xp = coord(ielem+1)
     alp2 = 1.0 - uprimi(1)
     arho1 = uconsi(2)
     arho2 = uconsi(3)
     rhomix = uconsi(2) + uconsi(3)
     arhoe1 = uconsi(5)
     arhoe2 = uconsi(6)
     pmix = uprimi(1)*uprimi(2) + alp2*uprimi(3)
     tmix = uprimi(1)*uprimi(5) + alp2*uprimi(6)
     e_mix = ((arhoe1 - 0.5*arho1*uprimi(4)*uprimi(4)) + &
              (arhoe2 - 0.5*arho2*uprimi(4)*uprimi(4))) / rhomix

     if ( (uprimi(1) .gt. 10.0*g_alphamin) &
         .and. (uprimi(1) .lt. 1.0-10.0*g_alphamin) ) then
       trcell = 1.0
     else
       trcell = 0.0
     end if

     write(24,'(12E16.6)') xp, &              !1
                           uprimi(1), &       !2
                           rhomix*rho_nd, &   !3
                           uprimi(4)*a_nd , & !4
                           pmix*p_nd, &       !5
                           tmix*t_nd, &       !6
                           uprimi(2)*p_nd, &  !7
                           uprimi(3)*p_nd, &  !8
                           uprimi(5)*t_nd, &  !9
                           uprimi(6)*t_nd, &  !10
                           e_mix, &           !11
                           trcell             !12

     write(24,*) " "

  end do !ielem

  close(24)

end subroutine gnuplot_flow_p1_mm6eq

!----------------------------------------------------------------------------------------------

subroutine gnuplot_diagnostics_mm6eq(cons_err, itstep)

integer, intent(in) :: itstep
real*8,  intent(in) :: cons_err(6)

integer :: ielem

  do ielem = 1,imax

    write(33,*) itstep, &       !1
                cons_err(1), &  !2
                cons_err(2), &  !3
                cons_err(3), &  !4
                cons_err(4), &  !5
                cons_err(5), &  !6
                cons_err(6)     !7

  end do !ielem

end subroutine gnuplot_diagnostics_mm6eq

!----------------------------------------------------------------------------------------------

subroutine errorcalc_p1(ucons, t, err_log)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), t

integer :: ig, ie, ieqn, ngauss
data       ngauss/3/

real*8  :: dx, xg, wi, &
           u(g_neqns), up(g_neqns), &
           carea(3), weight(3), &
           s(g_neqns), err, err_log

  call rutope(1, ngauss, carea, weight)

  err_log = 0.0

  if (iprob .eq. -1) then

    do ie = 1,imax
    do ig = 1,ngauss

      dx = coord(ie+1)-coord(ie)
      wi = 0.5 * weight(ig) * dx
      xg = carea(ig) * 0.5*dx + 0.5 * ( coord(ie+1)+coord(ie) )

      ! basis function

      if (g_nsdiscr .ge. 1) then
        do ieqn = 1,g_neqns
          u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie)
        end do !ieqn
      else
        u(:) = ucons(1,:,ie)
      end if

      s = gaussian(xg,t)

      err = u(1) - s(1)

      err_log = err_log + wi*err*err

    end do !ig
    end do !ie

    err_log = dsqrt(err_log)
    err_log = dlog10(err_log)

  else

    write(*,*) "  WARNING: Error computation not configured for iprob = ", iprob

  end if

end subroutine errorcalc_p1

!----------------------------------------------------------------------------------------------

END MODULE preprocess
