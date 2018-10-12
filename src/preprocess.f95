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
        read(12,*) g_gam1, g_pc1, g_cp1
        read(12,*) g_gam2, g_pc2, g_cp2
        read(12,*) ! blank line
        read(12,*) u_fs
        read(12,*) pr1_fs, pr2_fs
        read(12,*) ! blank line
        read(12,*) nsdiscr
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

        if (iprob.eq.2) then
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
real*8  :: uprim(ndof,g_neqns,0:imax+1), uprimn(ndof,g_neqns,0:imax+1), &
           ucons(g_neqns,0:imax+1), uconsn(g_neqns,0:imax+1)
real*8  :: xf, pl, pr

        !--- SCD
        if (iprob .eq. 0) then

           alphamin = 1.d-6

           alpha1_fs = alphamin
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
                 ucons(1,ielem) = alphamin * eos1_density(pl)
                 ucons(2,ielem) = ucons(1,ielem) * u1_fs
                 ucons(3,ielem) = (1.0-alphamin) * eos2_density(pl)
                 ucons(4,ielem) = ucons(3,ielem) * u2_fs
                 uprim(1,2,ielem) = pl
              else
                 ucons(1,ielem) = (1.0-alphamin) * eos1_density(pr)
                 ucons(2,ielem) = ucons(1,ielem) * u1_fs
                 ucons(3,ielem) = alphamin * eos2_density(pr)
                 ucons(4,ielem) = ucons(3,ielem) * u2_fs
                 uprim(1,2,ielem) = pr
              end if

           end do !ielem

        !--- Saurel 2001 water-air shocktube
        elseif (iprob .eq. 1) then

           alphamin = 1.d-6

           do ielem = 0,imax+1

              xf = coord(ielem)
              pl = 1.d9
              pr = 101325.0

              if (xf .le. 0.7) then
                 ucons(1,ielem) = alphamin * eos1_density(pl)
                 ucons(2,ielem) = 0.0
                 ucons(3,ielem) = (1.0-alphamin) * eos2_density(pl)
                 ucons(4,ielem) = 0.0
                 uprim(1,2,ielem) = pl
              else
                 ucons(1,ielem) = (1.0-alphamin) * eos1_density(pr)
                 ucons(2,ielem) = 0.0
                 ucons(3,ielem) = alphamin * eos2_density(pr)
                 ucons(4,ielem) = 0.0
                 uprim(1,2,ielem) = pr
              end if

           end do !ielem

        !--- Ransom's faucet   
        else if (iprob .eq. 2) then

           alphamin = 1.d-6

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

              ucons(1,ielem) = alpha1_fs * eos1_density(pl)
              ucons(2,ielem) = ucons(1,ielem) * u1_fs
              ucons(3,ielem) = (1.0-alpha1_fs) * eos2_density(pl)
              ucons(4,ielem) = ucons(3,ielem) * u2_fs
              uprim(1,2,ielem) = pl

           end do !ielem

        else
           write(*,*) "Incorrect problem setup code!"

        end if

        ! boundary conditions:
        call get_bc_4eq(ucons)
        call decode_uprim(ucons,uprim)

        uprimn(:,:,:) = uprim(:,:,:)
        uconsn(:,:)   = ucons(:,:)

        call gnuplot_flow_4eq(uprim, ucons, 0)

end subroutine init_soln_4eq

!----------------------------------------------------------------------------------------------
!----- Solution initialization for multimaterial 6-eq pressure non-equilibrium 
!----- velocity equilibrium 2fluid model:
!----------------------------------------------------------------------------------------------

subroutine init_soln_mm6eq(uprim, uprimn, ucons, uconsn)

integer :: i, ielem
real*8  :: uprim(ndof,g_neqns,0:imax+1), uprimn(ndof,g_neqns,0:imax+1), &
           ucons(g_neqns,0:imax+1), uconsn(g_neqns,0:imax+1)
real*8  :: xf, p1l, p1r, t1l, t1r, &
           ul, ur, p2l, p2r, t2l, t2r, rho1, rho2

  !--- SCD
  !----------
  if (iprob .eq. 0) then

     alphamin = 1.d-10

     alpha1_fs = alphamin
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
           ucons(1,ielem) = alphamin
           ucons(2,ielem) = alphamin * rho1
           ucons(3,ielem) = (1.0-alphamin) * rho2
           ucons(4,ielem) = (ucons(2,ielem)+ucons(3,ielem)) * u_fs
           ucons(5,ielem) = alphamin * eos3_rhoe(g_gam1, g_pc1, p1l, rho1, ul)
           ucons(6,ielem) = (1.0-alphamin) * eos3_rhoe(g_gam2, g_pc2, p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1r, t1r)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2r, t2r)
           ucons(1,ielem) = 1.0-alphamin
           ucons(2,ielem) = (1.0-alphamin) * rho1
           ucons(3,ielem) = alphamin * rho2
           ucons(4,ielem) = (ucons(2,ielem)+ucons(3,ielem)) * u_fs
           ucons(5,ielem) = (1.0-alphamin) * eos3_rhoe(g_gam1, g_pc1, p1r, rho1, ur)
           ucons(6,ielem) = alphamin * eos3_rhoe(g_gam2, g_pc2, p2r, rho2, ur)
        end if

     end do !ielem

  !--- Shocktube
  !----------
  else if (iprob .eq. 1) then

     alphamin = 1.d-10

     alpha1_fs = alphamin
     u_fs = 0.0
     pr1_fs = 1.d5
     pr2_fs = 1.d5
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
           ucons(1,ielem) = alphamin
           ucons(2,ielem) = alphamin * rho1
           ucons(3,ielem) = (1.0-alphamin) * rho2
           ucons(4,ielem) = (ucons(2,ielem)+ucons(3,ielem)) * u_fs
           ucons(5,ielem) = alphamin * eos3_rhoe(g_gam1, g_pc1, p1l, rho1, ul)
           ucons(6,ielem) = (1.0-alphamin) * eos3_rhoe(g_gam2, g_pc2, p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam1, g_cp1, g_pc1, p1r, t1r)
           rho2 = eos3_density(g_gam2, g_cp2, g_pc2, p2r, t2r)
           ucons(1,ielem) = 1.0-alphamin
           ucons(2,ielem) = (1.0-alphamin) * rho1
           ucons(3,ielem) = alphamin * rho2
           ucons(4,ielem) = (ucons(2,ielem)+ucons(3,ielem)) * u_fs
           ucons(5,ielem) = (1.0-alphamin) * eos3_rhoe(g_gam1, g_pc1, p1r, rho1, ur)
           ucons(6,ielem) = alphamin * eos3_rhoe(g_gam2, g_pc2, p2r, rho2, ur)
        end if

     end do !ielem

  else
     write(*,*) "Incorrect problem setup code!"

  end if

  ! boundary conditions:
  call get_bc_mm6eq(ucons)
  !call decode_uprim(ucons,uprim)

  uprimn(:,:,:) = uprim(:,:,:)
  uconsn(:,:)   = ucons(:,:)

  call gnuplot_flow_mm6eq(ucons, 0)

end subroutine init_soln_mm6eq

!----------------------------------------------------------------------------------------------
!----- GNUPLOT outputs:
!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow_4eq(uprim, ucons, itstep)

integer, intent(in) :: itstep
real*8,  intent(in) :: uprim(ndof,g_neqns,0:imax+1), ucons(g_neqns,0:imax+1)

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
                                ucons(1,ielem), ucons(2,ielem), ucons(3,ielem), ucons(4,ielem)

        end do !ielem

        close(23)

end subroutine gnuplot_flow_4eq

!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow_mm6eq(ucons, itstep)

integer, intent(in) :: itstep
real*8,  intent(in) :: ucons(g_neqns,0:imax+1)

integer :: ielem
real*8  :: xcc, p1, p2, t1, t2, pmix, tmix, rhomix, umix, rhoe_mix, &
           rho1, rho2, rhoe1, rhoe2, alp1, alp2

character(len=100) :: filename2,filename3

  write(filename2,'(1I50)')itstep
  filename3 = trim(adjustl(filename2)) // '.twofluid.'//'dat'
  open(23,file=trim(adjustl(filename3)),status='unknown')

  write(23,'(11A8)') "# xcc,", &   !1
                     "alp1,", &    !2
                     "rhomix,", &  !3
                     "umix," , &   !4
                     "pmix,", &    !5
                     "tmix,", &    !6
                     "p1,", &      !7
                     "p2,", &      !8
                     "t1,", &      !9
                     "t2,", &      !10
                     "rhoe_m"      !11

  do ielem = 1,imax

     xcc = 0.5d0 * (coord(ielem) + coord(ielem+1))
     alp1 = ucons(1,ielem)
     alp2 = 1.0 - alp1
     rho1 = ucons(2,ielem)/alp1
     rho2 = ucons(3,ielem)/alp2
     rhomix = ucons(2,ielem) + ucons(3,ielem) !alp1*rho1 + alp2*rho2
     umix = ucons(4,ielem)/rhomix
     rhoe1 = ucons(5,ielem)/alp1
     rhoe2 = ucons(6,ielem)/alp2
     p1   = eos3_pr(g_gam1, g_pc1, rho1, rhoe1, umix)
     p2   = eos3_pr(g_gam2, g_pc2, rho2, rhoe2, umix)
     t1   = eos3_t(g_gam1, g_cp1, g_pc1, rho1, rhoe1, umix)
     t2   = eos3_t(g_gam2, g_cp2, g_pc2, rho2, rhoe2, umix)
     pmix = alp1*p1 + alp2*p2
     tmix = alp1*t1 + alp2*t2
     rhoe_mix = (alp1 * (rhoe1 - 0.5*rho1*umix*umix) + &
                 alp2 * (rhoe2 - 0.5*rho2*umix*umix)) / rhomix

     write(23,'(11E16.6)') xcc, &            !1
                           alp1, &           !2
                           rhomix*rho_nd, &  !3
                           umix*a_nd , &     !4
                           pmix*p_nd, &      !5
                           tmix*t_nd, &      !6
                           p1*p_nd, &        !7
                           p2*p_nd, &        !8
                           t1*t_nd, &        !9
                           t2*t_nd, &        !10
                           rhoe_mix          !11

  end do !ielem

  close(23)

end subroutine gnuplot_flow_mm6eq

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

END MODULE preprocess
