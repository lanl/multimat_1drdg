!!---------------------------------------------------------------------------------------
!!----- Multi-material and multi-phase (Preprocessing module)
!!----- by
!!----- Aditya K Pandare
!!---------------------------------------------------------------------------------------

MODULE preprocess

USE rhs_flux_mm6eq

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Screen output:
!----------------------------------------------------------------------------------------------

subroutine screen_output()

write(*,*) " "
if (i_system .eq. -1) then
  write(*,*) " k-exactness check for least-squares reconstruction "
else if (i_system .eq. 1) then
  write(*,*) " Multi-material system: "
  write(*,*) " Pressure non-equilibrium " 
  write(*,*) " Velocity equilibrium "
  write(*,*) " "
  write(*,*) "  Indices used: "
  write(*,*) "    nummat:", g_mmi%nummat
  write(*,*) "    iamin: ", g_mmi%iamin
  write(*,*) "    iamax: ", g_mmi%iamax
  write(*,*) "    irmin: ", g_mmi%irmin
  write(*,*) "    irmax: ", g_mmi%irmax
  write(*,*) "    imome: ", g_mmi%imome
  write(*,*) "    iemin: ", g_mmi%iemin
  write(*,*) "    iemax: ", g_mmi%iemax
  write(*,*) " "
  write(*,*) "  # equations:", g_neqns
else
  write(*,*) " System not configured; i_system: ", i_system
end if

if (g_nprelx .eq. 1) then
  write(*,*) " "
  write(*,*) " Pressure relaxation time-scale: ", g_prelct
end if

write(*,*) " "
write(*,*) " nelem = ", imax
write(*,*) " dt =    ", dt_u
write(*,*) " "
write(*,*) " Spatial discretization: "
if (g_nsdiscr .eq. 0) then
  write(*,*) "   DG(P0)"
else if (g_nsdiscr .eq. 1) then
  write(*,*) "   rDG(P0P1)"
else if (g_nsdiscr .eq. 11) then
  write(*,*) "   DG(P1)"
else if (g_nsdiscr .eq. 12) then
  write(*,*) "   rDG(P1P2)"
end if

if (i_system > -1) then
  write(*,*) " "
  write(*,*) " Using appr. Riemann solver: "
  if (i_flux .eq. 1) then
    write(*,*) "   Lax-Friedrichs flux."
  else if (i_flux .eq. 2) then
    write(*,*) "   AUSM+ flux."
  else if (i_flux .eq. 3) then
    write(*,*) "   HLL flux."
  else
    write(*,*) "Invalid flux scheme."
    stop
  end if

  if (g_nsdiscr > 0) then
  write(*,*) " "
  write(*,*) " Material interface detection: "
  if (g_nmatint .eq. 1) then
    write(*,*) "   ON"
  else
    write(*,*) "   OFF"
  end if

  write(*,*) " "
  write(*,*) " 'Redundant' reconst/limiting for a*p, u: "
  if (g_pureco .eq. 1) then
    write(*,*) "   ON"
  else
    write(*,*) "   OFF"
  end if

  write(*,*) " "
  write(*,*) " Limiter: "
  if (g_nlim .eq. 0) then
    write(*,*) "   None."
  else if (g_nlim .eq. 1) then
    write(*,*) "   Min: Vertex-based."
  else if (g_nlim .eq. 2) then
    write(*,*) "   Vertex-based."
  else if (g_nlim .eq. 3) then
    write(*,*) "   Overbee+Vertex-based."
  else if (g_nlim .eq. 4) then
    write(*,*) "   WENO."
  else if (g_nlim .eq. 5) then
    write(*,*) "   Vertex-based+WENO."
  else if (g_nlim .eq. 6) then
    write(*,*) "   LINC+Vertex-based."
  else if (g_nlim .eq. 7) then
    write(*,*) "   LINC+Vertex-based+WENO."
  else if (g_nlim .eq. 8) then
    write(*,*) "   Vertex-based+Positivity."
  else
    write(*,*) "Invalid limiter."
    stop
  end if
  end if
end if

write(*,*) " "

end subroutine screen_output

!----------------------------------------------------------------------------------------------
!----- Read the control file:
!----------------------------------------------------------------------------------------------

subroutine read_cntl()

integer :: imat
logical :: file_exists

        inquire(file='setflow.cntl', exist=file_exists)
        if (.not.file_exists) then
          write(*,*) "Error: Input file 'setflow.cntl' does not exist"
          stop
        end if

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
        read(12,*) g_mmi%nummat

        ! allocate mat-property arrays
        allocate(g_gam(g_mmi%nummat), g_cp(g_mmi%nummat), g_pc(g_mmi%nummat), &
                 alpha_fs(g_mmi%nummat), rhomat_fs(g_mmi%nummat))

        do imat = 1,g_mmi%nummat
          read(12,*) g_gam(imat), g_pc(imat), g_cp(imat)
        end do !imat
        read(12,*) ! blank line
        read(12,*) u_fs
        read(12,*) pr_fs
        read(12,*) ! blank line
        read(12,*) g_nsdiscr
        read(12,*) g_nlim, g_nmatint, g_pureco
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

        if ((iprob==8) .and. (i_system==1)) then
          ldomn = 10.0
          coord(1) = -5.0
        else if ((iprob==9) .and. (i_system==1)) then
          ldomn = 9.0
          coord(1) = 0.0
        else
          ldomn = 1.0
          coord(1)  = 0.d0
        end if

        dx = ldomn/imax
        coord(0)  = coord(1)-dx

        do ipoin = 2,imax+1
           coord(ipoin) = coord(ipoin-1) + dx
        end do !ipoin

        coord(imax+2) = coord(imax+1) + dx

end subroutine gen_mesh

!----------------------------------------------------------------------------------------------
!----- Reference quantities and nondimensionalization:
!----------------------------------------------------------------------------------------------

subroutine nondimen_mm6eq()

integer :: imat

associate (nummat=>g_mmi%nummat)

  a_nd = dsqrt(pr_fs/rhomat_fs(1))
  t_nd = t_fs
  p_nd = pr_fs
  rho_nd = rhomat_fs(1)

  write(*,*)" Reference quantities used in non-dimensionalization:"
  write(*,*)"  Speed of sound: ", a_nd
  write(*,*)"  Temperature:    ", t_nd
  write(*,*)"  Pressure:       ", p_nd
  write(*,*)"  Density:        ", rho_nd
  write(*,*) "-----------------------------------------------"
  write(*,*) "-----------------------------------------------"
  write(*,*)" "

  !--- nondimensionalization
  do imat = 1,nummat
    rhomat_fs(imat) = rhomat_fs(imat)/rho_nd
    g_pc(imat) = g_pc(imat)/p_nd
    g_cp(imat) = g_cp(imat)/ (a_nd*a_nd/t_nd)
  end do !imat
  pr_fs = pr_fs/p_nd
  t_fs = t_fs/t_nd
  u_fs = u_fs/a_nd
  dt_u = dt_u*a_nd

end associate

end subroutine nondimen_mm6eq

!----------------------------------------------------------------------------------------------
!----- Solution initialization for k-exactness check:
!----------------------------------------------------------------------------------------------

subroutine init_soln_kex(ucons, ndof_el)

integer :: ig, ie, ngauss, ieqn, ndof_el(2,0:imax+1)
data       ngauss/3/
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), rhs(g_gdof,g_neqns), s(g_neqns)
real*8  :: dx, xg, wi, xc, &
           carea(3), weight(3)

  call rutope(1, ngauss, carea, weight)

  !--- 2-exactness
  !----------
  if (iprob .eq. -1) then

    do ie = 0,imax+1

      xc = 0.5*(coord(ie+1)+coord(ie))
      dx = coord(ie+1)-coord(ie)

      rhs = 0.0

      do ig = 1,ngauss

        wi = 0.5 * weight(ig) * dx
        xg = carea(ig) * 0.5*dx + 0.5 * ( coord(ie+1)+coord(ie) )

        s = quadraticfn(xg, 0.0)

        do ieqn = 1,g_neqns
          rhs(1,ieqn) = rhs(1,ieqn) + wi * s(ieqn)
          rhs(2,ieqn) = rhs(2,ieqn) + wi * s(ieqn) * carea(ig)
        end do !ieqn

      end do !ig

      do ieqn = 1,g_neqns
        ucons(1,ieqn,ie) = rhs(1,ieqn) / dx
        ucons(2,ieqn,ie) = rhs(2,ieqn) / (dx/3.0)
        ucons(3,ieqn,ie) = 0.0
      end do !ieqn

    end do !ie

  else
    write(*,*) "Incorrect problem setup code!"
    call exit

  end if

  ndof_el(1,:) = g_gdof
  ndof_el(2,:) = 0

end subroutine init_soln_kex

!----------------------------------------------------------------------------------------------
!----- Solution initialization for multimaterial 6-eq pressure non-equilibrium 
!----- velocity equilibrium 2fluid model:
!----------------------------------------------------------------------------------------------

subroutine init_soln_mm6eq(reconst_mm6eq, ucons, uprim, matint_el, ndof_el)

procedure(), pointer :: reconst_mm6eq
integer :: matint_el(0:imax+1), ndof_el(2,0:imax+1), imat, ielem
real*8  :: ucons(g_tdof,g_neqns,0:imax+1), uprim(g_tdof,g_nprim,0:imax+1)
real*8  :: s(g_neqns), xf, p1l, p1r, t1l, t1r, &
           ul, ur, p2l, p2r, t2l, t2r, rho1, rho2

  !--- Gaussian
  !----------
  if (iprob .eq. -1) then

     g_alphamin = 1.d-10

     alpha_fs(1) = g_alphamin
     alpha_fs(2) = 1.0-alpha_fs(1)
     t_fs = 0.003484320557
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)

     call nondimen_mm6eq()

     if (g_nsdiscr .ge. 11) then
       call weakinit_p1(ucons)
       if (g_nsdiscr .eq. 12) ucons(3,:,:) = 0.0

     else
       do ielem = 1,imax
         xf = 0.5*(coord(ielem) + coord(ielem+1))

         s = gaussian(xf,0.0)
         ucons(1,:,ielem) = s(:)
       end do !ielem
       if (g_nsdiscr .eq. 1) ucons(2,:,:) = 0.0

     end if

  !--- MMS 1: advection of equilibrium interface
  !----------
  else if (iprob .eq. -2) then

     g_alphamin = 1.d-12

     alpha_fs(:) = g_alphamin
     alpha_fs(1) = 1.0 - dble(g_mmi%nummat-1)*g_alphamin
     if (g_mmi%nummat > 2) then
       write(*,*) " Error: 'tanh manufactured problem' not configured for &
         more than two materials"
       stop
     end if

     ! Quantities below calculated based on following:
     ! rhoe_1=1/rho_1, rhoe_2=1/rho_2, gam_1=gam_2, rho_1=rho_2=rho_01+k1*x
     ! rho_01=5
     u_fs = 1.0
     pr_fs = 0.4
     t_fs = 1.0/(5.0*g_cp(1)/g_gam(1))
     do imat = 1, g_mmi%nummat
       rhomat_fs(imat) = eos3_density(g_gam(imat), g_cp(imat), g_pc(imat), &
         pr_fs, t_fs)
     end do !imat

     call nondimen_mm6eq()

     if (g_nsdiscr .ge. 11) then
       call weakinit_p1(ucons)
       if (g_nsdiscr .eq. 12) ucons(3,:,:) = 0.0

     else
       do ielem = 1,imax
         xf = 0.5*(coord(ielem) + coord(ielem+1))

         s = mms_tanh(xf,0.0)
         ucons(1,:,ielem) = s(:)
       end do !ielem
       if (g_nsdiscr .eq. 1) ucons(2,:,:) = 0.0

     end if

  !--- MMS 2: two-material nonlinear energy growth
  !----------
  else if (iprob .eq. -3) then

     g_alphamin = 1.d-12

     alpha_fs(:) = g_alphamin
     alpha_fs(1) = 1.0 - dble(g_mmi%nummat-1)*g_alphamin
     if (g_mmi%nummat > 2) then
       write(*,*) " Error: 'nleg manufactured problem' not configured for &
         more than two materials"
       stop
     end if

     ! Quantities below calculated based on following:
     ! \rho_10 = \rho_20 = 30.0
     ! e_20 = 5.0
     ! c_1 = 0.5
     ! c_2 = 1.0
     ! c_3 = -1.0
     ! beta = 1.0
     ! k = 0.8
     u_fs = 1.0
     pr_fs = 1.0
     t_fs = 1.0
     rhomat_fs = 1.0
     call nondimen_mm6eq()

     pr_fs = eos3_pr(g_gam(1), g_pc(1), 30.0, 150.0, 0.0)
     t_fs = eos3_t(g_gam(1), g_cp(1), g_pc(1), 30.0, 150.0, 0.0)
     do imat = 1, g_mmi%nummat
       rhomat_fs(imat) = eos3_density(g_gam(imat), g_cp(imat), g_pc(imat), &
         pr_fs, t_fs)
     end do !imat

     if (g_nsdiscr .ge. 11) then
       call weakinit_p1(ucons)
       if (g_nsdiscr .eq. 12) ucons(3,:,:) = 0.0

     else
       do ielem = 1,imax
         xf = 0.5*(coord(ielem) + coord(ielem+1))

         s = mms_nleg(xf,0.0)
         ucons(1,:,ielem) = s(:)
       end do !ielem
       if (g_nsdiscr .eq. 1) ucons(2,:,:) = 0.0

     end if

  !--- SCD
  !----------
  else if (iprob .eq. 0) then

     g_alphamin = 1.d-12

     alpha_fs(1) = g_alphamin
     alpha_fs(2) = 1.0-alpha_fs(1)
     t_fs = 300.0
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     p2l = pr_fs
     t1l = t_fs
     t2l = t_fs
     ul  = u_fs
     ! right state
     p1r = pr_fs
     p2r = pr_fs
     t1r = t_fs
     t2r = t_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.5) then
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1l, t1l)
           rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), p2l, t2l)
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = 1.0-alpha_fs(1)
           ucons(1,3,ielem) = alpha_fs(1) * rho1
           ucons(1,4,ielem) = (1.0-alpha_fs(1)) * rho2
           ucons(1,5,ielem) = (ucons(1,3,ielem)+ucons(1,4,ielem)) * u_fs
           ucons(1,6,ielem) = alpha_fs(1) * eos3_rhoe(g_gam(1), g_pc(1), p1l, rho1, ul)
           ucons(1,7,ielem) = (1.0-alpha_fs(1)) * eos3_rhoe(g_gam(2), g_pc(2), p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1r, t1r)
           rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), p2r, t2r)
           ucons(1,1,ielem) = 1.0-alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(1)
           ucons(1,3,ielem) = (1.0-alpha_fs(1)) * rho1
           ucons(1,4,ielem) = alpha_fs(1) * rho2
           ucons(1,5,ielem) = (ucons(1,3,ielem)+ucons(1,4,ielem)) * u_fs
           ucons(1,6,ielem) = (1.0-alpha_fs(1)) * eos3_rhoe(g_gam(1), g_pc(1), p1r, rho1, ur)
           ucons(1,7,ielem) = alpha_fs(1) * eos3_rhoe(g_gam(2), g_pc(2), p2r, rho2, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- Sod Shocktube
  !----------
  else if (iprob .eq. 1) then

     g_alphamin = 1.d-12

     alpha_fs(1) = g_alphamin
     alpha_fs(2) = 1.0-alpha_fs(1)
     u_fs = 0.0
     pr_fs = 1.0
     t_fs = 3.484321d-3
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     p2l = pr_fs
     t1l = t_fs
     t2l = t_fs
     ul  = u_fs
     ! right state
     p1r = 0.1*pr_fs
     p2r = 0.1*pr_fs
     t1r = 0.8*t_fs
     t2r = 0.8*t_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.5) then
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1l, t1l)
           rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), p2l, t2l)
           ucons(1,1,ielem) = 1.0-alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(1)
           ucons(1,3,ielem) = ucons(1,1,ielem) * rho1
           ucons(1,4,ielem) = ucons(1,2,ielem) * rho2
           ucons(1,5,ielem) = (ucons(1,3,ielem)+ucons(1,4,ielem)) * u_fs
           ucons(1,6,ielem) = ucons(1,1,ielem) * eos3_rhoe(g_gam(1), g_pc(1), p1l, rho1, ul)
           ucons(1,7,ielem) = ucons(1,2,ielem) * eos3_rhoe(g_gam(2), g_pc(2), p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1r, t1r)
           rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), p2r, t2r)
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = 1.0-alpha_fs(1)
           ucons(1,3,ielem) = ucons(1,1,ielem) * rho1
           ucons(1,4,ielem) = ucons(1,2,ielem) * rho2
           ucons(1,5,ielem) = (ucons(1,3,ielem)+ucons(1,4,ielem)) * u_fs
           ucons(1,6,ielem) = ucons(1,1,ielem) * eos3_rhoe(g_gam(1), g_pc(1), p1r, rho1, ur)
           ucons(1,7,ielem) = ucons(1,2,ielem) * eos3_rhoe(g_gam(2), g_pc(2), p2r, rho2, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- Abgrall's water-air Shocktube
  !----------
  else if (iprob .eq. 2) then

     g_alphamin = 1.d-12

     alpha_fs(1) = g_alphamin
     alpha_fs(2) = 1.0-alpha_fs(1)
     u_fs = 0.0
     pr_fs = 1.0d9 !2.0d8
     t_fs = 494.646 !247.323
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     p2l = pr_fs
     t1l = t_fs
     t2l = t_fs
     ul  = u_fs
     ! right state
     p1r = pr_fs/1.0d4 !pr_fs/2.0d3
     p2r = pr_fs/1.0d4 !pr_fs/2.0d3
     t1r = 0.070441*t_fs !0.028176*t_fs
     t2r = 0.070441*t_fs !0.028176*t_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.75) then
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1l, t1l)
           rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), p2l, t2l)
           ucons(1,1,ielem) = 1.0-g_alphamin
           ucons(1,2,ielem) = g_alphamin
           ucons(1,3,ielem) = (1.0-g_alphamin) * rho1
           ucons(1,4,ielem) = g_alphamin * rho2
           ucons(1,5,ielem) = (ucons(1,3,ielem)+ucons(1,4,ielem)) * u_fs
           ucons(1,6,ielem) = (1.0-g_alphamin) * eos3_rhoe(g_gam(1), g_pc(1), p1l, rho1, ul)
           ucons(1,7,ielem) = g_alphamin * eos3_rhoe(g_gam(2), g_pc(2), p2l, rho2, ul)
        else
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1r, t1r)
           rho2 = eos3_density(g_gam(2), g_cp(2), g_pc(2), p2r, t2r)
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = 1.0-g_alphamin
           ucons(1,3,ielem) = g_alphamin * rho1
           ucons(1,4,ielem) = (1.0-g_alphamin) * rho2
           ucons(1,5,ielem) = (ucons(1,3,ielem)+ucons(1,4,ielem)) * u_fs
           ucons(1,6,ielem) = g_alphamin * eos3_rhoe(g_gam(1), g_pc(1), p1r, rho1, ur)
           ucons(1,7,ielem) = (1.0-g_alphamin) * eos3_rhoe(g_gam(2), g_pc(2), p2r, rho2, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- 3-material SCD
  !----------
  else if (iprob .eq. 3) then

     g_alphamin = 1.d-12

     alpha_fs(1) = 1.0-2.0*g_alphamin
     alpha_fs(2) = g_alphamin
     alpha_fs(3) = g_alphamin
     t_fs = 300.0
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)
     rhomat_fs(3) = eos3_density(g_gam(3), g_cp(3), g_pc(3), pr_fs, t_fs)

     call nondimen_mm6eq()

     p1l = pr_fs
     t1l = t_fs
     ul  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.25) then
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(2)
           ucons(1,3,ielem) = alpha_fs(3)

        elseif (xf .le. 0.35) then
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = 1.0-2.0*g_alphamin
           ucons(1,3,ielem) = g_alphamin

        else
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = g_alphamin
           ucons(1,3,ielem) = 1.0-2.0*g_alphamin

        end if

        do imat = 1,g_mmi%nummat
           ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rhomat_fs(imat)
           ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
             * eos3_rhoe(g_gam(imat), g_pc(imat), p1l, rhomat_fs(imat), ul)
        end do !imat
        ucons(1,g_mmi%imome,ielem) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ielem)) * u_fs

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- 3-material shocked problem
  !----------
  else if (iprob .eq. 4) then

     g_alphamin = 1.d-12

     alpha_fs(1) = 1.0-2.0*g_alphamin
     alpha_fs(2) = g_alphamin
     alpha_fs(3) = g_alphamin
     t_fs = 0.0348432
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)
     rhomat_fs(3) = eos3_density(g_gam(3), g_cp(3), g_pc(3), pr_fs, t_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     t1l = t_fs
     ul  = u_fs
     ! right state
     p1r = 0.1*pr_fs
     t1r = 0.1*t_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.4) then
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(2)
           ucons(1,3,ielem) = alpha_fs(3)
           do imat = 1,g_mmi%nummat
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rhomat_fs(imat)
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1l, rhomat_fs(imat), ul)
           end do !imat

        elseif (xf .le. 0.6) then
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = 1.0-2.0*g_alphamin
           ucons(1,3,ielem) = g_alphamin
           do imat = 1,g_mmi%nummat
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rhomat_fs(imat)
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1r, rhomat_fs(imat), ur)
           end do !imat

        else
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = g_alphamin
           ucons(1,3,ielem) = 1.0-2.0*g_alphamin
           do imat = 1,g_mmi%nummat
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rhomat_fs(imat)
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1r, rhomat_fs(imat), ur)
           end do !imat

        end if

        ucons(1,g_mmi%imome,ielem) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ielem)) * u_fs

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- 3-material shocked problem
  !----------
  else if (iprob .eq. 5) then

     g_alphamin = 1.d-08

     alpha_fs(1) = 1.0-2.0*g_alphamin
     alpha_fs(2) = g_alphamin
     alpha_fs(3) = g_alphamin
     u_fs = 0.0
     pr_fs = 1.0d9
     t_fs = 494.646
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)
     rhomat_fs(3) = eos3_density(g_gam(3), g_cp(3), g_pc(3), pr_fs, t_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     t1l = t_fs
     ul  = u_fs
     ! right state
     p1r = pr_fs/1.0d4
     t1r = 0.070441*t_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.6) then
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(2)
           ucons(1,3,ielem) = alpha_fs(3)
           do imat = 1,g_mmi%nummat
              rho1 = eos3_density(g_gam(imat), g_cp(imat), g_pc(imat), p1l, t1l)
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rho1
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1l, rho1, ul)
           end do !imat

        elseif (xf .le. 0.75) then
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = 1.0-2.0*g_alphamin
           ucons(1,3,ielem) = g_alphamin
           do imat = 1,g_mmi%nummat
              rho1 = eos3_density(g_gam(imat), g_cp(imat), g_pc(imat), p1r, t1r)
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rho1
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1r, rho1, ur)
           end do !imat

        else
           ucons(1,1,ielem) = g_alphamin
           ucons(1,2,ielem) = g_alphamin
           ucons(1,3,ielem) = 1.0-2.0*g_alphamin
           do imat = 1,g_mmi%nummat
              rho1 = eos3_density(g_gam(imat), g_cp(imat), g_pc(imat), p1r, t1r)
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rho1
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1r, rho1, ur)
           end do !imat

        end if

        ucons(1,g_mmi%imome,ielem) = sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ielem)) * u_fs

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- single-material Sod Shocktube
  !----------
  else if (iprob .eq. 6) then

     g_alphamin = 1.d-14

     alpha_fs(1) = 1.0
     u_fs = 0.0
     pr_fs = 1.0
     t_fs = 3.484321d-3
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     t1l = t_fs
     ul  = u_fs
     ! right state
     p1r = 0.1*pr_fs
     t1r = 0.8*t_fs
     ur  = u_fs

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 0.5) then
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1l, t1l)
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(1) * rho1
           ucons(1,3,ielem) = ucons(1,2,ielem) * u_fs
           ucons(1,4,ielem) = alpha_fs(1) * eos3_rhoe(g_gam(1), g_pc(1), p1l, rho1, ul)
        else
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1r, t1r)
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(1) * rho1
           ucons(1,3,ielem) = ucons(1,2,ielem) * u_fs
           ucons(1,4,ielem) = alpha_fs(1) * eos3_rhoe(g_gam(1), g_pc(1), p1r, rho1, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- 3-material cavitation tube problem
  !----------
  else if (iprob .eq. 7) then

     g_alphamin = 1.d-12

     alpha_fs(2) = 0.01
     alpha_fs(3) = 0.1
     alpha_fs(1) = 1.0-sum(alpha_fs(2:3))
     t_fs = 354.0
     rhomat_fs(1) = eos3_density(g_gam(1), g_cp(1), g_pc(1), pr_fs, t_fs)
     rhomat_fs(2) = eos3_density(g_gam(2), g_cp(2), g_pc(2), pr_fs, t_fs)
     rhomat_fs(3) = eos3_density(g_gam(3), g_cp(3), g_pc(3), pr_fs, t_fs)

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     t1l = t_fs
     ul  = -20.0/a_nd
     ! right state
     p1r = pr_fs
     t1r = t_fs
     ur  = 20.0/a_nd

     do ielem = 0,imax+1

        xf = coord(ielem)

        ucons(1,1,ielem) = alpha_fs(1)
        ucons(1,2,ielem) = alpha_fs(2)
        ucons(1,3,ielem) = alpha_fs(3)

        if (xf .le. 0.5) then
           do imat = 1,g_mmi%nummat
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rhomat_fs(imat)
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1l, rhomat_fs(imat), ul)
           end do !imat
           ucons(1,g_mmi%imome,ielem) = &
             sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ielem)) * ul

        else
           do imat = 1,g_mmi%nummat
              ucons(1,g_mmi%irmin+imat-1,ielem) = ucons(1,imat,ielem) * rhomat_fs(imat)
              ucons(1,g_mmi%iemin+imat-1,ielem) = ucons(1,imat,ielem) &
                * eos3_rhoe(g_gam(imat), g_pc(imat), p1r, rhomat_fs(imat), ur)
           end do !imat
           ucons(1,g_mmi%imome,ielem) = &
             sum(ucons(1,g_mmi%irmin:g_mmi%irmax,ielem)) * ur

        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  !--- Shu-Osher shock-entropy wave interaction
  !----------
  else if (iprob .eq. 8) then

     g_alphamin = 1.d-12

     alpha_fs(:) = g_alphamin
     alpha_fs(1) = 1.0 - dble(g_mmi%nummat-1)*g_alphamin
     if (g_mmi%nummat > 2) then
       write(*,*) " Error: shock-entropy wave interaction not configured for &
         more than two materials"
       stop
     end if
     u_fs = 2.629369
     pr_fs = 10.3333
     t_fs = 9.33450733d-3
     do imat = 1, g_mmi%nummat
       rhomat_fs(imat) = eos3_density(g_gam(imat), g_cp(imat), g_pc(imat), &
         pr_fs, t_fs)
     end do !imat

     call nondimen_mm6eq()

     if (g_nsdiscr .ge. 11) then
       call weakinit_p1(ucons)
       if (g_nsdiscr .eq. 12) ucons(3,:,:) = 0.0

     else
       do ielem = 1,imax
         xf = 0.5*(coord(ielem) + coord(ielem+1))

         s = shockentropywave(xf,0.0)
         ucons(1,:,ielem) = s(:)
       end do !ielem
       if (g_nsdiscr .eq. 1) ucons(2,:,:) = 0.0

     end if

  !--- LeBlanc shocktube
  !----------
  else if (iprob .eq. 9) then

     g_alphamin = 0.0 !1.d-12

     alpha_fs(:) = g_alphamin
     alpha_fs(1) = 1.0 - dble(g_mmi%nummat-1)*g_alphamin
     if (g_mmi%nummat > 1) then
       write(*,*) " Error: LeBlanc shocktube not configured for &
         more than one material"
       stop
     end if
     u_fs = 0.0
     pr_fs = eos3_pr(g_gam(1), g_pc(1), 1.0, 0.1, 0.0)
     t_fs = eos3_t(g_gam(1), g_cp(1), g_pc(1), 1.0, 0.1, 0.0)
     do imat = 1, g_mmi%nummat
       rhomat_fs(imat) = eos3_density(g_gam(imat), g_cp(imat), g_pc(imat), &
         pr_fs, t_fs)
     end do !imat

     call nondimen_mm6eq()

     ! left state
     p1l = pr_fs
     t1l = t_fs
     ul  = 0.0
     ! right state (need to non-dimensionalize)
     p1r = eos3_pr(g_gam(1), g_pc(1), 0.001/rho_nd, 1.0d-10/(rho_nd*p_nd), 0.0)
     t1r = eos3_t(g_gam(1), g_cp(1), g_pc(1), 0.001/rho_nd, &
       1.0d-10/(rho_nd*p_nd), 0.0)
     ur  = 0.0

     do ielem = 0,imax+1

        xf = coord(ielem)

        if (xf .le. 3.0) then
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1l, t1l)
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(1) * rho1
           ucons(1,3,ielem) = ucons(1,2,ielem) * u_fs
           ucons(1,4,ielem) = alpha_fs(1) * eos3_rhoe(g_gam(1), g_pc(1), p1l, rho1, ul)
        else
           rho1 = eos3_density(g_gam(1), g_cp(1), g_pc(1), p1r, t1r)
           ucons(1,1,ielem) = alpha_fs(1)
           ucons(1,2,ielem) = alpha_fs(1) * rho1
           ucons(1,3,ielem) = ucons(1,2,ielem) * u_fs
           ucons(1,4,ielem) = alpha_fs(1) * eos3_rhoe(g_gam(1), g_pc(1), p1r, rho1, ur)
        end if

        if (g_nsdiscr .ge. 1) then
          ucons(2,:,ielem) = 0.0
            if (g_nsdiscr .ge. 12) then
              ucons(3,:,ielem) = 0.0
            end if
        end if

     end do !ielem

  else
     write(*,*) "Incorrect problem setup code!"

  end if

  ! set ndof indicators
  ndof_el(1,:) = g_gdof
  ndof_el(2,:) = 0
  if (g_intreco > 0) call fill_ndofel(ucons, ndof_el)

  ! boundary conditions:
  call get_bc_mm6eq(ucons)
  call fill_matintel(ucons, matint_el)
  call weak_recons_primitives(ucons, uprim, ndof_el)
  call ignore_tinyphase_mm6eq(ucons, uprim)
  call reconst_mm6eq(ucons, uprim)

  call gnuplot_flow_mm6eq(ucons, uprim, matint_el, 0)
  call gnuplot_flow_p1_mm6eq(ucons, uprim, matint_el, 0)

end subroutine init_soln_mm6eq

!-------------------------------------------------------------------------------
!----- Weak initialization for DG(P1) for smooth problems:
!-------------------------------------------------------------------------------

subroutine weakinit_p1(ucons)

integer :: ig, ie, ieqn, ngauss
data       ngauss/5/

real*8  :: wi, vol, xc, x, s(g_neqns), rhs(g_tdof,g_neqns), carea(5), weight(5)
real*8  :: ucons(g_tdof,g_neqns,0:imax+1)

  call rutope(1, ngauss, carea, weight)

  do ie = 1,imax

    vol = coord(ie+1)-coord(ie)
    xc  = 0.5*(coord(ie+1)+coord(ie))

    rhs(:,:) = 0.0

    do ig = 1,ngauss

      wi = 0.5 * weight(ig) * vol
      x = carea(ig) * 0.5 * vol + xc
      if (iprob == -1) then
        s = gaussian(x,0.0)
      else if (iprob == -2) then
        s = mms_tanh(x,0.0)
      else if (iprob == -3) then
        s = mms_nleg(x,0.0)
      else if (iprob == 8) then
        s = shockentropywave(x, 0.0)
      end if

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

function quadraticfn(x,t)

real*8, intent(in)  :: x, t
real*8  :: xc, quadraticfn(g_neqns)

  xc  = 0.25
  quadraticfn(1) = 10.0 - 20.0*x*x + 10.0*x
  quadraticfn(2) = 100.0*x*x*x - 5.0*x - 100.0*x*x+ 10.0
  quadraticfn(3) = dexp( -(x-xc)*(x-xc)/(2.0 * 0.002) )

end function

!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow_mm6eq(ucons, uprim, matint_el, itstep)

integer, intent(in) :: matint_el(0:imax+1), itstep
real*8,  intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                       uprim(g_tdof,g_nprim,0:imax+1)

integer :: ielem, imat
real*8  :: xcc, pmix, tmix, rhomix, emix, temix, trcell
real*8  :: uconsi(g_neqns), uprimi(g_nprim), tmat(g_mmi%nummat)

character(len=100) :: filename2,filename3

associate (nummat=>g_mmi%nummat)

  write(filename2,'(1I50)')itstep
  filename3 = trim(adjustl(filename2)) // '.twofluid.'//'dat'
  open(23,file=trim(adjustl(filename3)),status='unknown')

  !--- write material and bulk meta-data to gnuplot file
  write(23,'(A8)',advance='no') "# xcc, "
  do imat = 1,nummat
     write(23,'(A12)',advance='no') "alphamat, "
  end do !imat
  write(23,'(4A8)',advance='no') "rhomix, ", &
                                 "umix, ", &
                                 "pmix, ", &
                                 "tmix, "
  do imat = 1,nummat
     write(23,'(A8)',advance='no') "pmat, "
  end do !imat
  do imat = 1,nummat
     write(23,'(A8)',advance='no') "tmat, "
  end do !imat
  write(23,'(5A15)') "e_m, ", &
                     "te_m," , &
                     "limcons_cell, ", &
                     "limprim_cell, ", &
                     "int_cell"

  do ielem = 1,imax

     uconsi = ucons(1,:,ielem)
     if (g_pureco == 1) then
       uprimi = uprim(1,:,ielem)
     else
       call get_uprim_mm6eq(uconsi, uprimi)
     end if

     xcc = 0.5d0 * (coord(ielem) + coord(ielem+1))

     rhomix = 0.0
     pmix = 0.0
     tmix = 0.0
     emix = 0.0
     temix = 0.0
     do imat = 1,nummat
        rhomix = rhomix + uconsi(g_mmi%irmin+imat-1)
        pmix = pmix + uprimi(apr_idx(nummat,imat))
        tmat(imat) = eos3_t(g_gam(imat), g_cp(imat), g_pc(imat), &
          uconsi(g_mmi%irmin+imat-1)/uconsi(imat), &
          uconsi(g_mmi%iemin+imat-1)/uconsi(imat), uprimi(vel_idx(nummat,0)))
        tmix = tmix + uconsi(imat)*tmat(imat)
        emix = emix + (uconsi(g_mmi%iemin+imat-1) &
                       - 0.5*uconsi(g_mmi%irmin+imat-1) &
                         *uprimi(vel_idx(nummat,0))*uprimi(vel_idx(nummat,0)))
        temix = temix + uconsi(g_mmi%iemin+imat-1)
     end do !imat
     emix = emix/rhomix
     temix = temix/rhomix

     !--- write material and bulk data to gnuplot file
     write(23,'(E16.6)',advance='no') xcc
     do imat = 1,nummat
        write(23,'(E16.6)',advance='no') uconsi(imat)
     end do !imat
     write(23,'(4E16.6)',advance='no') rhomix*rho_nd, &
                                       uprimi(vel_idx(nummat,0))*a_nd , &
                                       pmix*p_nd, &
                                       tmix*t_nd
     do imat = 1,nummat
        write(23,'(E16.6)',advance='no') uprimi(apr_idx(nummat,imat))*p_nd
     end do !imat
     do imat = 1,nummat
        write(23,'(E16.6)',advance='no') tmat(imat)*t_nd
     end do !imat
     write(23,'(5E16.6)') emix*p_nd/rho_nd, &
                          temix*p_nd/rho_nd, &
                          g_limcell(1,ielem), &
                          g_limcell(2,ielem), &
                          dble(matint_el(ielem))

  end do !ielem

  close(23)

end associate

end subroutine gnuplot_flow_mm6eq

!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow_p1_mm6eq(ucons, uprim, matint_el, itstep)

integer, intent(in) :: matint_el(0:imax+1), itstep
real*8,  intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), &
                       uprim(g_tdof,g_nprim,0:imax+1)

integer :: ielem, imat
real*8  :: dx, xc, xp, pmix, tmix, rhomix, emix, temix, trcell
real*8  :: uconsi(g_neqns), uprimi(g_nprim), basis(g_tdof), tmat(g_mmi%nummat)

character(len=100) :: filename2,filename3

associate (nummat=>g_mmi%nummat)

  write(filename2,'(1I50)')itstep
  filename3 = trim(adjustl(filename2)) // '.dgtwofluid.'//'dat'
  open(24,file=trim(adjustl(filename3)),status='unknown')

  write(filename2,'(1I50)')itstep
  filename3 = trim(adjustl(filename2)) // '.conserved.'//'dat'
  open(25,file=trim(adjustl(filename3)),status='unknown')

  !--- write material and bulk meta-data to gnuplot file
  write(24,'(A8)',advance='no') "# xcc, "
  do imat = 1,nummat
     write(24,'(A12)',advance='no') "alphamat, "
  end do !imat
  write(24,'(4A8)',advance='no') "rhomix, ", &
                                 "umix, ", &
                                 "pmix, ", &
                                 "tmix, "
  do imat = 1,nummat
     write(24,'(A8)',advance='no') "pmat, "
  end do !imat
  do imat = 1,nummat
     write(24,'(A8)',advance='no') "tmat, "
  end do !imat
  write(24,'(3A8)') "e_m, ", &
                    "te_m, ", &
                    "int_cell"

  write(25,'(A8)',advance='no') "# xcc, "
  do imat = 1,nummat
     write(25,'(A12)',advance='no') "rhomat, "
  end do !imat
  do imat = 1,nummat
     write(25,'(A12)',advance='no') "temat, "
  end do !imat
  write(25,*) " "

  do ielem = 1,imax

     ! cell geometry
     dx = coord(ielem+1) - coord(ielem)
     xc = coord(ielem) + 0.5 * dx

     ! left face
     xp = coord(ielem)
     call get_basisfns(xp, xc, dx, basis)
     call ho_reconstruction(g_neqns, ucons(:,:,ielem), basis, uconsi(:))
     if (g_pureco == 1) then
       call ho_reconstruction(g_nprim, uprim(:,:,ielem), basis, uprimi(:))
     else
       call get_uprim_mm6eq(uconsi, uprimi)
     end if

     rhomix = 0.0
     pmix = 0.0
     tmix = 0.0
     emix = 0.0
     temix = 0.0
     do imat = 1,nummat
        rhomix = rhomix + uconsi(g_mmi%irmin+imat-1)
        pmix = pmix + uprimi(apr_idx(nummat,imat))
        tmat(imat) = eos3_t(g_gam(imat), g_cp(imat), g_pc(imat), &
          uconsi(g_mmi%irmin+imat-1)/uconsi(imat), &
          uconsi(g_mmi%iemin+imat-1)/uconsi(imat), uprimi(vel_idx(nummat,0)))
        tmix = tmix + uconsi(imat)*tmat(imat)
        emix = emix + (uconsi(g_mmi%iemin+imat-1) &
                       - 0.5*uconsi(g_mmi%irmin+imat-1) &
                         *uprimi(vel_idx(nummat,0))*uprimi(vel_idx(nummat,0)))
        temix = temix + uconsi(g_mmi%iemin+imat-1)
     end do !imat
     emix = emix/rhomix
     temix = temix/rhomix

     trcell = dble(matint_el(ielem))

     !--- write material and bulk data to gnuplot file
     write(24,'(E16.6)',advance='no') xp
     do imat = 1,nummat
        write(24,'(E16.6)',advance='no') uconsi(imat)
     end do !imat
     write(24,'(4E16.6)',advance='no') rhomix*rho_nd, &
                                       uprimi(vel_idx(nummat,0))*a_nd , &
                                       pmix*p_nd, &
                                       tmix*t_nd
     do imat = 1,nummat
        write(24,'(E16.6)',advance='no') uprimi(apr_idx(nummat,imat))*p_nd
     end do !imat
     do imat = 1,nummat
        write(24,'(E16.6)',advance='no') tmat(imat)*t_nd
     end do !imat
     write(24,'(3E16.6)') emix*p_nd/rho_nd, &
                          temix*p_nd/rho_nd, &
                          trcell

     write(25,'(E16.6)',advance='no') xp
     do imat = 1,nummat
        write(25,'(E16.6)',advance='no') uconsi(g_mmi%irmin+imat-1)*rho_nd
     end do !imat
     do imat = 1,nummat
        write(25,'(E16.6)',advance='no') uconsi(g_mmi%iemin+imat-1)*p_nd
     end do !imat
     write(25,*) " "

     ! right face
     xp = coord(ielem+1)
     call get_basisfns(xp, xc, dx, basis)
     call ho_reconstruction(g_neqns, ucons(:,:,ielem), basis, uconsi(:))
     if (g_pureco == 1) then
       call ho_reconstruction(g_nprim, uprim(:,:,ielem), basis, uprimi(:))
     else
       call get_uprim_mm6eq(uconsi, uprimi)
     end if

     rhomix = 0.0
     pmix = 0.0
     tmix = 0.0
     emix = 0.0
     temix = 0.0
     do imat = 1,nummat
        rhomix = rhomix + uconsi(g_mmi%irmin+imat-1)
        pmix = pmix + uprimi(apr_idx(nummat,imat))
        tmat(imat) = eos3_t(g_gam(imat), g_cp(imat), g_pc(imat), &
          uconsi(g_mmi%irmin+imat-1)/uconsi(imat), &
          uconsi(g_mmi%iemin+imat-1)/uconsi(imat), uprimi(vel_idx(nummat,0)))
        tmix = tmix + uconsi(imat)*tmat(imat)
        emix = emix + (uconsi(g_mmi%iemin+imat-1) &
                       - 0.5*uconsi(g_mmi%irmin+imat-1) &
                         *uprimi(vel_idx(nummat,0))*uprimi(vel_idx(nummat,0)))
        temix = temix + uconsi(g_mmi%iemin+imat-1)
     end do !imat
     emix = emix/rhomix
     temix = temix/rhomix

     trcell = dble(matint_el(ielem))

     !--- write material and bulk data to gnuplot file
     write(24,'(E16.6)',advance='no') xp
     do imat = 1,nummat
        write(24,'(E16.6)',advance='no') uconsi(imat)
     end do !imat
     write(24,'(4E16.6)',advance='no') rhomix*rho_nd, &
                                       uprimi(vel_idx(nummat,0))*a_nd , &
                                       pmix*p_nd, &
                                       tmix*t_nd
     do imat = 1,nummat
        write(24,'(E16.6)',advance='no') uprimi(apr_idx(nummat,imat))*p_nd
     end do !imat
     do imat = 1,nummat
        write(24,'(E16.6)',advance='no') tmat(imat)*t_nd
     end do !imat
     write(24,'(3E16.6)') emix*p_nd/rho_nd, &
                          temix*p_nd/rho_nd, &
                          trcell

     write(25,'(E16.6)',advance='no') xp
     do imat = 1,nummat
        write(25,'(E16.6)',advance='no') uconsi(g_mmi%irmin+imat-1)*rho_nd
     end do !imat
     do imat = 1,nummat
        write(25,'(E16.6)',advance='no') uconsi(g_mmi%iemin+imat-1)*p_nd
     end do !imat
     write(25,*) " "

     write(24,*) " "

     write(25,*) " "

  end do !ielem

  close(24)
  close(25)

end associate

end subroutine gnuplot_flow_p1_mm6eq

!----------------------------------------------------------------------------------------------

subroutine gnuplot_diagnostics_mm6eq(ucons, uprim, cons_err, time)

real*8,  intent(in) :: time, ucons(g_tdof,g_neqns,0:imax+1), &
                       uprim(g_tdof,g_nprim,0:imax+1), cons_err(2)

logical :: mixed_cell
integer :: ie, k
real*8 :: alk, ap(g_nprim), pb, denob, pavg(g_mmi%nummat), deno(g_mmi%nummat)

associate (nummat=>g_mmi%nummat)

  pavg = 0.0
  deno = 0.0
  pb = 0.0
  denob = 0.0
  do ie=1,imax
    mixed_cell = .false.
    if (g_pureco == 1) then
      ap = uprim(1,:,ie)
    else
      call get_uprim_mm6eq(ucons(1,:,ie), ap)
    end if
    do k = 1,nummat
      alk = ucons(1,g_mmi%iamin+k-1,ie)
      if ((alk > 1d-2) .and. (alk < 1.0-1d-2)) then
        mixed_cell = .true.
        pavg(k) = pavg(k) + ap(k)
        deno(k) = deno(k) + alk
      end if
    end do !k
    if (mixed_cell) then
      pb = pb + sum(ap(apr_idx(nummat,1):apr_idx(nummat,nummat)))
      denob = denob + 1.0
    end if
  end do !ie

  do k = 1,nummat
  pavg(k) = pavg(k) / max(1d-14, deno(k))
  end do !k
  pb = pb / max(1d-14, denob)

  write(33,'(3E16.6)',advance='no') time, & !1
    cons_err(1), &                          !2
    cons_err(2)                             !3
  do k = 1,nummat
     write(33,'(E16.6)',advance='no') pavg(k)*p_nd
  end do !k
  write(33,*) pb*p_nd

end associate

end subroutine gnuplot_diagnostics_mm6eq

!----------------------------------------------------------------------------------------------

subroutine errorcalc_p1(ucons, t, err_log, linfty)

real*8, intent(in) :: ucons(g_tdof,g_neqns,0:imax+1), t

integer :: ig, ie, ieqn, imat, ngauss
logical :: file_exists
data       ngauss/3/

real*8  :: dx, xg, wi, xc, b3i, &
           u(g_neqns), &
           carea(3), weight(3), &
           s(g_neqns), &
           err, linfty, err_log(2,g_neqns), &
           errrho, errte

  call rutope(1, ngauss, carea, weight)

  err_log = 0.0
  linfty = 0.d0

  if ( (iprob .eq. -1) .or. (iprob .eq. -2) .or. (iprob .eq. -3) &
    .or. (iprob .eq. 0) ) then

    do ie = 1,imax
    do ig = 1,ngauss

      dx = coord(ie+1)-coord(ie)
      xc = 0.5*(coord(ie+1)+coord(ie))
      wi = 0.5 * weight(ig) * dx
      xg = carea(ig) * 0.5*dx + 0.5 * ( coord(ie+1)+coord(ie) )

      ! basis function

      if (g_nsdiscr .ge. 12) then
        b3i = p2basis(xg,xc,dx)
        do ieqn = 1,g_neqns
          u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie) + b3i*ucons(3,ieqn,ie)
        end do !ieqn
      else if (g_nsdiscr .ge. 1) then
        do ieqn = 1,g_neqns
          u(ieqn) = ucons(1,ieqn,ie) + carea(ig) * ucons(2,ieqn,ie)
        end do !ieqn
      else
        u(:) = ucons(1,:,ie)
      end if

      if (i_system .eq. 1) then
        if (iprob.eq.-1) s = gaussian(xg,t)
        if (iprob.eq.-2) s = mms_tanh(xg,t)
        if (iprob.eq.-3) s = mms_nleg(xg,t)
        if (iprob.eq.0) s = contact(xg,t)
      else if (i_system .eq. -1) then
        s = quadraticfn(xg,t)
      end if

      do ieqn = 1,g_neqns
        err = u(ieqn) - s(ieqn)
        err_log(1,ieqn) = err_log(1,ieqn) + wi*dabs(err)
        err_log(2,ieqn) = err_log(2,ieqn) + wi*err*err
      end do !ieqn

      linfty = max(linfty, dabs(u(1) - s(1)))

    end do !ig
    end do !ie

    !--- final error computation
    do ieqn = 1,g_neqns
      err_log(1,ieqn) = dlog10(err_log(1,ieqn))
      err_log(2,ieqn) = dsqrt(err_log(2,ieqn))
      err_log(2,ieqn) = dlog10(err_log(2,ieqn))
    end do !ieqn

    !--- write errors to file
    inquire(file='logl1errors.dat', exist=file_exists)
    if (file_exists) then
      open(41,file='logl1errors.dat',status='old', position='append')
    else
      open(41,file='logl1errors.dat',status='new')
      write(41,'(A13)') "#  l1-errors: "
      write(41,'(A31)') "#  dx,  ndof,  errors(nunk), ..."
    end if
    write(41,'(2F16.6)',advance='no') dlog10(coord(2)-coord(1)), &
      dlog10(dble(imax*g_gdof))
    do ieqn = 1,g_neqns
      write(41,'(F16.6)',advance='no') err_log(1,ieqn)
    end do !ieqn
    close(41)

    inquire(file='logl2errors.dat', exist=file_exists)
    if (file_exists) then
      open(42,file='logl2errors.dat',status='old', position='append')
    else
      open(42,file='logl2errors.dat',status='new')
      write(42,'(A13)') "#  l2-errors: "
      write(42,'(A31)') "#  dx,  ndof,  errors(nunk), ..."
    end if
    write(42,'(2F16.6)',advance='no') dlog10(coord(2)-coord(1)), &
      dlog10(dble(imax*g_gdof))
    do ieqn = 1,g_neqns
      write(42,'(F16.6)',advance='no') err_log(2,ieqn)
    end do !ieqn
    close(42)

    inquire(file='absl1errors.dat', exist=file_exists)
    if (file_exists) then
      open(43,file='absl1errors.dat',status='old', position='append')
    else
      open(43,file='absl1errors.dat',status='new')
      write(43,'(A13)') "#  l1-errors: "
      write(43,'(A31)') "#  dx,  ndof,  errors(nunk), ..."
    end if
    write(43,'(E16.6)',advance='no') coord(2)-coord(1)
    write(43,'(I16)',advance='no') imax*g_gdof
    do ieqn = 1,g_neqns
      write(43,'(E16.6)',advance='no') 10.0**err_log(1,ieqn)
    end do !ieqn
    close(43)

    inquire(file='absl2errors.dat', exist=file_exists)
    if (file_exists) then
      open(44,file='absl2errors.dat',status='old', position='append')
    else
      open(44,file='absl2errors.dat',status='new')
      write(44,'(A13)') "#  l2-errors: "
      write(44,'(A31)') "#  dx,  ndof,  errors(nunk), ..."
    end if
    write(44,'(E16.6)',advance='no') coord(2)-coord(1)
    write(44,'(I16)',advance='no') imax*g_gdof
    do ieqn = 1,g_neqns
      write(44,'(E16.6)',advance='no') 10.0**err_log(2,ieqn)
    end do !ieqn
    close(44)

  else

    write(*,*) "  WARNING: Error computation not configured for iprob = ", iprob

  end if

end subroutine errorcalc_p1

!----------------------------------------------------------------------------------------------

END MODULE preprocess
