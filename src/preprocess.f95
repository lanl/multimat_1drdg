!!---------------------------------------------------------------------------------------
!!----- Non-equilibrium Isothermal Multiphase (Preprocessing Module)
!!----- by
!!----- Aditya K Pandare
!!----- Department of Mechanical and Aerospace Engineering,
!!----- North Carolina State University.
!!---------------------------------------------------------------------------------------

MODULE preprocess

USE glob_var
USE rhs_flux
USE eos

implicit none

CONTAINS

!----------------------------------------------------------------------------------------------
!----- Read the control file:
!----------------------------------------------------------------------------------------------

subroutine read_cntl()

integer :: i

        open(12, file = 'setflow.cntl')

        read(12,*) iprob
        read(12,*) imax
        read(12,*) i_restart
        read(12,*) i_flux
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
!----- Initialize the solution:
!----------------------------------------------------------------------------------------------

subroutine init_soln(uprim, uprimn, ucons, uconsn)

integer :: i, ielem
real*8  :: uprim(ndof,neqns,0:imax+1), uprimn(ndof,neqns,0:imax+1), &
           ucons(neqns,0:imax+1), uconsn(neqns,0:imax+1)
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
        call get_bc(ucons)
        call decode_uprim(ucons,uprim)

        uprimn(:,:,:) = uprim(:,:,:)
        uconsn(:,:)   = ucons(:,:)

        call gnuplot_flow(uprim, ucons, 0)

end subroutine init_soln

!----------------------------------------------------------------------------------------------
!----- GNUPLOT outputs:
!----------------------------------------------------------------------------------------------

subroutine gnuplot_flow(uprim, ucons, itstep)

integer, intent(in) :: itstep
real*8,  intent(in) :: uprim(ndof,neqns,0:imax+1), ucons(neqns,0:imax+1)

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

end subroutine gnuplot_flow

!----------------------------------------------------------------------------------------------

END MODULE preprocess
