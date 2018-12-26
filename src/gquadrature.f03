MODULE gquadrature

implicit none

CONTAINS
!
!----------------------------------------------------------------------------
!
subroutine rutope(ndimn, ngaus, posgp, weigp)

integer, intent(in) :: ndimn, ngaus
real*8 :: posgp(ndimn,*), weigp(*)
!
!     This routine sets up the integration constants of open rules for
!     triangles and tetrahedra
!
!     NDIMN = 1             NDIME = 2             NDIME = 3
!
!     NGAUS  EXACT POL.    NGAUS  EXACT POL.     NGAUS  EXACT POL. 
!     -----  ---------     -----  ----------     -----  ----------
!         1         p1         1       p1            1       p1
!         2         p3         3       p2            4       p2
!         3         p5         4       p3            5       p3
!                              6       p4           11       p4
!                              7       p5           14       p5
!                             13       p9
!
!...  Line integral (the same as for brick elements)
!
  if(ngaus .eq. 1) then
    posgp(1,1) = 0.0
    weigp(  1) = 2.0
  else if(ngaus .eq. 2) then
    posgp(1,1) =-0.577350269189626
    posgp(1,2) = 0.577350269189626
    weigp(  1) = 1.0
    weigp(  2) = 1.0
  else if(ngaus .eq. 3) then
    posgp(1,1) =-0.774596669241483
    posgp(1,2) = 0.0
    posgp(1,3) = 0.774596669241483
    weigp(  1) = 0.555555555555556
    weigp(  2) = 0.888888888888889
    weigp(  3) = 0.555555555555556
  else
    write(*,*) "Quadrature unavailable for number of points:", ngaus
  end if

end subroutine rutope

!----------------------------------------------------------------------

END MODULE gquadrature
