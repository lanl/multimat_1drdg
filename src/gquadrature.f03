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
  else if(ngaus .eq. 4) then
    posgp(1,1) =-0.861136311594053
    posgp(1,2) =-0.339981043584856
    posgp(1,3) = 0.339981043584856
    posgp(1,4) = 0.861136311594053
    weigp(  1) = 0.347854845137454
    weigp(  2) = 0.652145154862546
    weigp(  3) = 0.652145154862546
    weigp(  4) = 0.347854845137454
  else if(ngaus .eq. 5) then
    posgp(1,1) =-0.906179845938666
    posgp(1,2) =-0.538469310105681
    posgp(1,3) = 0.0
    posgp(1,4) = 0.538469310105681
    posgp(1,5) = 0.906179845938666
    weigp(  1) = 0.236926885056190
    weigp(  2) = 0.478628670499366
    weigp(  3) = 0.568888888888889
    weigp(  4) = 0.478628670499366
    weigp(  5) = 0.236926885056190
  else
    write(*,*) "Quadrature unavailable for number of points:", ngaus
  end if

end subroutine rutope

!----------------------------------------------------------------------

END MODULE gquadrature
