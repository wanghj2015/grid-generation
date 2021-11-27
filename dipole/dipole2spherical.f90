subroutine dipole2spherical(npts, nlp, nmp, mu, chi, r, theta)

use mod_precision, only : wp=>dp
use mod_const,     only : pi, dtr

implicit none

integer, intent(in) :: npts, nlp, nmp

! dipole coordinates (mu,chi,phi)

real(wp), dimension(npts,nlp,nmp), intent(in) :: mu, chi

! spherical coordinates (r,theta,phi)

real(wp), dimension(npts,nlp,nmp), intent(out) :: r, theta

! local variables

integer :: i, j, k, midpoint

real(wp) :: zeta, gama, c1, c2, w, u


midpoint = (npts+1)/2

! transformation to spherical coordinates

! r is distance from Earth's center (normalized by its radius 1Re)
! theta is colatitude
c1 = 2.0_wp**(7.0_wp/3.0_wp)*3.0_wp**(-1.0_wp/3.0_wp)
c2 = 2.0_wp**(1.0_wp/3.0_wp)*3.0_wp**(2.0_wp/3.0_wp)
do i = 1, nmp
do j = 1, nlp
do k = 1, npts
zeta = (mu(k,j,i)/chi(k,j,i)**2)**2
gama = (9.0_wp*zeta + sqrt(3.0_wp) * &
       sqrt(27.0_wp*zeta**2+256.0_wp*zeta**3))**(1.0_wp/3.0_wp)
w    = -c1/gama + gama/c2/zeta
u    = -sqrt(w)/2.0_wp + sqrt(-w+2.0_wp/zeta/sqrt(w))/2.0_wp

r(k,j,i)     = u/chi(k,j,i)
theta(k,j,i) = asin(sqrt(u))

if (k == midpoint) then
   r(k,j,i)     = 1.0_wp/chi(k,j,i)
   theta(k,j,i) = asin(sqrt(1.0_wp))
endif

if (k > midpoint) then
   theta(k,j,i) = pi - theta(k,j,i)
endif

enddo
enddo
enddo


!print *, ' r = ', r(1,:,1)
!print *, ' r [midpoint] = ', r(midpoint,:,1)
!print *, ' r = ', r(npts,:,1)

!print *, ' theta = ', theta(:,1,1)/dtr
!print *, ' theta = ', theta(1,:,1)/dtr
!print *, ' theta [midpoint] = ', theta(midpoint,:,1)/dtr
!print *, ' theta = ', theta(npts,:,1)/dtr


end subroutine dipole2spherical

