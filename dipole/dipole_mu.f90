subroutine dipole_mu(npts, nlp, nmp, theta1, theta2, mu, chi, phi, dmu, dchi)

! standard orthogonal dipole coordinates (mu, chi, phi):
! Kageyama, A., Sugiyama, T., Watanabe, K., & Sato, T. (2006). 
! A note on the dipole coordinates. Computers and Geosciences, 32(2), 265–269. 
! https://doi.org/10.1016/j.cageo.2005.06.006

use mod_precision, only : wp=>dp
use mod_const,     only : dtr

implicit none

integer, intent(in) :: npts, nlp, nmp

real(wp), intent(in) :: theta1, theta2

! dipole coordinates (mu,chi,phi)

real(wp), dimension(npts,nlp,nmp), intent(out) :: mu, chi, phi

real(wp), intent(out) :: dmu, dchi

! local variables

integer :: i, j, k, midpoint

real(wp) :: r1, r2, mu1, mu2, chi1, chi2


midpoint = (npts+1)/2

! dipole coordinates

! phi is longitude 

phi = 0.

! assume that at r=1 the field lines (const chi) are between theta1 and theta2
! then r at the equator (theta=90) are (from constancy of chi) respectively

r1 = 1./sin(theta1*dtr)**2
r2 = 1./sin(theta2*dtr)**2


!print *, ' r1, r2 = ', r1, r2


! const-chi curves in a meridian plane (phi=const), denoting dipole field lines

chi1 = 1./r1
chi2 = 1./r2
dchi = (chi2-chi1)/(nlp-1)
do i = 1, nmp
do k = 1, npts
chi(k,1,i) = chi1
do j = 2, nlp
chi(k,j,i) = chi(k,j-1,i) + dchi
enddo
enddo
enddo


!print *, ' chi = ', chi(midpoint,:,1)
!print *, ' chi1, chi2, dchi = ', chi1, chi2, dchi


! mu is a potential function of a dipole field (B ~ grad mu)
mu1 = -cos(theta1*dtr)
mu2 = 0.
dmu = (mu2-mu1)/((npts-1)/2.)
do i = 1, nmp
do j = 1, nlp
! at equator
mu(midpoint,j,i) = 0.
do k = midpoint-1, 1, -1
mu(k,j,i) = mu(k+1,j,i) - dmu
enddo
do k = midpoint+1, npts
mu(k,j,i) = mu(k-1,j,i) + dmu
enddo
enddo
enddo


!print *, ' mu = ', mu(midpoint,:,1)
!print *, ' mu = ', mu(:,1,1)
!print *, ' mu1, mu2, dmu = ', mu1, mu2, dmu

end subroutine dipole_mu

