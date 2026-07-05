subroutine dipole_psi(npts, nlp, nmp, theta1, theta2, aa, ab, &
                      psi, mu, chi, phi, dpsi, dchi)

! modified orthogonal dipole coordinates (psi, chi, phi):
! Kageyama, A., Sugiyama, T., Watanabe, K., & Sato, T. (2006). 
! A note on the dipole coordinates. Computers and Geosciences, 32(2), 265–269. 
! https://doi.org/10.1016/j.cageo.2005.06.006

use mod_precision, only : wp=>dp
use mod_const,     only : dtr

implicit none

integer, intent(in) :: npts, nlp, nmp

real(wp), intent(in) :: theta1, theta2, aa, ab

! dipole coordinates (mu,chi,phi)

real(wp), dimension(npts,nlp,nmp), intent(out) :: mu, chi, phi

! modified dipole coordinates (psi,chi,phi)

real(wp), dimension(npts,nlp,nmp), intent(out) :: psi

real(wp), intent(out) :: dpsi, dchi

! local variables

integer :: i, j, k, midpoint

real(wp) :: r1, r2, mu1, mu2, chi1, chi2
real(wp) :: psi1, psi2


midpoint = (npts+1)/2

! dipole coordinates

! phi is longitude 
phi = 0.0d0

! assume that at r=1 the field lines (const chi) are between theta1 and theta2
! then r at the equator (theta=90) are (from constancy of chi) respectively

r1 = 1.0d0 / sin(theta1*dtr)**2
r2 = 1.0d0 / sin(theta2*dtr)**2


!print *, ' r1, r2 = ', r1, r2


! const-chi curves in a meridian plane (phi=const), denoting dipole field lines

chi1 = 1.0d0/r1
chi2 = 1.0d0/r2
dchi = (chi2-chi1)/(nlp-1)
do i = 1, nmp
do k = 1, npts
chi(k,1,i) = chi1
do j = 2, nlp
chi(k,j,i) = chi(k,j-1,i) + dchi
enddo
enddo
enddo


!print *, ' chi [midpoint] = ', chi(midpoint,:,1)
!print *, ' chi1, chi2, dchi = ', chi1, chi2, dchi


! modified orthogonal dipole coordinates (psi, chi, phi)

mu1  = -cos(theta1*dtr)
psi1 = log(aa*mu1 + sqrt(1.0d0 + (aa*mu1)**2))/ab
psi2 = 0.0d0
dpsi = (psi2-psi1)/((npts-1)/2)
do i = 1, nmp
do j = 1, nlp
! at equator
psi(midpoint,j,i) = 0.0d0
do k = midpoint-1, 1, -1
psi(k,j,i) = psi(k+1,j,i) - dpsi
enddo
do k = midpoint+1, npts
psi(k,j,i) = psi(k-1,j,i) + dpsi
enddo
enddo
enddo

print *, ' mu1 = ', mu1

!print *, ' psi = ', psi(:,1,1)
!print *, ' psi [midpoint] = ', psi(midpoint,:,1)
 print *, ' psi1, psi2, dpsi = ', psi1, psi2, dpsi


! mu is a potential function of a dipole field (B ~ grad mu)
do i = 1, nmp
do j = 1, nlp
do k = 1, npts
mu(k,j,i) = sinh(ab*psi(k,j,i))/aa
enddo
enddo
enddo


!print *, ' mu = ', mu(:,1,1)


end subroutine dipole_psi

