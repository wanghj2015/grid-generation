program dipole_coordinates

! Grid generation for an axial-centered dipole magnetic field
!
! Author: Houjun Wang, wanghoujun@gmail.com, 2012-2021

use mod_precision, only : wp=>dp
use mod_const,     only : pi, dtr
use mod_igrf,      only : get_igrf_coef, gx

implicit none

! namelist parameters

! use_psi = 0: use standard dipole coordinate 'mu' 
!           1: use modified dipole coordinate 'psi'
integer :: nptx, nlp, nmp, use_psi, use_dipole

! at r=1 field lines are between colatitudes 'theta1' & 'theta2'
real(wp) :: epoch, theta1, theta2, hmin, ds_step, aa

integer :: midpoint

! Earth’s mean reference spherical radius [km]
real(wp), parameter :: Re = 6371.2_wp

! dipole coordinates (mu,chi,phi)
real(wp), dimension(:,:,:), allocatable :: mu, chi, phi

! modified dipole coordinates (psi,chi,phi)
real(wp), dimension(:,:,:), allocatable :: psi

! spherical coordinates (r,theta,phi)

real(wp), dimension(:,:,:), allocatable :: r, theta

! start (nh) and end (sh) points on each fluxtube
integer, dimension(:,:), allocatable :: jn, js

! end/number of 'l' points on the constant mu/psi curves
integer, dimension(:,:), allocatable :: lptx

integer, allocatable :: lpts(:,:)

! magnetic field 'b' and arc length 's' along fluxtubes
type vc
   real(wp), allocatable :: vc1d(:)
end type vc

type(vc), dimension(:,:), allocatable :: b, s, h1, h2, h3, x, y
type(vc), dimension(:,:), allocatable :: para, perp


real(wp) :: dmu, dchi
real(wp) :: ab, dx, mindx, dpsi

real(wp), allocatable :: ell(:)

real(wp) :: gm

integer :: i, j, k, kk, npts
integer :: l, m

integer, allocatable :: klm(:,:)

namelist /grid_params/ epoch, nptx, nlp, nmp, theta1, theta2, hmin, &
                       ds_step, use_psi, aa, use_dipole


! get namelist

call read_namelist

! allocate space

allocate (mu(nptx,nlp,nmp), chi(nptx,nlp,nmp), phi(nptx,nlp,nmp))
allocate (psi(nptx,nlp,nmp))
allocate (r(nptx,nlp,nmp), theta(nptx,nlp,nmp))
allocate (jn(nlp,nmp), js(nlp,nmp))
allocate (lptx(nptx,nmp), ell(nlp), klm(nlp,nmp))

midpoint=(nptx+1)/2


! get the IGRF coefficients for the epoch
call get_igrf_coef (epoch)

! dipole moment for the corresponding epoch
gm = abs(gx(1,0))

print *, ' gm = ', gm

 
if (use_psi == 0) then

   ! dipole coordinates

   call dipole_mu(nptx, nlp, nmp, theta1, theta2, mu, chi, phi, dmu, dchi)

else

   ! modified dipole coordinates

   ab = log(aa + sqrt(1.0d0 + aa**2))

   print *, ' aa, ab = ', aa, ab

   call dipole_psi(nptx, nlp, nmp, theta1, theta2, aa, ab, &
                   psi, mu, chi, phi, dpsi, dchi)

   print *, ' dpsi, dchi = ', dpsi, dchi

endif


! (inverse) transformation to spherical coordinates

! r is distance from Earth's center (normalized by its radius 1Re)
! theta is colatitude

call dipole2spherical(nptx, nlp, nmp, mu, chi, r, theta)

!do m = 1, nmp
!do l = 1, nlp
!print *, ' l, r     = ', l, pack(r(:,l,m), r(:,l,m) /= 0)*371.2
!print *, ' l, theta = ', l, pack(theta(:,l,m)/dtr, theta(:,l,m) /= 0)
!print *, ' l, phi   = ', l, pack(phi(:,l,m)/dtr, phi(:,l,m) /= 0)
!enddo
!enddo


! r in [km]
r = r * Re


lptx = 0

do i = 1, nmp
do k = 1, nptx

do j = nlp, 1, -1
if (r(k,j,i) > Re + hmin) then
   lptx(k,i) = j

   !print *, ' k, lptx = ', k, lptx(k,i), r(k,lptx(k,i),i)-6371.2

   exit
endif
enddo

!print *, ' k, lptx = ', k, lptx(k,i), r(k,lptx(k,i),i)-6371.2

enddo
enddo


do i = 1, nmp
do j = 1, nlp

do k = 1, midpoint
if (r(k,j,i) > Re + hmin .and. lptx(k,i) > 1) then
   jn(j,i) = k
   exit
endif
enddo

do k = nptx, midpoint, -1
if (r(k,j,i) > Re + hmin .and. lptx(k,i) > 1) then
   js(j,i) = k
   exit
endif
enddo

!print *, ' j, jn, js = ', j, jn(j,i), js(j,i), r(jn(j,i),j,i)-6371.2, r(js(j,i),j,i)-6371.2

enddo
enddo



! shift 'jn' & 'js' so that 'jn' starts from '1' on the first field line

npts = maxval(js(:,:)-jn(:,:)+1)

print *, ' npts = ', npts

allocate (lpts(npts,nmp))

do m = 1, nmp
do j = 1, npts
i = j + jn(1,1) - 1

lpts(j,m) = lptx(i,m)

!print *, ' j, lpts(j,m) = ', j, lpts(j,m)

enddo
enddo

print *, ' lpts = ', minval(lpts), maxval(lpts)



! magentic field, arc length, and metric terms

allocate(b(nlp,nmp), s(nlp,nmp))
allocate(h1(nlp,nmp), h2(nlp,nmp), h3(nlp,nmp))

do m = 1, nmp
do l = 1, nlp

klm(l,m) = js(l,m) - jn(l,m) + 1

!print *, ' l, m, klm(l,m) = ', l, m, klm(l,m)

allocate(b(l,m)%vc1d(klm(l,m)), s(l,m)%vc1d(klm(l,m)))
allocate(h1(l,m)%vc1d(klm(l,m)), h2(l,m)%vc1d(klm(l,m)), h3(l,m)%vc1d(klm(l,m)))
enddo
enddo


! 'r' from 'km' to to 'meters'
r = r * 1.0d3


! 'b' in 'tesla'
do m = 1, nmp
do l = 1, nlp
do k = 1, klm(l,m)
kk = k + jn(l,m) - 1

!b(l,m)%vc1d(k) = 7.9e15_wp*sqrt(1.0d0+3.0d0*cos(theta(kk,l,m))**2)/r(kk,l,m)**3

! use the dipole moment for the corresponding epoch
b(l,m)%vc1d(k) = (gm*Re**3)*sqrt(1.0d0+3.0d0*cos(theta(kk,l,m))**2)/r(kk,l,m)**3

enddo

!print *, ' l, m, b = ', l, m, b(l,m)%vc1d(:)

 print *, ' b = ', minval(b(l,m)%vc1d(:)), maxval(b(l,m)%vc1d(:))

enddo
enddo


! metric terms

do m = 1, nmp
do l = 1, nlp

do k = 1, klm(l,m)
kk = k + jn(l,m) - 1

if (use_psi == 0) then
   h1(l,m)%vc1d(k) = r(kk,l,m)**3 /sqrt(1.0d0+3.0d0*cos(theta(kk,l,m))**2) / &
                     (6.3712d6)**2
else 
   h1(l,m)%vc1d(k) = ab*r(kk,l,m)**3*cosh(ab*psi(kk,l,m)) / &
                     (aa*sqrt(1.0d0+3.0d0*cos(theta(kk,l,m))**2)) / &
                     (6.3712d6)**2
endif

h2(l,m)%vc1d(k) = r(kk,l,m)**2 / sqrt(1.0d0+3.0d0*cos(theta(kk,l,m))**2) / &
                  sin(theta(kk,l,m)) / 6.3712d6

h3(l,m)%vc1d(k) = r(kk,l,m) * sin(theta(kk,l,m))

enddo


 print *, ' l, h1 = ', l, minval(h1(l,m)%vc1d(:)), maxval(h1(l,m)%vc1d(:))
!print *, ' l, m, h1 = ', l, m, h1(l,m)%vc1d(:)

 print *, ' l, h2 = ', l, minval(h2(l,m)%vc1d(:)), maxval(h2(l,m)%vc1d(:))
!print *, ' l, m, h2 = ', l, m, h2(l,m)%vc1d(:)

!print *, ' l, m, h2*dchi = ', l, m, h2(l,m)%vc1d(:)*dchi

 print *, ' l, h3 = ', l, minval(h3(l,m)%vc1d(:)), maxval(h3(l,m)%vc1d(:))
!print *, ' l, m, h3 = ', l, m, h3(l,m)%vc1d(:)


enddo
enddo



! 's' in 'meters' along field line (from north footpoint to south footpoint)

do m = 1, nmp
do l = 1, nlp

s(l,m)%vc1d(1) = 0.0d0

do k = 2, klm(l,m)

kk = k + jn(l,m) - 1

if (use_psi == 0) then
   s(l,m)%vc1d(k) = s(l,m)%vc1d(k-1) + &
                    (h1(l,m)%vc1d(k)+h1(l,m)%vc1d(k-1)) * 0.5d0 * dmu
else
   s(l,m)%vc1d(k) = s(l,m)%vc1d(k-1) + &
                    (h1(l,m)%vc1d(k)+h1(l,m)%vc1d(k-1)) * 0.5d0 * dpsi
endif

 enddo


!print *, ' l, m, s = ', l, m, s(l,m)%vc1d(:)
!print *, ' s = ', minval(s(l,m)%vc1d(:))/1000.0, maxval(s(l,m)%vc1d(:))/1000.0


ell(l) = 2.0_wp*dipole_arc_length (r(midpoint,l,m), theta(jn(l,m),l,m))


print *, ' l, ell, s = ', l, ell(l)/1000.0d0, maxval(s(l,m)%vc1d(:))/1000.0d0


mindx = 1.0d36
do k = 2, klm(l,m)
dx = abs(s(l,m)%vc1d(k) - s(l,m)%vc1d(k-1))
if (mindx > dx) mindx = dx
enddo

!print *, ' min s = ', mindx/1.0d3
!print *, ' ds = ', (abs(s(l,m)%vc1d(k) - s(l,m)%vc1d(k-1))/1.0d3, k = 2, klm(l,m))

enddo
enddo


open (51, file="fort.51", status="unknown", form="unformatted")

write(51) ell

close (51)


!dx = 0
!do m = 1, nmp
!do l = 1, nlp 
!dx = dx + klm(l,m)
!enddo
!enddo

!print *, ' total grid points 1 = ', dx

!dx = 0
!do m = 1, nmp
!do k = 1, npts
!dx = dx + lpts(k,m)
!enddo
!enddo

!print *, ' total grid points 2 = ', dx


allocate(x(npts,nmp), y(npts,nmp))

do m = 1, nmp
do k = 1, npts
allocate(x(k,m)%vc1d(lpts(k,m)), y(k,m)%vc1d(lpts(k,m)))
enddo
enddo


! 'x' in 'meters' along 'chi' coordinate (from outer to inner field line)

do m = 1, nmp
do k = 1, npts

x(k,m)%vc1d(1) = 0.0d0

do l = 2, lpts(k,m)

kk = k - (klm(1,m) - klm(l,m))/2


!print *, ' k, l, kk, klm = ', k, l, kk, klm(1,m), klm(l,m)


x(k,m)%vc1d(l) = x(k,m)%vc1d(l-1) + &
                 (h2(l,m)%vc1d(kk) + h2(l-1,m)%vc1d(kk)) * 0.5d0 * dchi

enddo


!print *, ' k, m, x = ', k, m, x(k,m)%vc1d(:)
!print *, ' x = ', minval(x(k,m)%vc1d(:)), maxval(x(k,m)%vc1d(:))


mindx = 1.0d36
do i = 2, lpts(k,m)
dx = abs(x(k,m)%vc1d(i) - x(k,m)%vc1d(i-1))
if (mindx > dx) mindx = dx
enddo

!print *, ' k, mindx = ', k, mindx/1.0d3
!print *, ' dx = ', (abs(x(k,m)%vc1d(i) - x(k,m)%vc1d(i-1))/1.0d3, i = 2, lpts(k,m))


enddo
enddo



allocate(para(nlp,nmp))

do m = 1, nmp
do l = 1, nlp
allocate(para(l,m)%vc1d(klm(l,m)))
enddo
enddo

allocate(perp(npts,nmp))

do m = 1, nmp
do k = 1, npts

allocate(perp(k,m)%vc1d(lpts(k,m)))

perp(k,m)%vc1d(:) = x(k,m)%vc1d(:)

!print *, ' perp 1 = ', perp(k,m)%vc1d(:)

enddo
enddo


! mapping from 'perp' grid to 'para' grid

do m = 1, nmp
do k = 1, npts
do l = 1, lpts(k,m)

kk = k - (klm(1,m) - klm(l,m))/2

para(l,m)%vc1d(kk) = perp(k,m)%vc1d(l)

enddo

!print *, ' para = ', para(l,m)%vc1d(:)

enddo
enddo


! mapping back from 'para' grid to 'perp' grid

do m = 1, nmp
do k = 1, npts
do l = 1, lpts(k,m)

kk = k - (klm(1,m) - klm(l,m))/2

perp(k,m)%vc1d(l) = para(l,m)%vc1d(kk)

enddo

!print *, ' perp 2 = ', perp(k,m)%vc1d(:)

enddo
enddo

print *, ' klm = ', minval(klm), maxval(klm)
print *, ' lpts = ', minval(lpts), maxval(lpts)


! output grid

open (20, file="fort.20", status="unknown", form="unformatted")

write(20) nlp, nmp, npts
write(20) klm, lpts

do m = 1, nmp
do l = 1, nlp
write(20) r(jn(l,m):js(l,m),l,m)
write(20) theta(jn(l,m):js(l,m),l,m)
write(20) b(l,m)%vc1d(:)
write(20) s(l,m)%vc1d(:)
enddo
enddo

close (20)


open (21, file="fort.21", status="unknown", form="unformatted")
do m = 1, nmp
do l = 1, nlp
write(21) h1(l,m)%vc1d(:)
write(21) h2(l,m)%vc1d(:)
write(21) h3(l,m)%vc1d(:)
enddo
enddo
close (21)

open (22, file="fort.22", status="unknown", form="unformatted")
do m = 1, nmp
do k = 1, npts
write(22) x(k,m)%vc1d(:)
enddo
enddo
close (22)

! output mu, chi, phi coordinates
open (24, file="fort.24", status="unknown", form="unformatted")
do m = 1, nmp
do l = 1, nlp
write(24) mu (jn(l,m):js(l,m),l,m)
write(24) chi(jn(l,m):js(l,m),l,m)
write(24) phi(jn(l,m):js(l,m),l,m)
enddo
enddo
close (24)


deallocate (b, s, h1, h2, h3, x, y)


contains



subroutine read_namelist

implicit none

open (unit=11, file='grid_params.namelist')

read (unit=11, nml=grid_params)

close (11)

end subroutine read_namelist



function dipole_arc_length (r_a, theta) result (s)

! dipole arc length from apex

implicit none

real(wp), intent(in) :: r_a, theta

real(wp) :: s

real(wp) :: x

x = cos(theta)

s = r_a/2.0_wp * (x*sqrt(1.0_wp+3.0_wp*x*x) + &
    1.0_wp/sqrt(3.0_wp)*log(sqrt(1.0_wp+3.0_wp*x*x) + sqrt(3.0_wp)*x))

end function dipole_arc_length


end program dipole_coordinates

