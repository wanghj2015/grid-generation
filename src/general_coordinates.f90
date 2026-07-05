program general_coordinates

! Grid generation for a general coordinate system using IGRF
!
! Author: Houjun Wang, wanghoujun@gmail.com, 2020-2021

use mod_precision, only : wp=>dp
use mod_const,     only : pi, dtr
use mod_igrf,      only : get_igrf_coef, gx
use mod_fieldline, only : find_apex_theta, tracing_down, fieldline, &
                          metric, dipole_arc_length, cross_prod
use mod_function,  only : sgn

implicit none


! namelist parameters

integer :: nptx, nlp, nmp, use_psi, use_dipole

real(wp) :: epoch, theta1, theta2, hmin, ds_step, aa

! at r=1 field lines are between colatitudes 'theta1' & 'theta2'
! ds_step: step size [km] for solving field line equations: ds_step = 5.0_wp


! Earth’s mean reference spherical radius [km]
real(wp), parameter :: a = 6371.2_wp


integer :: npts, nptz, midpoint

integer, dimension(:), allocatable :: nnx, nsx, nptw

! magnetic coordinates
real(wp), dimension(:,:,:), allocatable :: mu_m, chi_m, phi_m, psi_m

! geographic coordinates (r, theta, phi)
real(wp), dimension(:,:,:), allocatable :: r, theta, phi

! geographic coordinates at apex
real(wp), dimension(:), allocatable :: chi_apx
real(wp), dimension(:,:), allocatable :: r_apx, theta_apx, phi_apx, Vm_apx


! start (nh) and end (sh) points on each field line
integer, dimension(:,:), allocatable :: jn, js
integer, dimension(:,:), allocatable :: klm, kapx


! end/number of 'l' points on the constant mu/psi curves
!integer, dimension(nptx,nmp) :: lptx
integer, dimension(:,:), allocatable :: lptx

integer, dimension(:,:), allocatable :: lpts

! magnetic field 'b' and arc length 's' along field line
type vc1d
   real(wp), allocatable :: vc1d(:)
end type vc1d

type(vc1d), dimension(:,:), allocatable :: b, s
type(vc1d), dimension(:,:), allocatable :: h1, h2, h3
type(vc1d), dimension(:,:), allocatable :: x, y

type vc2d
   real(wp), allocatable :: vc2d(:,:)
end type vc2d

type(vc2d), dimension(:,:), allocatable :: ct_mu, ct_chi, ct_phi
type(vc2d), dimension(:,:), allocatable :: co_mu, co_chi, co_phi

type vc3d
   real(wp), allocatable :: vc3d(:,:,:)
end type vc3d

type(vc3d), dimension(:,:), allocatable :: ct_g, co_g


! some diagnostic variables
type(vc1d), dimension(:,:), allocatable :: dot12, dot13, dot23
type(vc1d), dimension(:,:), allocatable :: ang12, ang13, ang23
type(vc1d), dimension(:,:), allocatable :: xang12, xang13, xang23
type(vc2d), dimension(:,:), allocatable :: crs23
type(vc1d), dimension(:,:), allocatable :: dotxx, angxx, b1, b2

real(wp) :: a1, a2, a3


real(wp) :: r1, r2, chi1, chi2, dchi

real(wp) :: mu1, mu2, dmu, psi1, psi2, dpsi

real(wp) :: ab, g10, gm

real(wp) :: ell, ds

real(wp), dimension(:), allocatable :: dmux, dpsix, ellx

real(wp) :: r_n, theta_n, phi_n, Vm_n, ell_n, r_s, theta_s, phi_s, Vm_s, ell_s

real(wp) :: r0, theta0, phi0, rx, thetax, phix

real(wp), dimension(:), allocatable :: wrk1, wrk2, wrk3, wrk4

integer :: i, j, k, l, m, nn, ns, kk, Kp

namelist /grid_params/ epoch, nptx, nlp, nmp, theta1, theta2, hmin, &
                       ds_step, use_psi, aa, use_dipole


! get namelist
call read_namelist


! allocate space
allocate (jn(nlp,nmp), js(nlp,nmp))
allocate (klm(nlp,nmp), kapx(nlp,nmp))

allocate (nnx(nmp), nsx(nmp), nptw(nmp))
allocate (dmux(nmp), dpsix(nmp), ellx(nmp))

allocate (chi_apx(nlp))
allocate (r_apx(nlp,nmp), theta_apx(nlp,nmp), phi_apx(nlp,nmp))
allocate (Vm_apx(nlp,nmp))


! get the IGRF coefficients for the epoch
call get_igrf_coef (epoch, use_dipole)


r1 = a/sin(theta1*dtr)**2
r2 = a/sin(theta2*dtr)**2

!print *, ' r1, r2 = ', r1-a, r2-a

! const-chi curves in a meridian plane (phi=const), denoting dipole field lines

chi1 = 1.0d0/r1
chi2 = 1.0d0/r2
dchi = (chi2-chi1)/(nlp-1)

chi_apx(1) = chi1
do l = 2, nlp
chi_apx(l) = chi_apx(l-1) + dchi
enddo

! find geographic coordinates at apex: (r_apx, theta_apx, phi_apx)

do m = 1, nmp

do l = 1, nlp

r_apx(l,m)   = 1.0d0/chi_apx(l)

!print *, ' r_apx = ', r_apx - a

 phi_apx(l,m) = (m-1) * (2.0d0*pi/nmp)
!phi_apx(l,m) = pi * 0.5d0
!phi_apx(l,m) = pi
!phi_apx(l,m) = pi * 1.5d0


! initial guess for co-latitude
theta_apx(l,m) = pi/2.0d0

call find_apex_theta (r_apx(l,m), phi_apx(l,m), theta_apx(l,m), Vm_apx(l,m))

!print *, ' l, theta_apx, Vm_apx = ', l, theta_apx(l,m)/dtr, Vm_apx(l,m)

enddo ! do l = 1, nlp


! from apex tracing down (north and south) to locate (r, theta, phi) 
! where the field line cross the Earth's surface 

! apex location of the outer-most field line
r0     = r_apx(1,m)
theta0 = theta_apx(1,m)
phi0   = phi_apx(1,m)

! estimate field line arc length from a dipole

ell = dipole_arc_length (r_apx(1,m), theta1*dtr)

ell = 1.5d0 * ell

!print *, ' arc length = ', ell

! step size [km]
ds = ds_step

! to the north

sgn = +1.0d0

call tracing_down (r0, theta0, phi0, ell, ds, &
                   r_n, theta_n, phi_n, Vm_n, ell_n)

!print *, ' r_n, theta_n, phi_n, Vm_n, ell_n = ', &
!           r_n, theta_n/dtr, phi_n/dtr, Vm_n, ell_n

! to the south

sgn = -1.0d0

call tracing_down (r0, theta0, phi0, ell, ds, &
                   r_s, theta_s, phi_s, Vm_s, ell_s)

!print *, ' r_s, theta_s, phi_s, Vm_s, ell_s = ', &
!           r_s, theta_s/dtr, phi_s/dtr, Vm_s, ell_s

! updated arc length
ell = max(ell_n, ell_s)

print *, ' updated arc length = ', ell

ellx(m) = ell


if (use_psi == 0) then

   ! divide the potential equally along the field line

   mu1 = Vm_n; mu2 = Vm_s
   dmu = (mu2-mu1)/(nptx-1)

   print *, ' dmu = ', dmu

   dmux(m) = dmu

   mu1 = Vm_n; mu2 = Vm_apx(1,m)
   nn = ceiling((mu2-mu1)/dmu + 1.0d-20)

   mu1 = Vm_s; mu2 = Vm_apx(1,m)
   ns = ceiling((mu1-mu2)/dmu + 1.0d-20)

   nptz = 2 * max(nn, ns) - 1

   print *, ' nn, ns, nptz = ', nn, ns, nptz

   midpoint = (nptz+1)/2


   nnx(m) = nn
   nsx(m) = ns
   nptw(m) = nptz

else

   ! 'psi' grid

   ! mapping mu_m to psi_m

   ! dipole moment for the corresponding epoch 
   g10 = gx(1,0)
   gm  = abs(a*g10)


   ab = log(aa + sqrt(1.0d0 + aa**2))

  !mu1 = Vm_n/gm; mu2 = Vm_s/gm
   mu1 = (Vm_n - Vm_apx(1,m))/gm
   mu2 = (Vm_s - Vm_apx(1,m))/gm

   psi1 = log(aa*mu1 + sqrt(1.0d0 + (aa*mu1)**2))/ab
   psi2 = log(aa*mu2 + sqrt(1.0d0 + (aa*mu2)**2))/ab

   dpsi = (psi2-psi1)/(nptx-1)

   print *, ' aa, ab = ', aa, ab
   print *, ' mu1, mu2 = ', mu1, mu2
   print *, ' psi1, psi2, dpsi = ', psi1, psi2, dpsi

   dpsix(m) = dpsi

  !mu1 = Vm_n/gm; mu2 = Vm_apx(1,m)/gm
   mu1 = (Vm_n - Vm_apx(1,m))/gm
  !mu2 = (Vm_apx(1,m) - Vm_apx(1,m))/gm
   mu2 = 0.0d0

   psi1 = log(aa*mu1 + sqrt(1.0d0 + (aa*mu1)**2))/ab
   psi2 = log(aa*mu2 + sqrt(1.0d0 + (aa*mu2)**2))/ab

   print *, ' psi1, psi2 = ', psi1, psi2

   nn = ceiling((psi2-psi1)/dpsi + 1.0d-20)

  !mu1 = Vm_s/gm; mu2 = Vm_apx(1,m)/gm
   mu1 = (Vm_s - Vm_apx(1,m))/gm 
  !mu2 = (Vm_apx(1,m) - Vm_apx(1,m))/gm
   mu2 = 0.0d0

   psi1 = log(aa*mu1 + sqrt(1.0d0 + (aa*mu1)**2))/ab
   psi2 = log(aa*mu2 + sqrt(1.0d0 + (aa*mu2)**2))/ab

   print *, ' psi1, psi2 = ', psi1, psi2

   ns = ceiling((psi1-psi2)/dpsi + 1.0d-20)

   nptz = 2 * max(nn, ns) - 1

   print *, ' nn, ns, nptz = ', nn, ns, nptz

   midpoint = (nptz+1)/2

   print *, ' midpoint = ', midpoint

   nnx(m) = nn
   nsx(m) = ns
   nptw(m) = nptz

endif

enddo ! m = 1, nmp


! allocate space

nptz = maxval(nptw)

allocate (mu_m(nptz,nlp,nmp), chi_m(nptz,nlp,nmp), phi_m(nptz,nlp,nmp))
allocate (psi_m(nptz,nlp,nmp))


if (use_psi == 0) then

   ! 'potential' equally distributed along the field line

   do i = 1, nmp
   midpoint = (nptw(i)+1)/2
   dmu = dmux(i)
   do j = 1, nlp
   ! at apex
   mu_m(midpoint,j,i) = Vm_apx(j,i)
   do k = midpoint-1, 1, -1
   mu_m(k,j,i) = mu_m(k+1,j,i) - dmu
   enddo
   do k = midpoint+1, nptz
   mu_m(k,j,i) = mu_m(k-1,j,i) + dmu
   enddo
   enddo
   enddo

  !print *, ' mu = ', mu_m(midpoint,:,1)
  !print *, ' mu = ', mu_m(:,1,1) / gm
  !print *, ' mu = ', mu_m(:,1,1)


else

   ! 'psi' grid

   do i = 1, nmp
   midpoint = (nptw(i)+1)/2
   dpsi = dpsix(i)
   do j = 1, nlp

   ! at apex
  !mu2 = Vm_apx(j,i)/gm
   mu2 = 0.0d0

   psi2 = log(aa*mu2 + sqrt(1.0d0 + (aa*mu2)**2))/ab
   psi_m(midpoint,j,i) = psi2

   do k = midpoint-1, 1, -1
   psi_m(k,j,i) = psi_m(k+1,j,i) - dpsi
   enddo
   do k = midpoint+1, nptw(i)
   psi_m(k,j,i) = psi_m(k-1,j,i) + dpsi
   enddo
   enddo
   enddo

   ! convert 'psi' to 'mu'
   do i = 1, nmp
   do j = 1, nlp
   do k = 1, nptw(i)
   mu_m(k,j,i) = sinh(ab*psi_m(k,j,i))/aa

   ! get the potential at the grid points,
   ! which is used in the following computation
   mu_m(k,j,i) = mu_m(k,j,i) * gm + Vm_apx(j,i)
   enddo
   enddo
   enddo

  !mu_m = mu_m * gm

  !print *, ' mu = ', mu_m(midpoint,:,1)
  !print *, ' mu = ', mu_m(:,1,1)

endif


! set magnetic phi_m to phi at apex

do i = 1, nmp
do j = 1, nlp
!do k = 1, nptz
do k = 1, nptw(i)
!phi_m(k,j,i) = (i-1) * (2.0d0*pi/nmp)
phi_m(k,j,i) = phi_apx(j,i)
enddo
enddo
enddo

! const-chi curves in a meridian plane (phi=const), denoting dipole field lines

do i = 1, nmp
!do k = 1, nptz
do k = 1, nptw(i)
chi_m(k,1,i) = chi1
do j = 2, nlp
chi_m(k,j,i) = chi_m(k,j-1,i) + dchi
enddo
enddo
enddo

! from (mu_m, chi_m, phi_m) find (r, theta, phi) tracing field line
! with fixed dmu, starting from apex until below 90km altitude

allocate (r   (nptz,nlp,nmp), theta(nptz,nlp,nmp), phi  (nptz,nlp,nmp))
!allocate (wrk1(midpoint), wrk2(midpoint), wrk3(midpoint), wrk4(midpoint))

ds = ds_step

do m = 1, nmp

midpoint = (nptw(m)+1)/2

ell = ellx(m)

do l = 1, nlp

!print *, ' l = ', l, ' to the north'

! to the north
sgn = +1.0d0

nn = midpoint

allocate (wrk1(nn), wrk2(nn), wrk3(nn), wrk4(nn))

wrk1(:) = mu_m(nn:1:-1,l,m)

call fieldline (r_apx(l,m), theta_apx(l,m), phi_apx(l,m), ell, ds, &
                wrk1, wrk2, wrk3, wrk4)

r    (1:nn,l,m) = wrk2(nn:1:-1)
theta(1:nn,l,m) = wrk3(nn:1:-1)
phi  (1:nn,l,m) = wrk4(nn:1:-1)

deallocate (wrk1, wrk2, wrk3, wrk4)

!print *, ' l = ', l, ' to the south'

! to the south
sgn = -1.0d0

ns = nptw(m)-midpoint+1

allocate (wrk1(ns), wrk2(ns), wrk3(ns), wrk4(ns))

wrk1(:) = mu_m(nn:nn+ns-1,l,m)

call fieldline (r_apx(l,m), theta_apx(l,m), phi_apx(l,m), ell, ds, &
                wrk1, wrk2, wrk3, wrk4)

r    (nn:nn+ns-1,l,m) = wrk2(1:ns)
theta(nn:nn+ns-1,l,m) = wrk3(1:ns)
phi  (nn:nn+ns-1,l,m) = wrk4(1:ns)

deallocate (wrk1, wrk2, wrk3, wrk4)

!print *, ' m, l, r     = ', m, l, pack(r(:,l,m)-a, r(:,l,m) /= 0)
!print *, ' m, l, theta = ', m, l, pack(theta(:,l,m)/dtr, theta(:,l,m) /= 0)
!print *, ' m, l, phi   = ', m, l, pack(phi(:,l,m)/dtr, phi(:,l,m) /= 0)

enddo
enddo


! tidy up the grids


allocate (lptx(nptz,nmp))

lptx = 0

do i = 1, nmp

do k = 1, nptw(i)

do j = nlp, 1, -1
if (r(k,j,i) > a + hmin) then
   lptx(k,i) = j

   !print *, ' k, lptx = ', k, lptx(k,i), r(k,lptx(k,i),i)-a

   exit
endif
enddo

!print *, ' m, k, lptx = ', i, k, lptx(k,i), r(k,lptx(k,i),i)-a

enddo

enddo



do i = 1, nmp

midpoint = (nptw(i)+1)/2

do j = 1, nlp

do k = 1, midpoint
!if (r(k,j,i) > a + hmin) then
if (r(k,j,i) > a + hmin .and. lptx(k,i) > 1) then
   jn(j,i) = k
   exit
endif
enddo

!do k = nptz, midpoint, -1
do k = nptw(i), midpoint, -1
!if (r(k,j,i) > a + hmin) then
if (r(k,j,i) > a + hmin .and. lptx(k,i) > 1) then
   js(j,i) = k
   exit
endif
enddo


! adjust grid
!if (i == 3 .or. i == 4) then
!   jn(1,i) = jn(1,i) + 1
!   js(1,i) = js(1,i) - 1
!endif


! location of apex point
kapx(j,i) = midpoint - jn(j,i) + 1

!print *, ' m, l, jn, js = ', i, j, jn(j,i), js(j,i), &
!           kapx(j,i)
!           r(jn(j,i),j,i)-a, r(js(j,i),j,i)-a

enddo ! nlp

enddo ! nmp



! shift 'jn' & 'js' so that 'jn' starts from '1' on the first field line

do m = 1, nmp
!print *, ' m, npts = ', m, maxval(js(:,m)-jn(:,m)+1)
print *, ' m, npts = ', m, (js(1,m)-jn(1,m)+1)
print *, ' m, lptx = ', m, lptx(jn(1,m),m), lptx(js(1,m),m)
enddo


!npts = maxval(js(:,:)-jn(:,:)+1)
!npts = maxval(js(1,:)-jn(1,:)+1)
npts = minval(js(1,:)-jn(1,:)+1)

print *, ' npts = ', npts

! reset js
js(1,:) = npts + jn(1,:) - 1


allocate (lpts(npts,nmp))

do m = 1, nmp

do j = 1, npts

i = j+jn(1,m)-1

lpts(j,m) = lptx(i,m)

!print *, ' m, j, lpts(j,m) = ', m, j, lpts(j,m)

enddo

 print *, ' m, lpts = ', m, minval(lpts(:,m), mask=lpts(:,m)>0), maxval(lpts(:,m))
!print *, ' m, lpts = ', m, minval(lpts(:,m)                  ), maxval(lpts(:,m))

enddo

!print *, ' lpts = ', minval(lpts), maxval(lpts)


! magentic field, arc length, and metric terms

! allocate space

allocate(b(nlp,nmp), s(nlp,nmp))
allocate(ct_mu(nlp,nmp), ct_chi(nlp,nmp), ct_phi(nlp,nmp))
allocate(co_mu(nlp,nmp), co_chi(nlp,nmp), co_phi(nlp,nmp))
allocate(h1(nlp,nmp), h2(nlp,nmp), h3(nlp,nmp))
allocate(ct_g(nlp,nmp), co_g(nlp,nmp))

do m = 1, nmp
do l = 1, nlp

klm(l,m) = js(l,m) - jn(l,m) + 1

!print *, ' l, m, klm(l,m) = ', l, m, klm(l,m)

Kp = klm(l,m)

allocate(b(l,m)%vc1d(Kp), s(l,m)%vc1d(Kp))
allocate(h1(l,m)%vc1d(Kp), h2(l,m)%vc1d(Kp), h3(l,m)%vc1d(Kp))
allocate(ct_mu(l,m)%vc2d(3,Kp), ct_chi(l,m)%vc2d(3,Kp), ct_phi(l,m)%vc2d(3,Kp))
allocate(co_mu(l,m)%vc2d(3,Kp), co_chi(l,m)%vc2d(3,Kp), co_phi(l,m)%vc2d(3,Kp))
allocate(ct_g(l,m)%vc3d(3,3,Kp), co_g(l,m)%vc3d(3,3,Kp))

enddo
enddo

! magentic field, basis vectors, arc length, and metric terms

ell = 1.2d0 * ell
ds  = ds_step

do m = 1, nmp
do l = 1, nlp

Kp = klm(l,m)

allocate (wrk1(Kp), wrk2(Kp), wrk3(Kp))

wrk1(:) = r    (jn(l,m):js(l,m),l,m)
wrk2(:) = theta(jn(l,m):js(l,m),l,m)
wrk3(:) = phi  (jn(l,m):js(l,m),l,m)

call metric(wrk1, wrk2, wrk3, ell, ds, b(l,m)%vc1d(:),                   &
     ct_mu(l,m)%vc2d(:,:), ct_chi(l,m)%vc2d(:,:), ct_phi(l,m)%vc2d(:,:), &
     co_mu(l,m)%vc2d(:,:), co_chi(l,m)%vc2d(:,:), co_phi(l,m)%vc2d(:,:), &
     ct_g(l,m)%vc3d(:,:,:), co_g(l,m)%vc3d(:,:,:),                       &
     h1(l,m)%vc1d(:), h2(l,m)%vc1d(:), h3(l,m)%vc1d(:))
 
  
 print *, ' l, b = ', l, minval(b(l,m)%vc1d(:)*1.0d-9), maxval(b(l,m)%vc1d(:)*1.0d-9)
!print *, ' l, m, b = ', l, m, b(l,m)%vc1d(:)


if (use_psi == 1) then
   do k = 1, klm(l,m)
   kk = k + jn(l,m) - 1
   h1(l,m)%vc1d(k) = h1(l,m)%vc1d(k) * ab*cosh(ab*psi_m(kk,l,m)) / aa
   enddo
endif


 g10 = gx(1,0)
 gm  = abs(a*g10)

!print *, ' h1 = ', minval(h1(l,m)%vc1d(:)*a**2), maxval(h1(l,m)%vc1d(:)*a**2)
 print *, ' l, h1 = ', l, minval(h1(l,m)%vc1d(:)*gm*1000.0d0), &
                          maxval(h1(l,m)%vc1d(:)*gm*1000.0d0)
!print *, ' l, m, h1 = ', l, m, h1(l,m)%vc1d(:)

 print *, ' l, h2 = ', l, minval(h2(l,m)%vc1d(:)/a*1.0d3), maxval(h2(l,m)%vc1d(:)/a*1.0d3)
!print *, ' l, m, h2 = ', l, m, h2(l,m)%vc1d(:)

!print *, ' l, m, h2*dchi = ', l, m, h2(l,m)%vc1d(:)*dchi

 print *, ' l, h3 = ', l, minval(h3(l,m)%vc1d(:)*1.0d3), maxval(h3(l,m)%vc1d(:)*1.0d3)
!print *, ' l, m, h3 = ', l, m, h3(l,m)%vc1d(:)

deallocate (wrk1, wrk2, wrk3)
enddo
enddo

! arc length distribution along field line

do m = 1, nmp
do l = 1, nlp
s(l,m)%vc1d(1) = 0.0d0
do k = 2, klm(l,m)
kk = k + jn(l,m) - 1

if (use_psi == 0) then
   s(l,m)%vc1d(k) = s(l,m)%vc1d(k-1) + &
                    (h1(l,m)%vc1d(k)+h1(l,m)%vc1d(k-1)) * 0.5d0 * &
                    (mu_m(kk,l,m)-mu_m(kk-1,l,m)) 
else
   s(l,m)%vc1d(k) = s(l,m)%vc1d(k-1) + &
                    (h1(l,m)%vc1d(k)+h1(l,m)%vc1d(k-1)) * 0.5d0 * &
                    (psi_m(kk,l,m)-psi_m(kk-1,l,m)) * gm
endif

enddo

!print *, ' l, m, s = ', l, m, s(l,m)%vc1d(:)
 print *, ' s = ', minval(s(l,m)%vc1d(:)), maxval(s(l,m)%vc1d(:))

enddo
enddo


! 'x' distance along 'chi' coordinate (from outer to inner field line)

allocate(x(npts,nmp), y(npts,nmp))

do m = 1, nmp
do k = 1, npts
allocate(x(k,m)%vc1d(lpts(k,m)), y(k,m)%vc1d(lpts(k,m)))
enddo
enddo

do m = 1, nmp
do k = 1, npts

x(k,m)%vc1d(1) = 0.0

do l = 2, lpts(k,m)

!kk = k - (klm(1,m) - klm(l,m))/2
kk = k - (kapx(1,m) - kapx(l,m))

!print *, ' k, l, kk, klm = ', k, l, kk, klm(1,m), klm(l,m)

x(k,m)%vc1d(l) = x(k,m)%vc1d(l-1) + &
                 (h2(l,m)%vc1d(kk) + h2(l-1,m)%vc1d(kk)) * 0.5d0 * dchi
enddo

!print *, ' k, m, x = ', k, m, x(k,m)%vc1d(:)
!print *, ' x = ', minval(x(k,m)%vc1d(:)), maxval(x(k,m)%vc1d(:))

enddo
enddo



! some diagnostics

allocate(ang12(nlp,nmp), ang13(nlp,nmp), ang23(nlp,nmp))
allocate(xang12(nlp,nmp), xang13(nlp,nmp), xang23(nlp,nmp))
allocate(dot12(nlp,nmp), dot13(nlp,nmp), dot23(nlp,nmp), crs23(nlp,nmp))
allocate(dotxx(nlp,nmp), angxx(nlp,nmp), b1   (nlp,nmp), b2   (nlp,nmp))

do m = 1, nmp
do l = 1, nlp

Kp = klm(l,m)

allocate(ang12(l,m)%vc1d(Kp), ang13(l,m)%vc1d(Kp), ang23(l,m)%vc1d(Kp))
allocate(xang12(l,m)%vc1d(Kp), xang13(l,m)%vc1d(Kp), xang23(l,m)%vc1d(Kp))
allocate(dot12(l,m)%vc1d(Kp), dot13(l,m)%vc1d(Kp), dot23(l,m)%vc1d(Kp))
allocate(crs23(l,m)%vc2d(3,Kp))
allocate(dotxx(l,m)%vc1d(Kp), angxx(l,m)%vc1d(Kp))
allocate(b1   (l,m)%vc1d(Kp), b2   (l,m)%vc1d(Kp))

enddo
enddo


do m = 1, nmp
do l = 1, nlp

Kp = klm(l,m)


! are they orthogonal?

do k = 1, Kp

dot12(l,m)%vc1d(k) = dot_product(ct_mu (l,m)%vc2d(:,k), &
                                 ct_chi(l,m)%vc2d(:,k))
dot13(l,m)%vc1d(k) = dot_product(ct_mu (l,m)%vc2d(:,k), &
                                 ct_phi(l,m)%vc2d(:,k))
dot23(l,m)%vc1d(k) = dot_product(ct_chi(l,m)%vc2d(:,k), &
                                 ct_phi(l,m)%vc2d(:,k))

a1 = sqrt(dot_product(ct_mu (l,m)%vc2d(:,k), ct_mu (l,m)%vc2d(:,k)))
a2 = sqrt(dot_product(ct_chi(l,m)%vc2d(:,k), ct_chi(l,m)%vc2d(:,k)))
a3 = sqrt(dot_product(ct_phi(l,m)%vc2d(:,k), ct_phi(l,m)%vc2d(:,k)))

ang12(l,m)%vc1d(k) = acos(dot12(l,m)%vc1d(k) / (a1 * a2)) / dtr 
ang13(l,m)%vc1d(k) = acos(dot13(l,m)%vc1d(k) / (a1 * a3)) / dtr 
ang23(l,m)%vc1d(k) = acos(dot23(l,m)%vc1d(k) / (a2 * a3)) / dtr 



! use covariant basis vectors

dot12(l,m)%vc1d(k) = dot_product(co_mu (l,m)%vc2d(:,k), &
                                 co_chi(l,m)%vc2d(:,k))
dot13(l,m)%vc1d(k) = dot_product(co_mu (l,m)%vc2d(:,k), &
                                 co_phi(l,m)%vc2d(:,k))
dot23(l,m)%vc1d(k) = dot_product(co_chi(l,m)%vc2d(:,k), &
                                 co_phi(l,m)%vc2d(:,k))

a1 = sqrt(dot_product(co_mu (l,m)%vc2d(:,k), co_mu (l,m)%vc2d(:,k)))
a2 = sqrt(dot_product(co_chi(l,m)%vc2d(:,k), co_chi(l,m)%vc2d(:,k)))
a3 = sqrt(dot_product(co_phi(l,m)%vc2d(:,k), co_phi(l,m)%vc2d(:,k)))

xang12(l,m)%vc1d(k) = acos(dot12(l,m)%vc1d(k) / (a1 * a2)) / dtr
xang13(l,m)%vc1d(k) = acos(dot13(l,m)%vc1d(k) / (a1 * a3)) / dtr
xang23(l,m)%vc1d(k) = acos(dot23(l,m)%vc1d(k) / (a2 * a3)) / dtr

enddo


print *, ' l, dot12 = ', l, minval(dot12(l,m)%vc1d(:)), &
                            maxval(dot12(l,m)%vc1d(:))
print *, ' l, dot13 = ', l, minval(dot13(l,m)%vc1d(:)), &
                            maxval(dot13(l,m)%vc1d(:))
print *, ' l, dot23 = ', l, minval(dot23(l,m)%vc1d(:)), &
                            maxval(dot23(l,m)%vc1d(:))

print *, ' l, ang12 = ', l, minval(ang12(l,m)%vc1d(:)), &
                            maxval(ang12(l,m)%vc1d(:))
print *, ' l, ang13 = ', l, minval(ang13(l,m)%vc1d(:)), &
                            maxval(ang13(l,m)%vc1d(:))
print *, ' l, ang23 = ', l, minval(ang23(l,m)%vc1d(:)), &
                            maxval(ang23(l,m)%vc1d(:))


print *, ' l, xang12 = ', l, minval(xang12(l,m)%vc1d(:)), &
                             maxval(xang12(l,m)%vc1d(:))
print *, ' l, xang13 = ', l, minval(xang13(l,m)%vc1d(:)), &
                             maxval(xang13(l,m)%vc1d(:))
print *, ' l, xang23 = ', l, minval(xang23(l,m)%vc1d(:)), &
                             maxval(xang23(l,m)%vc1d(:))


! are they Euler potentials?


! get the l-dependent normalization constant

crs23(l,m)%vc2d(:,1) = cross_prod(ct_chi(l,m)%vc2d(:,1), &
                                  ct_phi(l,m)%vc2d(:,1))
a3 = sqrt(dot_product(crs23(l,m)%vc2d(:,1), crs23(l,m)%vc2d(:,1))) / &
     b(l,m)%vc1d(1)


do k = 1, Kp

crs23(l,m)%vc2d(:,k) = cross_prod(ct_chi(l,m)%vc2d(:,k), &
                                  ct_phi(l,m)%vc2d(:,k))

 crs23(l,m)%vc2d(:,k) = crs23(l,m)%vc2d(:,k) * 6371.2d0**3 * 29619.4d0
!crs23(l,m)%vc2d(:,k) = crs23(l,m)%vc2d(:,k) / a3


dotxx(l,m)%vc1d(k) = dot_product(ct_mu(l,m)%vc2d(:,k), crs23(l,m)%vc2d(:,k))

a1 = sqrt(dot_product(ct_mu(l,m)%vc2d(:,k), ct_mu(l,m)%vc2d(:,k)))
a2 = sqrt(dot_product(crs23(l,m)%vc2d(:,k), crs23(l,m)%vc2d(:,k)))

angxx(l,m)%vc1d(k) = acos(dotxx(l,m)%vc1d(k) / (a1 * a2)) / dtr

b1(l,m)%vc1d(k) = a1
b2(l,m)%vc1d(k) = a2


enddo


print *, ' l, crs23 = ', l, minval(crs23(l,m)%vc2d(:,:)), & 
                            maxval(crs23(l,m)%vc2d(:,:))

print *, ' l, ct_mu = ', l, minval(ct_mu(l,m)%vc2d(:,:)), & 
                            maxval(ct_mu(l,m)%vc2d(:,:))

print *, ' l, diffx = ', l, minval((crs23(l,m)%vc2d(:,:) - &
                                    ct_mu(l,m)%vc2d(:,:))), &
                            maxval((crs23(l,m)%vc2d(:,:) - &
                                    ct_mu(l,m)%vc2d(:,:)))

 print *, ' l, r_err = ', l, minval((crs23(l,m)%vc2d(:,:) - &
                                     ct_mu(l,m)%vc2d(:,:)) / &
                                     ct_mu(l,m)%vc2d(:,:)), &
                             maxval((crs23(l,m)%vc2d(:,:) - &
                                     ct_mu(l,m)%vc2d(:,:)) / &
                                     ct_mu(l,m)%vc2d(:,:))

 print *, ' l, angxx = ', l, minval(angxx(l,m)%vc1d(:)), &
                             maxval(angxx(l,m)%vc1d(:))

 print *, ' l, ratio = ', l, minval(b1(l,m)%vc1d(:)/b2(l,m)%vc1d(:)), &
                             maxval(b1(l,m)%vc1d(:)/b2(l,m)%vc1d(:))

 print *, ' l, ratio2= ', l, minval(b1(l,m)%vc1d(:)/b (l,m)%vc1d(:)), &
                             maxval(b1(l,m)%vc1d(:)/b (l,m)%vc1d(:))

enddo
enddo




! output grid

open (20, file="fort.20", status="unknown", form="unformatted")

write(20) nlp, nmp, npts
write(20) klm, lpts

do m = 1, nmp
do l = 1, nlp
write(20) r    (jn(l,m):js(l,m),l,m) * 1.0d3
write(20) theta(jn(l,m):js(l,m),l,m)
write(20) phi  (jn(l,m):js(l,m),l,m)
write(20) b(l,m)%vc1d(:) * 1.0d-9
write(20) s(l,m)%vc1d(:) * 1.0d3
enddo
enddo

close (20)


open (21, file="fort.21", status="unknown", form="unformatted")

do m = 1, nmp
do l = 1, nlp
write(21) h1(l,m)%vc1d(:) * 1.0d3
write(21) h2(l,m)%vc1d(:) * 1.0d3
write(21) h3(l,m)%vc1d(:) * 1.0d3
enddo
enddo
close (21)

open (22, file="fort.22", status="unknown", form="unformatted")

do m = 1, nmp
do k = 1, npts
write(22) x(k,m)%vc1d(:) * 1.0d3
enddo
enddo

write(22) kapx

close (22)

open (23, file="fort.23", status="unknown", form="unformatted")

do m = 1, nmp
do l = 1, nlp
write(23) ct_mu (l,m)%vc2d(:,:) / 1.0d3
write(23) ct_chi(l,m)%vc2d(:,:) / 1.0d3
write(23) ct_phi(l,m)%vc2d(:,:) / 1.0d3
write(23) co_mu (l,m)%vc2d(:,:) * 1.0d3
write(23) co_chi(l,m)%vc2d(:,:) * 1.0d3
write(23) co_phi(l,m)%vc2d(:,:) * 1.0d3
enddo
enddo

close (23)


do m = 1, nmp
do l = 1, nlp
do k = jn(l,m), js(l,m)
mu_m(k,l,m) = (mu_m(k,l,m) - Vm_apx(l,m)) / gm
enddo
enddo
enddo


open (24, file="fort.24", status="unknown", form="unformatted")

do m = 1, nmp
do l = 1, nlp
write(24) mu_m (jn(l,m):js(l,m),l,m)
write(24) chi_m(jn(l,m):js(l,m),l,m)
write(24) phi_m(jn(l,m):js(l,m),l,m)
enddo
enddo

close (24)


! output some diagnostic variables

open (25, file="fort.25", status="unknown", form="unformatted")

do m = 1, nmp
do l = 1, nlp
write(25) ang12(l,m)%vc1d(:)
write(25) ang13(l,m)%vc1d(:)
write(25) ang23(l,m)%vc1d(:)
write(25) b1(l,m)%vc1d(:)/b2(l,m)%vc1d(:)
enddo
enddo

close (25)

open (26, file="fort.26", status="unknown", form="unformatted")

do m = 1, nmp
do l = 1, nlp
write(26) xang12(l,m)%vc1d(:)
write(26) xang13(l,m)%vc1d(:)
write(26) xang23(l,m)%vc1d(:)
enddo
enddo

close (26)


! coordinate transformation matrix

!do m = 1, nmp
!do l = 1, nlp
!print *, ' s2d1 = ', 2.0*cos(theta(jn(l,m):js(l,m),l,m)) / &
!                     sqrt(1.0+3.0*cos(theta(jn(l,m):js(l,m),l,m))**2)

!print *, ' co_mu1 ', co_mu(l,m)%vc2d(1,:)/(h1(l,m)%vc1d(:))

!print *, ' s2d2 = ', sin(theta(jn(l,m):js(l,m),l,m)) / &
!                     sqrt(1.0+3.0*cos(theta(jn(l,m):js(l,m),l,m))**2)
!print *, ' co_mu2 ', co_mu(l,m)%vc2d(2,:)/(h1(l,m)%vc1d(:))

!print *, ' co_mu3 ', co_mu(l,m)%vc2d(3,:)/(h1(l,m)%vc1d(:))

!enddo
!enddo


! check B/B0

!do m = 1, nmp
!do l = 1, nlp
!print *, ' b/gm = ', b(l,m)%vc1d(:)/abs(g10)
!print *, ' b/b0 = ', sqrt(1.0+3.0*cos(theta(jn(l,m):js(l,m),l,m))**2) / &
!                     (r(jn(l,m):js(l,m),l,m)/a)**3
!enddo
!enddo



deallocate (jn, js, klm, kapx)
deallocate (nnx, nsx, nptw, dmux, dpsix, ellx)
deallocate (chi_apx, r_apx, theta_apx, phi_apx, Vm_apx)

deallocate (b, s)
deallocate (ct_mu, ct_chi, ct_phi)
deallocate (co_mu, co_chi, co_phi)
deallocate (h1, h2, h3, ct_g, co_g)
deallocate (mu_m, chi_m, phi_m, psi_m, r, theta, phi)
deallocate (lptx, lpts)
deallocate (x, y)
deallocate (ang12, ang13, ang23)
deallocate (dot12, dot13, dot23, crs23)
deallocate (dotxx, angxx, b1, b2)
deallocate (xang12, xang13, xang23)



contains



subroutine read_namelist

implicit none

open (unit=11, file='grid_params.namelist')

read (unit=11, nml=grid_params)

close (11)

end subroutine read_namelist


end program general_coordinates


