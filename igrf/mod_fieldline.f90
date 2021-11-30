module mod_fieldline

use mod_precision, only : wp=>dp
use mod_const,     only : pi, dtr
use mod_igrf,      only : igrf

implicit none

! Earth’s mean reference spherical radius [km] 
real(wp), parameter :: a = 6371.2_wp

!real(wp), parameter :: pi = 3.14159265358979_wp

real(wp), parameter :: tol = 1.0d-12


contains


subroutine find_apex_theta (r_a, phi_a, theta_a, Vm_a)

! given geographic coordinates at apex (r_a, phi_a) and a guess of 'theta_a', 
! use newton-raphson method to find theta_a (theta at apex)

implicit none

integer, parameter :: kmax = 100

!real(wp), parameter :: tol = 1.0-15

real(wp), intent(in) :: r_a, phi_a

real(wp), intent(inout) :: theta_a, Vm_a

real(wp) :: x, dx, fx, dfx

real(wp) :: Br, Bt, Bp, Vm, dBrdt

integer :: k, ierror
 
! initial guess
x = theta_a
 
do k = 1, kmax

   call igrf (r_a, theta_a, phi_a, Br, Bt, Bp, Vm, dBrdt)

   !print *, ' Br, Bt, Bp, Vm, dBrdt = ', Br, Bt, Bp, Vm, dBrdt

   Vm_a = Vm

   fx  = Br
   dfx = dBrdt
 
   ! if the error tolerance is satisfied, then jump out of do loop
   if ( abs (fx) < tol ) exit

   if ( k > kmax .and. abs (fx) > tol ) then
      print *, '*** Warning: has not yet converged'
   endif

   if ( dfx == 0.0d+00 ) then
      print *, '*** Warning: dfx too small'
      exit
   endif
 
   ! compute increment
   dx = - fx / dfx
 
   ! update x
   x = x + dx
    
   theta_a = x

enddo

end subroutine find_apex_theta



subroutine tracing_down (r0, theta0, phi0, s, ds, rx, thetax, phix, Vmx, ell)

! tracing down to locate where field line cross Earth's aurface

! Input:
! (r0, theta0, phi0): starting point (e.g., apex)
! s: estimate of arc legnth (based on dipole field)
! ds: step size using for tracing field line
!
! Output:
! (rx, thetax, phix): location where field line cross Earth's surface 
! Vmx: magnetic poential
! ell: updated estimate of src length

use rksuite_90,   only : wp, rk_comm_real_1d, setup, range_integrate, &
                         collect_garbage
use mod_function, only : sgn, f

implicit none

integer, parameter :: neqn = 3

real(wp), intent(in) :: r0, theta0, phi0, s, ds

real(wp), intent(out) :: rx, thetax, phix, Vmx, ell

real(wp) :: r, theta, phi

real(wp) :: Brx, Btx, Bpx, Bx, dBrdt

real(wp) :: dx

type(rk_comm_real_1d) :: comm
real(kind=wp), dimension(neqn) :: y_start, y_got, thresholds
real(kind=wp) :: s_start, s_end, tolerance, s_want, s_inc, s_got
integer :: flag

s_start = 0.0_wp
s_end   = s
s_inc   = ds

y_start(1) = r0
y_start(2) = theta0
y_start(3) = phi0

tolerance  = 1.0e-8_wp
tolerance  = tol

! 'thresholds' can be set for each components
thresholds = 1.0e-8_wp 
thresholds = tol

call setup(comm,s_start,y_start,s_end,tolerance,thresholds)

r     = r0
theta = theta0
phi   = phi0

s_want = s_start

do

   rx     = r
   thetax = theta
   phix   = phi

   !print *, ' rx, thetax, phix = ', rx, thetax/dtr, phix/dtr

   s_want = s_want + s_inc

   if (s_want > s_end) exit

   call range_integrate(comm,f,s_want,s_got,y_got=y_got,flag=flag)

   !print *, ' s = ', s_got, ' y = ', y_got

   if (flag /= 1) exit

   r     = y_got(1)
   theta = y_got(2)
   phi   = y_got(3)

   ! check if it cross the Earth's surface

   if (r < a) then

      call igrf (rx, thetax, phix, Brx, Btx, Bpx, Vmx, dBrdt)

      Bx = sqrt(Brx*Brx + Btx*Btx + Bpx*Bpx)
      dx = -sgn*(rx - a)/(Brx/Bx)

      phix   = phix   + sgn * dx * (Bpx/Bx) / (rx*sin(thetax))
      thetax = thetax + sgn * dx * (Btx/Bx) / rx
      rx     = rx     + sgn * dx * (Brx/Bx)

      ! update (Br, Bt, Bp) and Vm
      call igrf (rx, thetax, phix, Brx, Btx, Bpx, Vmx, dBrdt)

      ! output 'ell' too: ok larger than the actual arc length
      ell = s_got

      exit

   endif

enddo

! release space from rksuite integrator
call collect_garbage(comm)

end subroutine tracing_down



subroutine fieldline (r0, theta0, phi0, s, ds, Vm, r, theta, phi)

! set up grid along field line: starting from apex

! Input:
! (r0, theta0, phi0): starting point (e.g., apex)
! s:  arc legnth
! ds: step size using in field line tracing
! Vm(:): scalar potential distributed along the field line
!
! Output:
! (r, theta, phi): location at Vm(:) along field line

use rksuite_90,   only : wp, rk_comm_real_1d, setup, range_integrate, &
                         collect_garbage
use mod_function, only : sgn, f

implicit none

integer, parameter :: neqn=3, kmax=100

real(wp), parameter :: hmin=80.0d0

real(wp), intent(in) :: r0, theta0, phi0, s, ds

real(wp), dimension(:), intent(in) :: Vm

real(wp), dimension(:), intent(out) :: r, theta, phi

type(rk_comm_real_1d) :: comm
real(kind=wp), dimension(neqn) :: y_start, y_got, thresholds
real(kind=wp) :: s_start, s_end, tolerance, s_want, s_inc, s_got
integer :: flag

real(wp), dimension(:), allocatable :: rx, thetax, phix, Vmx

real(wp) :: r1, theta1, phi1, r2, theta2, phi2

real(wp) :: Brx, Btx, Bpx, Bx, Vx, dBrdt

real(wp) :: aa, bb, ss, eps, fval

integer :: npts, nstep, nx, nm, np

integer :: i, n, icount, istat

logical :: bracket


npts = size(Vm)

!print *, ' npts = ', npts

! form a first guess of grids

nstep = ceiling(s/ds)

!print *, ' nstep = ', nstep

allocate (rx(nstep), thetax(nstep), phix(nstep), Vmx(nstep))

rx(1)     = r0
thetax(1) = theta0
phix(1)   = phi0

Vmx(1)    = Vm(1)

s_start = 0.0_wp
s_end   = s
s_inc   = ds

y_start(1) = r0
y_start(2) = theta0
y_start(3) = phi0

tolerance  = 1.0e-8_wp
tolerance  = tol

! 'thresholds' can be set for each components
thresholds = 1.0e-8_wp 
thresholds = tol

call setup(comm,s_start,y_start,s_end,tolerance,thresholds)

s_want = s_start

do n = 2, nstep

   s_want = s_want + s_inc

   if (s_want > s_end) exit

   call range_integrate(comm,f,s_want,s_got,y_got=y_got,flag=flag)

   !print *, ' s = ', s_got, ' y = ', y_got

   if (flag /= 1) exit

   rx(n)     = y_got(1)
   thetax(n) = y_got(2)
   phix(n)   = y_got(3)

   call igrf (rx(n), thetax(n), phix(n), Brx, Btx, Bpx, Vmx(n), dBrdt)

   !print *, ' n, rx, thetax, phix, Vmx = ', &
   !           n, rx(n), thetax(n)/dtr, phix(n)/dtr, Vmx(n)

   ! check if crossing Earth's surface

   !if (rx(n) < a) then
   ! need to cover a little more range
    if (rx(n) < a - hmin) then
   !if (rx(n) < a - 1   ) then
      nx = n
      exit
   else
      nx = nstep
   endif

enddo

! release space from rksuite integrator
call collect_garbage(comm)


! locate (r, theta, phi) at Vm from (rx, thetax, phix) and Vmx

! initialze to zero first

r = 0.0d0; theta = 0.0d0; phi = 0.0d0

r(1)     = r0
theta(1) = theta0
phi(1)   = phi0


eps = epsilon(eps)


aa = 0.0d0
!bb = ds
! incrase 'bracket' slightly since we use low-order integrator next 
! (high-order integrator is too expensive)
bb = 1.1_wp * ds

!print *, ' Vmx = ', Vmx

do i = 2, npts

!print *, ' i, Vm = ', i, Vm(i)


! bracket Vm(i)

bracket = .False.

do n = 2, nx

if ((sgn < 0.d0 .and. Vmx(n-1) <= Vm(i) .and. Vm(i) <= Vmx(n)) .or. &
    (sgn > 0.d0 .and. Vmx(n-1) >= Vm(i) .and. Vm(i) >= Vmx(n))) then

   nm = n - 1
   np = n

   r1     = rx(nm)
   theta1 = thetax(nm)
   phi1   = phix(nm)

   !print *, ' Vm-, Vm, Vm+: ', Vmx(n-1), Vm(i), Vmx(n)
   !print *, ' r1, theta1, phi1 = ', r1, theta1/dtr, phi1/dtr

   bracket = .True.

   exit

endif

enddo


!if (.not. bracket) cycle
if (.not. bracket) exit


!if (abs(Vm(i) - Vmx(nm)) <= epsilon(Vm(i))) then
!   r(i)     = rx(nm)
!   theta(i) = thetax(nm)
!   phi(i)   = phix(nm)
!   cycle
!else if (abs(Vm(i) - Vmx(np)) <= epsilon(Vm(i))) then
!   r(i)     = rx(np)
!   theta(i) = thetax(np)
!   phi(i)   = phix(np)
!   cycle
!endif



r2     = r1
theta2 = theta1
phi2   = phi1

!print *, ' r2, theta2, phi2 = ', r2, theta2/dtr, phi2/dtr

call igrf (r2, theta2, phi2, Brx, Btx, Bpx, Vx, dBrdt)

istat  = 0
icount = 0

do

   icount = icount + 1

   call zero_rc ( aa, bb, eps, ss, fval, istat )

   if ( istat < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ZERO_RC returned an error flag!'

     !print *, ' Vx, Vm(i), Vmx(nm), Vmx(np) = ', Vx, Vm(i), Vmx(nm), Vmx(np)

      exit
   endif

   Bx = sqrt(Brx*Brx + Btx*Btx + Bpx*Bpx)

   phi2   = phi1   + sgn * (ss-aa) * (Bpx/Bx) / (r2*sin(theta2))
   theta2 = theta1 + sgn * (ss-aa) * (Btx/Bx) / r2
   r2     = r1     + sgn * (ss-aa) * (Brx/Bx)

   !print *, ' r2, theta2, phi2 = ', r2, theta2, phi2

   call igrf (r2, theta2, phi2, Brx, Btx, Bpx, Vx, dBrdt)

   fval = Vx - Vm(i)

   !print *, ' istat, aa, bb, ss, fval = ', istat, aa, bb, ss, fval

   r(i)     = r2
   theta(i) = theta2
   phi(i)   = phi2

   if ( istat == 0 ) then
      !print *, ' istat, aa, bb, ss, fval = ', istat, aa, bb, ss, fval
      exit
   endif

   if (icount > kmax) then
      print *, '*** Warning: exceeding maximum iterations: ', kmax
      exit
   endif

enddo 

! if r is already close to earth's surface, exit
if (r(i) < a + hmin) exit

enddo 

deallocate (rx, thetax, phix, Vmx)

end subroutine fieldline



subroutine metric (r, theta, phi, s, ds, B,                      &
                   ct_mu, ct_chi, ct_phi, co_mu, co_chi, co_phi, &
                   ct_g, co_g, h1, h2, h3)

! compute magnetic field, basis vectors, and metric terms

! Input:
! (r, theta, phi): location of grid points along field line
! s:  (total) arc legnth
! ds: step size using in field line tracing
!
! Output:
! B: magnitude of magentic field
! ct_*: contravariant-basis vectors
! co_*: covariant-basis vectors
! ct_g, co_g: metric coefficients
! h1, h2, h3: scale factors

implicit none

real(wp), intent(in) :: s, ds

real(wp), dimension(:), intent(in) :: r, theta, phi

real(wp), dimension(:  ), intent(out) :: B, h1, h2, h3
real(wp), dimension(:,:), intent(out) :: ct_mu, ct_chi, ct_phi
real(wp), dimension(:,:), intent(out) :: co_mu, co_chi, co_phi
real(wp), dimension(:,:,:), intent(out) :: ct_g, co_g

real(wp), dimension(:), allocatable :: Br, Bt, Bp, Jac

real(wp) :: Vm, dBrdt

real(wp) :: dr, dtheta, dphi
real(wp) :: r1, theta1, phi1
real(wp) :: r2, theta2, phi2

real(wp) :: chip, chim, phip, phim

integer :: npts, nstep

integer :: i, j, n, istat, icount

npts = size(r)

!print *, ' npts = ', npts

allocate (Br(npts), Bt(npts), Bp(npts), Jac(npts))

! compute magnetic field

do n = 1, npts
call igrf (r(n), theta(n), phi(n), Br(n), Bt(n), Bp(n), Vm, dBrdt)
B(n) = sqrt(Br(n)*Br(n) + Bt(n)*Bt(n) + Bp(n)*Bp(n))
enddo

! compute basis vectors

do n = 1, npts

!dr = 1.0d0
 dr = 2.0d1
 dr = 1.0d1
call find_apex (r(n)+dr, theta(n), phi(n), s, ds, r1, theta1, phi1)
call find_apex (r(n)-dr, theta(n), phi(n), s, ds, r2, theta2, phi2)

!print *, ' n, r1, r2 = ', n, r1, r2

chip = 1.0d0/r1; chim = 1.0d0/r2
ct_chi(1,n) = (chip - chim) / (2.0d0*dr)

!print *, ' n, ct_chi(1) = ', n, ct_chi(1,n)

phip = phi1; phim = phi2
ct_phi(1,n) = (phip - phim) / (2.0d0*dr)

!print *, ' n, ct_phi(1) = ', n, ct_phi(1,n)

dtheta = 0.25d0*pi/180.0d0
dtheta = 0.1250*pi/180.0d0
call find_apex (r(n), theta(n)+dtheta, phi(n), s, ds, r1, theta1, phi1)
call find_apex (r(n), theta(n)-dtheta, phi(n), s, ds, r2, theta2, phi2)

!print *, ' n, r1, r2 = ', n, r1, r2

chip = 1.0d0/r1; chim = 1.0d0/r2
ct_chi(2,n) = (chip - chim) / (2.0d0*dtheta*r(n))

!print *, ' n, ct_chi(2) = ', n, ct_chi(2,n)

phip = phi1; phim = phi2
ct_phi(2,n) = (phip - phim) / (2.0d0*dtheta*r(n))

!print *, ' n, ct_phi(2) = ', n, ct_phi(2,n)

dphi = 0.25d0*pi/180.0d0
dphi = 0.1250*pi/180.0d0
call find_apex (r(n), theta(n), phi(n)+dphi, s, ds, r1, theta1, phi1)
call find_apex (r(n), theta(n), phi(n)-dphi, s, ds, r2, theta2, phi2)

!print *, ' n, r1, r2 = ', n, r1, r2

chip = 1.0d0/r1; chim = 1.0d0/r2
ct_chi(3,n) = (chip - chim) / (2.0d0*dphi*r(n)*sin(theta(n)))

!print *, ' n, ct_chi(3) = ', n, ct_chi(3,n)

phip = phi1; phim = phi2
ct_phi(3,n) = (phip - phim) / (2.0d0*dphi*r(n)*sin(theta(n)))

!print *, ' n, ct_phi(3) = ', n, ct_phi(3,n)

ct_mu(1,n) = -Br(n)
ct_mu(2,n) = -Bt(n)
ct_mu(3,n) = -Bp(n)

! Jacobian

Jac(n) = dot_product(ct_mu(:,n), cross_prod(ct_chi(:,n),ct_phi(:,n)))
Jac(n) = 1.0d0 / Jac(n)

!print *, ' n, Jac   = ', n, Jac(n)/(29619.4d0*6371.2d0**3)
!print *, ' n, 1/B^2 = ', n, 1.0/(B(n)*B(n))

co_mu (:,n) = cross_prod(ct_chi(:,n), ct_phi(:,n))*Jac(n)
co_chi(:,n) = cross_prod(ct_phi(:,n), ct_mu (:,n))*Jac(n)
co_phi(:,n) = cross_prod(ct_mu (:,n), ct_chi(:,n))*Jac(n)

!print *, ' n, co_mu  = ', n, co_mu (:,n)
!print *, ' n, co_chi = ', n, co_chi(:,n)
!print *, ' n, co_phi = ', n, co_phi(:,n)

! metric terms

ct_g(1,1,n) = dot_product(ct_mu (:,n), ct_mu (:,n))
ct_g(1,2,n) = dot_product(ct_mu (:,n), ct_chi(:,n))
ct_g(1,3,n) = dot_product(ct_mu (:,n), ct_phi(:,n))
ct_g(2,1,n) = ct_g(1,2,n)
ct_g(2,2,n) = dot_product(ct_chi(:,n), ct_chi(:,n))
ct_g(2,3,n) = dot_product(ct_chi(:,n), ct_phi(:,n))
ct_g(3,1,n) = ct_g(1,3,n)
ct_g(3,2,n) = ct_g(2,3,n)
ct_g(3,3,n) = dot_product(ct_phi(:,n), ct_phi(:,n))

co_g(1,1,n) = dot_product(co_mu (:,n), co_mu (:,n))
co_g(1,2,n) = dot_product(co_mu (:,n), co_chi(:,n))
co_g(1,3,n) = dot_product(co_mu (:,n), co_phi(:,n))
co_g(2,1,n) = co_g(1,2,n)
co_g(2,2,n) = dot_product(co_chi(:,n), co_chi(:,n))
co_g(2,3,n) = dot_product(co_chi(:,n), co_phi(:,n))
co_g(3,1,n) = co_g(1,3,n)
co_g(3,2,n) = co_g(2,3,n)
co_g(3,3,n) = dot_product(co_phi(:,n), co_phi(:,n))

h1(n) = sqrt(co_g(1,1,n))
h2(n) = sqrt(co_g(2,2,n))
h3(n) = sqrt(co_g(3,3,n))

! reciprocity

!print *, ' n, reciprocity: mu  ', n, dot_product(co_mu (:,n),ct_mu (:,n)), &
!                                     dot_product(co_mu (:,n),ct_chi(:,n)), &
!                                     dot_product(co_mu (:,n),ct_phi(:,n))

!print *, ' n, reciprocity: chi ', n, dot_product(co_chi(:,n),ct_chi(:,n)), &
!                                     dot_product(co_chi(:,n),ct_phi(:,n)), &
!                                     dot_product(co_chi(:,n),ct_mu (:,n))

!print *, ' n, reciprocity: phi ', n, dot_product(co_phi(:,n),ct_phi(:,n)), &
!                                     dot_product(co_phi(:,n),ct_mu (:,n)), &
!                                     dot_product(co_phi(:,n),ct_chi(:,n))

! are they orthogonal?

!print *, ' n, othorgonality: ', n, dot_product(ct_chi(:,n), ct_phi(:,n)), &
!                                   dot_product(ct_phi(:,n), ct_mu (:,n)), &
!                                   dot_product(ct_mu (:,n), ct_chi(:,n))

! are they Euler potentials?

!print *, ' n, ct_mu =              ', n, ct_mu(:,n)

!print *, ' n, ct_chi x ct_phi =    ', n, 29619.4d0*(6371.2d0**3) * &
!           cross_prod(ct_chi(:,n), ct_phi(:,n))

!print *, ' n, relative error (%) = ', n, 100 * abs((29619.4d0*6371.2d0**3 * &
!           cross_prod(ct_chi(:,n), ct_phi(:,n)) - ct_mu(:,n)) / ct_mu(:,n))

enddo

deallocate (Br, Bt, Bp, Jac)

end subroutine metric



subroutine find_apex (r0, theta0, phi0, s, ds, r_a, theta_a, phi_a)

! find apex using rkf45 integrator

! Input:
! (r0, theta0, phi0): starting point
! s:  arc legnth
! ds: step size using in field line tracing
!
! Output:
! (r_a, theta_a, phi_a): coordinates at apex

use rksuite_90,   only : wp, rk_comm_real_1d, setup, range_integrate, &
                         collect_garbage
use mod_function, only : sgn, f

implicit none

integer, parameter :: neqn=3, kmax=100

real(wp), intent(in) :: r0, theta0, phi0, s, ds

real(wp), intent(out) :: r_a, theta_a, phi_a

real(wp) :: Br, Bt, Bp, B, Vm, dBrdt
real(wp) :: r, theta, phi

real(wp) :: aa, bb, ss, eps, fval

real(wp) :: rx, thetax, phix, Brx, Btx, Bpx, Bx

real(wp) :: r1, theta1, phi1

integer :: n, istat, icount

type(rk_comm_real_1d) :: comm
real(kind=wp), dimension(neqn) :: y_start, y_got, thresholds
real(kind=wp) :: s_start, s_end, tolerance, s_want, s_inc, s_got
integer :: flag

s_start = 0.0_wp
s_inc   = ds
s_end   = s

y_start(1) = r0
y_start(2) = theta0
y_start(3) = phi0

tolerance  = 1.0e-8_wp
tolerance  = tol

! 'thresholds' can be set for each components
thresholds = 1.0e-8_wp
thresholds = tol

r     = r0
theta = theta0
phi   = phi0

call igrf(r0, theta0, phi0, Br, Bt, Bp, Vm, dBrdt)

sgn = sign(1.0d0, Br)

call setup(comm,s_start,y_start,s_end,tolerance,thresholds)

s_want = s_start

icount = 0

do

   icount = icount + 1

   rx     = r
   thetax = theta
   phix   = phi

   r1     = r
   theta1 = theta
   phi1   = phi

   Brx = Br
   Btx = Bt
   Bpx = Bp

   s_want = s_want + s_inc

   if (s_want > s_end) exit

   call range_integrate(comm,f,s_want,s_got,y_got=y_got,flag=flag)

   if (flag /= 1) exit

   !print *, ' s = ', s_got, ' y = ', y_got

   r     = y_got(1)
   theta = y_got(2)
   phi   = y_got(3)

   call igrf(y_got(1), y_got(2), y_got(3), Br, Bt, Bp, Vm, dBrdt)

   !print *, ' Br, Bt, Bp, Vm, dBrdt = ', Br, Bt, Bp, Vm, dBrdt

   ! bracket 'apex'

   if (Br*Brx <= 0.0) then
      aa = 0.0_wp

      !bb = ds
      ! incrase 'bracket' slightly since we use low-order integrator next 
      ! (high-order integrator is too expensive)
      bb = 1.1_wp * ds

      exit
   endif

enddo

! release space from rksuite integrator
call collect_garbage(comm)

! find coordinates at apex

call igrf(rx, thetax, phix, Br, Bt, Bp, Vm, dBrdt)

eps = epsilon(eps)

istat  = 0
icount = 0

do

   icount = icount + 1

   call zero_rc ( aa, bb, eps, ss, fval, istat )

   if ( istat < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ZERO_RC returned an error flag!'
      exit
   endif

   Bx = sqrt(Brx*Brx + Btx*Btx + Bpx*Bpx)

   phix   = phi1   + sgn * (ss-aa) * (Bpx/Bx) / (rx*sin(thetax))
   thetax = theta1 + sgn * (ss-aa) * (Btx/Bx) / rx
   rx     = r1     + sgn * (ss-aa) * (Brx/Bx)

   call igrf(rx, thetax, phix, Brx, Btx, Bpx, Vm, dBrdt)

   fval = Brx

   !print *, ' istat, aa, bb, ss, fval = ', istat, aa, bb, ss, fval

   r_a     = rx
   theta_a = thetax
   phi_a   = phix

   !print *, ' rx, thetax, phix = ', rx, thetax, phix

   if ( istat == 0 ) then
      !print *, ' istat, aa, bb, ss, fval = ', istat, aa, bb, ss, fval
      exit
   endif

   if (icount > kmax) then
      print *, '*** Warning: exceeding maximum iterations: ', kmax
      exit
   endif

enddo

end subroutine find_apex



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



function cross_prod(a,b) result(axb)

implicit none

real(wp),dimension(3),intent(in) :: a
real(wp),dimension(3),intent(in) :: b

real(wp),dimension(3) :: axb

axb(1) = a(2)*b(3) - a(3)*b(2)
axb(2) = a(3)*b(1) - a(1)*b(3)
axb(3) = a(1)*b(2) - a(2)*b(1)

end function cross_prod



end module mod_fieldline



