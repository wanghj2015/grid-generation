module mod_function

implicit none

real(8), save :: sgn


contains


function f (s, y)

use rksuite_90_prec, only : wp
use mod_igrf,        only : igrf

implicit none

real(kind=wp), intent(in) :: s
real(kind=wp), dimension(:), intent(in) :: y ! = (r, theta, phi)
real(kind=wp), dimension(size(y))       :: f ! = dy/ds

 real(kind=wp) :: r, theta, phi
 real(kind=wp) :: Br, Bt, Bp, B, Vm, dBrdt
!real(8) :: r, theta, phi
!real(8) :: Br, Bt, Bp, B, Vm, dBrdt

r     = y(1)
theta = y(2)
phi   = y(3)

! compute Br, Bt, Bp & B from y = (r, theta, phi)
call igrf(r, theta, phi, Br, Bt, Bp, Vm, dBrdt)

B = sqrt(Br*Br + Bt*Bt + Bp*Bp)

!print *, ' Br, Bt, Bp = ', Br, Bt, Bp

f(1) = sgn*Br/B
f(2) = sgn*Bt/B/r
f(3) = sgn*Bp/B/r/sin(theta)

end function f


end module mod_function


