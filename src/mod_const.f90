module mod_const

use mod_precision, only : wp=>dp

implicit none

! physical constants

real(wp), parameter :: grav = 9.81_wp ! m s^-2
real(wp), parameter :: k_b = 1.381d-23 ! J K^-1
real(wp), parameter :: r_earth = 6.3712d6 ! m

! atomic/molecular weights

real(wp), parameter :: o_amu = 15.999_wp
real(wp), parameter :: h_amu = 1.008_wp
real(wp), parameter :: he_amu = 4.003_wp
real(wp), parameter :: n_amu = 14.007_wp
real(wp), parameter :: o2_amu = o_amu * 2.0_wp
real(wp), parameter :: n2_amu = n_amu * 2.0_wp
real(wp), parameter :: no_amu = n_amu + o_amu

! other constants

real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp), dtr = pi/180.0_wp

end module mod_const

