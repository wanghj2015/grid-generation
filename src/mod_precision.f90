module mod_precision

! see http://fortranwiki.org/fortran/show/Real+precision

use iso_fortran_env, only: real32, real64, real128

implicit none

integer, parameter :: sp = real32
integer, parameter :: dp = real64
integer, parameter :: qp = real128

!integer, parameter :: sp = selected_real_kind(6, 37)
!integer, parameter :: dp = selected_real_kind(15, 307)
!integer, parameter :: qp = selected_real_kind(33, 4931)

end module mod_precision
