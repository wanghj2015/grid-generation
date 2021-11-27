module mod_igrf

use mod_precision, only : wp=>dp

implicit none

integer, parameter :: Nmax=13, Mmax=13
integer, parameter :: Nepoch=25

real(wp), dimension(Nmax,0:Mmax,Nepoch) :: g, h
real(wp), dimension(Nmax,0:Mmax       ) :: sv

real(wp), dimension(Nmax,0:Mmax) :: gx, hx

real(wp) :: epochs(Nepoch)
data epochs/1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940, 1945, &
            1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995, &
            2000, 2005, 2010, 2015, 2020/

! Earth’s mean reference spherical radius [km] 
real(wp), parameter :: a = 6371.2_wp


contains


subroutine get_igrf_coef (epoch)

! read in the IGRF coefficients
! and interpolate to the right epoch

implicit none

real(wp), intent(in) :: epoch

real(wp) :: wgt(2)

character(1) :: c

integer :: i, j, k, n, m, im, ip

! initialize to zero

g(:,:,:) = 0.0_wp
h(:,:,:) = 0.0_wp
sv(:,: ) = 0.0_wp

! read in the IGRF coefficients

open (10, file='igrf13coeffs.txt', status='old')

! skip the first 3 lines
do i = 1, 4
read (10,*)
enddo

!print *, ' epochs = ', epochs

!stop

do n = 1, Nmax
read (10,*) c, j, k, g(n,0,:), sv(n,0)

!print *, c, j, k, g(n,0,:), sv(n,0)

do m = 1, n
read (10,*) c, j, k, g(n,m,:), sv(n,m)
read (10,*) c, j, k, h(n,m,:), sv(n,m)

!print *, c, j, k, g(n,m,:), sv(n,m)
!print *, c, j, k, h(n,m,:), sv(n,m)

enddo
enddo

close (10)

! interpolate to the right epoch

! need to interpolate to the right epoch

! bracket the epoch

if (epoch < epochs(1) .or. epoch > epochs(Nepoch) + 5.0_wp) then

   print *, ' epoch = ', epoch, ' is out of bounds, exiting now! '
   stop

elseif (epoch >= epochs(1) .and. epoch <= epochs(Nepoch)) then

   do i = 2, Nepoch
   if (epochs(i-1) <= epoch .and. epoch <= epochs(i)) then
      im = i-1
      ip = i
      exit
   endif 
   enddo

   wgt(2) = (epoch - epochs(im)) / (epochs(ip) - epochs(im))
   wgt(1) = 1.0_wp - wgt(2)

   gx(:,:) = wgt(1) * g(:,:,im) + wgt(2) * g(:,:,ip)
   hx(:,:) = wgt(1) * h(:,:,im) + wgt(2) * h(:,:,ip)

else

   do n = 1, Nmax
   do m = 0, n
   gx(n,m) = g(n,m,Nepoch) + sv(n,m) * (epoch - epochs(Nepoch))
   enddo
   enddo

endif


! for a dipole only
 do n = 1, Nmax
 do m = 0, n
 if (.not. (n == 1 .and. m == 0)) then
    gx(n,m) = 0.0_wp
 endif
 enddo
 enddo
 hx(:,:) = 0.0_wp


!do n = 1, Nmax
!do m = 0, n
!print *, 'n, m, gx = ', n, m, gx(n,m)
!print *, 'n, m, hx = ', n, m, hx(n,m)
!enddo
!enddo

end subroutine get_igrf_coef
  

end module mod_igrf


