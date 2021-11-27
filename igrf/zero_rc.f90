subroutine zero_rc ( a, b, t, arg, fval, istat )

!*****************************************************************************80
!
!! ZERO_RC seeks the root of a function F(X) using reverse communication.
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one fval C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    The routine is a revised version of the Brent zero finder 
!    algorithm, using reverse communication.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the fval of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
!    sign interval.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Output, real ( kind = 8 ) ARG, the currently considered point.  The user
!    does not need to initialize this fval.  On return with istat positive,
!    the user is requested to evaluate the function at ARG, and return
!    the fval in fval.  On return with istat zero, ARG is the routine's
!    estimate for the function's zero.
!
!    Input/output, integer ( kind = 4 ) istat, used to communicate between 
!    the user and the routine.  The user only sets istat to zero on the first 
!    call, to indicate that this is a startup call.  The routine returns istat
!    positive to request that the function be evaluated at ARG, or returns
!    istat as 0, to indicate that the iteration is complete and that
!    ARG is the estimated zero
!
!    Input, real ( kind = 8 ) fval, the function fval at ARG, as requested
!    by the routine on the previous call.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ), save :: c
  real ( kind = 8 ), save :: d
  real ( kind = 8 ), save :: e
  real ( kind = 8 ), save :: fa
  real ( kind = 8 ), save :: fb
  real ( kind = 8 ), save :: fc
  real ( kind = 8 ) m
  real ( kind = 8 ), save :: machep
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ), save :: sa
  real ( kind = 8 ), save :: sb
  integer ( kind = 4 ) istat
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) fval
!
!  Input istat = 0.
!  Initialize, request F(A).
!
  if ( istat == 0 ) then

    machep = epsilon ( a )

    sa = a
    sb = b
    e = sb - sa
    d = e

    istat = 1
    arg = a
    return
!
!  Input istat = 1.
!  Receive F(A), request F(B).
!
  else if ( istat == 1 ) then

    fa = fval

    istat = 2
    arg = sb
    return
!
!  Input istat = 2
!  Receive F(B).
!
  else if ( istat == 2 ) then

    fb = fval

    if ( 0.0D+00 < fa * fb ) then
      istat = -1
      return
    end if

    c = sa
    fc = fa

  else

    fb = fval

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if

  end if
!
!  Compute the next point at which a function fval is requested.
!
  if ( abs ( fc ) < abs ( fb ) ) then

    sa = sb
    sb = c
    c = sa
    fa = fb
    fb = fc
    fc = fa

  end if

  tol = 2.0D+00 * machep * abs ( sb ) + t
  m = 0.5D+00 * ( c - sb )

  if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
    istat = 0
    arg = sb
    return
  end if

  if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

    e = m
    d = e

  else

    s = fb / fa

    if ( sa == c ) then

      p = 2.0D+00 * m * s
      q = 1.0D+00 - s

    else

      q = fa / fc
      r = fb / fc
      p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
      q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

    end if

    if ( 0.0D+00 < p ) then
      q = - q
    else
      p = - p
    end if

    s = e
    e = d

    if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
      p < abs ( 0.5D+00 * s * q ) ) then
      d = p / q
    else
      e = m
      d = e
    end if

  end if

  sa = sb
  fa = fb

  if ( tol < abs ( d ) ) then
    sb = sb + d
  else if ( 0.0D+00 < m ) then
    sb = sb + tol
  else
    sb = sb - tol
  end if

  arg = sb
  istat = istat + 1

end subroutine zero_rc


