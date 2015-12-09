SUBROUTINE cdfnor ( which, p, q, x, mean, sd, status, bound )

!*****************************************************************************80
!
!! CDFNOR evaluates the CDF of the Normal distribution.
!
!  Discussion:
!
!    A slightly modified version of ANORM from SPECFUN
!    is used to calculate the cumulative standard normal distribution.
!
!    The rational functions from pages 90-95 of Kennedy and Gentle
!    are used as starting values to a Newton iteration which
!    compute the inverse standard normal.  Therefore no searches are
!    necessary for any parameter.
!
!    For X < -15, the asymptotic expansion for the normal is used  as
!    the starting value in finding the inverse standard normal.
!
!    The normal density is proportional to
!    exp ( - 0.5D+00 * (( X - MEAN)/SD)**2)
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions
!    1966, Formula 26.2.12.
!
!    William Cody,
!    Algorithm 715:
!    SPECFUN - A Portable FORTRAN Package of
!    Special Function Routines and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, pages 22-32, 1993.
!
!    William Kennedy, James Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) WHICH, indicates which argument is to be
!    calculated from the others.
!    1: Calculate P and Q from X, MEAN and SD;
!    2: Calculate X from P, Q, MEAN and SD;
!    3: Calculate MEAN from P, Q, X and SD;
!    4: Calculate SD from P, Q, X and MEAN.
!
!    Input/output, real ( kind = 8 ) P, the integral from -infinity to X
!    of the Normal density.  If this is an input or output value, it will
!    lie in the range [0,1].
!
!    Input/output, real ( kind = 8 ) Q, equal to 1-P.  If Q is an input
!    value, it should lie in the range [0,1].  If Q is an output value,
!    it will lie in the range [0,1].
!
!    Input/output, real ( kind = 8 ) X, the upper limit of integration of
!    the Normal density.
!
!    Input/output, real ( kind = 8 ) MEAN, the mean of the Normal density.
!
!    Input/output, real ( kind = 8 ) SD, the standard deviation of the
!    Normal density.  If this is an input value, it should lie in the
!    range (0,+infinity).
!
!    Output, integer ( kind = 4 ) STATUS, the status of the calculation.
!    0, if calculation completed correctly;
!    -I, if input parameter number I is out of range;
!    1, if answer appears to be lower than lowest search bound;
!    2, if answer appears to be higher than greatest search bound;
!    3, if P + Q /= 1.
!
!    Output, real ( kind = 8 ) BOUND, is only defined if STATUS is nonzero.
!    If STATUS is negative, then this is the value exceeded by parameter I.
!    if STATUS is 1 or 2, this is the search bound that was exceeded.
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) dinvnr
  real ( kind = 8 ) mean
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) sd
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) x
  real ( kind = 8 ) z

  status = 0
  bound = 0.0D+00
!
!  Check the arguments.
!
  status = 0

  if ( which < 1 ) then
    status = -1
    bound = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDFNOR - Fatal error!'
    write ( *, '(a)' ) '  The input parameter WHICH is out of range.'
    write ( *, '(a)' ) '  Legal values are between 1 and 4.'
    return
  else if ( 4 < which ) then
    status = -1
    bound = 4.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDFNOR - Fatal error!'
    write ( *, '(a)' ) '  The input parameter WHICH is out of range.'
    write ( *, '(a)' ) '  Legal values are between 1 and 4.'
    return
  end if
!
!  Unless P is to be computed, make sure it is legal.
!
  if ( which /= 1 ) then
    if ( p < 0.0D+00 ) then
      status = -2
      bound = 0.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CDFNOR - Fatal error!'
      write ( *, '(a)' ) '  Input parameter P is out of range.'
      return
    else if ( 1.0D+00 < p ) then
      status = -2
      bound = 1.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CDFNOR - Fatal error!'
      write ( *, '(a)' ) '  Input parameter P is out of range.'
      return
    end if
  end if
!
!  Unless Q is to be computed, make sure it is legal.
!
  if ( which /= 1 ) then
    if ( q < 0.0D+00 ) then
      status = -3
      bound = 0.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CDFNOR - Fatal error!'
      write ( *, '(a)' ) '  Input parameter Q is out of range.'
      return
    else if ( 1.0D+00 < q ) then
      status = -3
      bound = 1.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CDFNOR - Fatal error!'
      write ( *, '(a)' ) '  Input parameter Q is out of range.'
      return
    end if
  end if
!
!  Check that P + Q = 1.
!
  if ( which /= 1 ) then
    if ( 3.0D+00 * epsilon ( p ) < abs ( ( p + q ) - 1.0D+00 ) ) then
      status = 3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CDFNOR - Fatal error!'
      write ( *, '(a)' ) '  P + Q /= 1.'
      return
    end if
  end if

  if ( which /= 4 ) then
    if ( sd <= 0.0D+00 ) then
      bound = 0.0D+00
      status = -6
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CDFNOR - Fatal error!'
      write ( *, '(a)' ) '  Input parameter SD is out of range.'
      return
    end if
  end if
!
!  Calculate P and Q.
!
  if ( which == 1 ) then

    z = ( x - mean ) / sd
    call cumnor ( z, p, q )
!
!  Calculate X.
!
  else if ( which == 2 ) then

    z = dinvnr ( p, q )
    x = sd * z + mean
!
!  Calculate MEAN.
!
  else if ( which == 3 ) then

    z = dinvnr ( p, q )
    mean = x - sd * z
!
!  Calculate SD.
!
  else if ( which == 4 ) then

    z = dinvnr ( p, q )
    sd = ( x - mean ) / z

  end if

  return
end

SUBROUTINE cumnor ( arg, cum, ccum )

!*****************************************************************************80
!
!! CUMNOR computes the cumulative normal distribution.
!
!  Discussion:
!
!    This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!    This transportable program uses rational functions that
!    theoretically approximate the normal distribution function to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine dependent
!    constants.
!
!  Reference:
!
!    William Cody,
!    Rational Chebyshev approximations for the error function,
!    Mathematics of Computation,
!    1969, pages 631-637.
!
!    William Cody,
!    Algorithm 715:
!    SPECFUN - A Portable FORTRAN Package of Special Function Routines
!    and Test Drivers,
!    ACM Transactions on Mathematical Software,
!    Volume 19, Number 1, 1993, pages 22-32.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the upper limit of integration.
!
!    Output, real ( kind = 8 ) CUM, CCUM, the Normal density CDF and
!    complementary CDF.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) EPS, the argument below which anorm(x)
!    may be represented by 0.5 and above which  x*x  will not underflow.
!    A conservative value is the largest machine number X
!    such that   1.0D+00 + X = 1.0D+00   to machine precision.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287D+00, &
    1.6102823106855587881D+02, &
    1.0676894854603709582D+03, &
    1.8154981253343561249D+04, &
    6.5682337918207449113D-02 /)
  real ( kind = 8 ) arg
  real ( kind = 8 ), parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870D+01, &
    9.7609855173777669322D+02, &
    1.0260932208618978205D+04, &
    4.5507789335026729956D+04 /)
  real ( kind = 8 ), parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764D-01, &
    8.8831497943883759412D+00, &
    9.3506656132177855979D+01, &
    5.9727027639480026226D+02, &
    2.4945375852903726711D+03, &
    6.8481904505362823326D+03, &
    1.1602651437647350124D+04, &
    9.8427148383839780218D+03, &
    1.0765576773720192317D-08 /)
  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691D+01, &
    2.3538790178262499861D+02, &
    1.5193775994075548050D+03, &
    6.4855582982667607550D+03, &
    1.8615571640885098091D+04, &
    3.4900952721145977266D+04, &
    3.8912003286093271411D+04, &
    1.9685429676859990727D+04 /)
  real ( kind = 8 ) del
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699D-01, &
    1.274011611602473639D-01, &
    2.2235277870649807D-02, &
    1.421619193227893466D-03, &
    2.9112874951168792D-05, &
    2.307344176494017303D-02 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121D+00, &
    4.68238212480865118D-01, &
    6.59881378689285515D-02, &
    3.78239633202758244D-03, &
    7.29751555083966205D-05 /)
  real ( kind = 8 ), parameter :: root32 = 5.656854248D+00
  real ( kind = 8 ), parameter :: sixten = 16.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ), parameter :: sqrpi = 3.9894228040143267794D-01
  real ( kind = 8 ), parameter :: thrsh = 0.66291D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) xden
  real ( kind = 8 ) xnum
  real ( kind = 8 ) y
  real ( kind = 8 ) xsq
!
!  Machine dependent constants
!
  eps = epsilon ( 1.0D+00 ) * 0.5D+00

  x = arg
  y = abs ( x )

  if ( y <= thrsh ) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
    if ( eps < y ) then
      xsq = x * x
    else
      xsq = 0.0D+00
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do
    cum = x * ( xnum + a(4) ) / ( xden + b(4) )
    temp = cum
    cum = 0.5D+00 + temp
    ccum = 0.5D+00 - temp
!
!  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
!
  else if ( y <= root32 ) then

    xnum = c(9) * y
    xden = y
    do i = 1, 7
      xnum = ( xnum + c(i) ) * y
      xden = ( xden + d(i) ) * y
    end do
    cum = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call r8_swap ( cum, ccum )
    end if
!
!  Evaluate ANORM for sqrt(32) < |X|.
!
  else

    cum = 0.0D+00
    xsq = 1.0D+00 / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    do i = 1, 4
      xnum = ( xnum + p(i) ) * xsq
      xden = ( xden + q(i) ) * xsq
    end do

    cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    cum = ( sqrpi - cum ) / y
    xsq = aint ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    cum = exp ( - xsq * xsq * 0.5D+00 ) &
      * exp ( - del * 0.5D+00 ) * cum
    ccum = 1.0D+00 - cum

    if ( 0.0D+00 < x ) then
      call r8_swap ( cum, ccum )
    end if

  end if

  if ( cum < tiny ( cum ) ) then
    cum = 0.0D+00
  end if

  if ( ccum < tiny ( ccum ) ) then
    ccum = 0.0D+00
  end if

  return
end

subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end

function dinvnr ( p, q )

!*****************************************************************************80
!
!! DINVNR computes the inverse of the normal distribution.
!
!  Discussion:
!
!    This routine returns X such that
!
!      CUMNOR(X) = P,
!
!    that is, so that
!
!      P = integral ( -oo <= T <= X ) exp(-U*U/2)/sqrt(2*PI) dU
!
!    The rational function on page 95 of Kennedy and Gentle is used as a
!    starting value for the Newton method of finding roots.
!
!  Reference:
!
!    William Kennedy, James Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, Q, the probability, and the complementary
!    probability.
!
!    Output, real ( kind = 8 ) DINVNR, the argument X for which the
!    Normal CDF has the value P.
!
  implicit none

  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ) dinvnr
  real ( kind = 8 ) dx
  real ( kind = 8 ), parameter :: eps = 1.0D-13
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: maxit = 100
  real ( kind = 8 ) p
  real ( kind = 8 ) pp
  real ( kind = 8 ) q
  real ( kind = 8 ), parameter :: r2pi = 0.3989422804014326D+00
  real ( kind = 8 ) strtx
  real ( kind = 8 ) stvaln
  real ( kind = 8 ) xcur

  pp = min ( p, q )
  strtx = stvaln ( pp )
  xcur = strtx
!
!  Newton iterations.
!
  do i = 1, maxit

    call cumnor ( xcur, cum, ccum )
    dx = ( cum - pp ) / ( r2pi * exp ( -0.5D+00 * xcur * xcur ) )
    xcur = xcur - dx

    if ( abs ( dx / xcur ) < eps ) then
      if ( p <= q ) then
        dinvnr = xcur
      else
        dinvnr = -xcur
      end if
      return
    end if

  end do

  if ( p <= q ) then
    dinvnr = strtx
  else
    dinvnr = -strtx
  end if

  return
end

function stvaln ( p )

!*****************************************************************************80
!
!! STVALN provides starting values for the inverse of the normal distribution.
!
!  Discussion:
!
!    The routine returns an X for which it is approximately true that
!      P = CUMNOR(X),
!    that is,
!      P = Integral ( -infinity < U <= X ) exp(-U*U/2)/sqrt(2*PI) dU.
!
!  Reference:
!
!    William Kennedy, James Gentle,
!    Statistical Computing,
!    Marcel Dekker, NY, 1980, page 95,
!    QA276.4 K46
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the probability whose normal deviate
!    is sought.
!
!    Output, real ( kind = 8 ) STVALN, the normal deviate whose probability
!    is approximately P.
!
  implicit none

  real ( kind = 8 ) eval_pol
  real ( kind = 8 ) p
  real ( kind = 8 ) sgn
  real ( kind = 8 ) stvaln
  real ( kind = 8 ), parameter, dimension(0:4) :: xden = (/ &
    0.993484626060D-01, &
    0.588581570495D+00, &
    0.531103462366D+00, &
    0.103537752850D+00, &
    0.38560700634D-02 /)
  real ( kind = 8 ), parameter, dimension(0:4) :: xnum = (/ &
    -0.322232431088D+00, &
    -1.000000000000D+00, &
    -0.342242088547D+00, &
    -0.204231210245D-01, &
    -0.453642210148D-04 /)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  if ( p <= 0.5D+00 ) then

    sgn = -1.0D+00
    z = p

  else

    sgn = 1.0D+00
    z = 1.0D+00 - p

  end if

  y = sqrt ( -2.0D+00 * log ( z ) )
  stvaln = y + eval_pol ( xnum, 4, y ) / eval_pol ( xden, 4, y )
  stvaln = sgn * stvaln

  return
end

function eval_pol ( a, n, x )

!*****************************************************************************80
!
!! EVAL_POL evaluates a polynomial at X.
!
!  Discussion:
!
!    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
!
!  Modified:
!
!    15 December 1999
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(0:N), coefficients of the polynomial.
!
!    Input, integer ( kind = 4 ) N, length of A.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) EVAL_POL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  real ( kind = 8 ) eval_pol
  integer ( kind = 4 ) i
  real ( kind = 8 ) term
  real ( kind = 8 ) x

  term = a(n)
  do i = n - 1, 0, -1
    term = term * x + a(i)
  end do

  eval_pol = term

  return
end


