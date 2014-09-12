
subroutine naca4_cambered ( m, p, t, c, n, xc, xu, yu, xl, yl )

!*****************************************************************************80
!
!! NACA4_CAMBERED: (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
!    "The characteristics of 78 related airfoil sections from tests in
!    the variable-density wind tunnel",
!    NACA Report 460, 1933.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) M, the maximum camber.
!    0.0 < M.
!
!    Input, real ( kind = 8 ) P, the location of maximum camber.
!    0.0 < P < 1.0
!
!    Input, real ( kind = 8 ) T, the maximum relative thickness.
!    0.0 < T <= 1.0
!
!    Input, real ( kind = 8 ) C, the chord length.
!    0.0 < C.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!
!    Input, real ( kind = 8 ) XC(N), points along the chord length.  
!    0.0 <= XC(*) <= C.
!
!    Output, real ( kind = 8 ) XU(N), YU(N), XL(N), YL(N), for each value of 
!    XC, measured along the camber line, the corresponding values (XU,YU) 
!    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) divisor
  real ( kind = 8 ) dycdx
  integer ( kind = 4 ) i
  real ( kind = 8 ) m
  real ( kind = 8 ) p
  real ( kind = 8 ) t
  real ( kind = 8 ) theta
  real ( kind = 8 ) xc(n)
  real ( kind = 8 ) xl(n)
  real ( kind = 8 ) xu(n)
  real ( kind = 8 ) yc
  real ( kind = 8 ) yl(n)
  real ( kind = 8 ) yt
  real ( kind = 8 ) yu(n)

  do i = 1, n

    if ( 0.0D+00 <= xc(i) / c .and. xc(i) / c <= p ) then
      divisor = p ** 2
    else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0D+00 ) then
      divisor = ( 1.0D+00 - p ) ** 2
    else
      divisor = 1.0D+00
    end if

    dycdx = 2.0D+00 * m * ( p - xc(i) / c ) / divisor

    theta = atan ( dycdx )
   
    yt = 5.0D+00 * t * c * ( &
       0.2969D+00 * sqrt ( xc(i) / c ) &
       + (((( &
         - 0.1015D+00 ) * ( xc(i) / c ) &
         + 0.2843D+00 ) * ( xc(i) / c ) &
         - 0.3516D+00 ) * ( xc(i) / c ) &
         - 0.1260D+00 ) * ( xc(i) / c ) )

    if ( 0.0D+00 <= xc(i) / c .and. xc(i) / c <= p ) then
      yc = m * xc(i) * ( 2.0D+00 * p - xc(i) / c ) / p ** 2
    else if ( p <= xc(i) / c .and. xc(i) / c <= 1.0D+00 ) then
      yc = m * ( xc(i) - c ) * ( 2.0D+00 * p - xc(i) / c - 1.0D+00 ) &
        / ( 1.0D+00 - p ) ** 2
    else
      yc = 0.0D+00
    end if

    xu(i) = xc(i) - yt * sin ( theta )
    yu(i) = yc + yt * cos ( theta )
    xl(i) = xc(i) + yt * sin ( theta )
    yl(i) = yc - yt * cos ( theta )

  end do

  return
end

subroutine naca4_symmetric ( t, c, n, x, y )

!*****************************************************************************80
!
!! NACA4_SYMMETRIC evaluates y(x) for a NACA symmetric 4-digit airfoil.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
!    "The characteristics of 78 related airfoil sections from tests in
!    the variable-density wind tunnel",
!    NACA Report 460, 1933.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the maximum relative thickness.
!
!    Input, real ( kind = 8 ) C, the chord length.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!
!    Input, real ( kind = 8 ) X(N), points along the chord length.  
!    0.0 <= X(*) <= C.
!
!    Output, real ( kind = 8 ) Y(N), for each value of X, the corresponding
!    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
!    lower wing surface.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 5.0D+00 * t * c * ( &
    0.2969D+00 * sqrt ( x(1:n) / c ) &
    + (((( &
      - 0.1015D+00 ) * ( x(1:n) / c ) &
      + 0.2843D+00 ) * ( x(1:n) / c ) &
      - 0.3516D+00 ) * ( x(1:n) / c ) &
      - 0.1260D+00 ) * ( x(1:n) / c ) )

  return
end

