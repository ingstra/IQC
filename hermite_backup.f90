
MODULE hermite
 integer, parameter :: dp = selected_real_kind(15, 307)

CONTAINS
subroutine hn_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
!
!  Discussion:
!
!    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.
!
!    These polynomials satisfy the orthonormality condition:
!
!      Integral ( -oo < X < +oo ) 
!        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) P(M,0:N), the values of the polynomials of 
!    index 0 through N.
!
  implicit none

  integer, intent(in) :: m, n
  real(dp), intent(in) :: x(m)
  real(dp), intent(out) :: p(m,0:n)

  integer :: j
  real(dp) :: fact, two_power
  real(dp), parameter :: r8_pi = 3.141592653589793D+00


  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = 2.0D+00 * x(1:m)
 
  do j = 2, n
    p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) &
      - 2.0D+00 * real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
!
!  Normalize.
!
  fact = 1.0D+00
  two_power = 1.0D+00
  do j = 0, n
    p(1:m,j) = p(1:m,j) / sqrt ( fact * two_power * sqrt ( r8_pi ) )
    fact = fact * real ( j + 1, kind = 8 )
    two_power = two_power * 2.0D+00
  end do

  return
end

END MODULE hermite

