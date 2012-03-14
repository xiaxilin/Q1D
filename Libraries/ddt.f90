module ddt

  use set_precision, only : dp

  implicit none

  private

  public :: ddt3
  public :: ddt4
  public :: ddt5

  type ddt3
    real(dp)               :: f ! function
    real(dp), dimension(3) :: d ! the derivatives of f
  end type ddt3

  type ddt4
    real(dp)               :: f ! function
    real(dp), dimension(4) :: d ! the derivatives of f
  end type ddt4

  type ddt5
    real(dp)               :: f ! function
    real(dp), dimension(5) :: d ! the derivatives of f
  end type ddt5

  public :: ddt3_new
  public :: ddt4_new
  public :: ddt5_new

  public :: ddt3_identity
  public :: ddt4_identity
  public :: ddt5_identity

  public :: ddt_jacobian
  interface ddt_jacobian
    module procedure ddt3_jacobian
    module procedure ddt4_jacobian
    module procedure ddt5_jacobian
  end interface

  public :: assignment(=)
  interface assignment(=)
    module procedure ddt3_assign_dr
    module procedure ddt3_assign_rd
    module procedure ddt4_assign_dr
    module procedure ddt4_assign_rd
    module procedure ddt5_assign_dr
    module procedure ddt5_assign_rd
  end interface

  public :: operator(+)
  interface operator(+)
    module procedure ddt3_plus
    module procedure ddt3_plus_dr
    module procedure ddt3_plus_rd
    module procedure ddt4_plus
    module procedure ddt4_plus_dr
    module procedure ddt4_plus_rd
    module procedure ddt5_plus
    module procedure ddt5_plus_dr
    module procedure ddt5_plus_rd
  end interface

  public :: operator(-)
  interface operator(-)
    module procedure ddt3_minus
    module procedure ddt3_minus_rd
    module procedure ddt3_minus_dr
    module procedure ddt3_minus_unary
    module procedure ddt4_minus
    module procedure ddt4_minus_rd
    module procedure ddt4_minus_dr
    module procedure ddt4_minus_unary
    module procedure ddt5_minus
    module procedure ddt5_minus_rd
    module procedure ddt5_minus_dr
    module procedure ddt5_minus_unary
  end interface

  public :: operator(*)
  interface operator(*)
    module procedure ddt3_multiply
    module procedure ddt3_multiply_dr
    module procedure ddt3_multiply_rd
    module procedure ddt4_multiply
    module procedure ddt4_multiply_dr
    module procedure ddt4_multiply_rd
    module procedure ddt5_multiply
    module procedure ddt5_multiply_dr
    module procedure ddt5_multiply_rd
  end interface

  public :: operator(/)
  interface operator(/)
    module procedure ddt3_divide
    module procedure ddt3_divide_dr
    module procedure ddt3_divide_rd
    module procedure ddt4_divide
    module procedure ddt4_divide_dr
    module procedure ddt4_divide_rd
    module procedure ddt5_divide
    module procedure ddt5_divide_dr
    module procedure ddt5_divide_rd
  end interface

  public :: operator(**)
  interface operator(**)
    module procedure ddt3_power
    module procedure ddt3_power_dr
    module procedure ddt3_power_di
    module procedure ddt4_power
    module procedure ddt4_power_dr
    module procedure ddt4_power_di
    module procedure ddt5_power
    module procedure ddt5_power_dr
    module procedure ddt5_power_di
  end interface

  public :: operator(>)
  interface operator(>)
    module procedure ddt3_greater_than
    module procedure ddt3_greater_than_dr
    module procedure ddt3_greater_than_rd
    module procedure ddt4_greater_than
    module procedure ddt4_greater_than_dr
    module procedure ddt4_greater_than_rd
    module procedure ddt5_greater_than
    module procedure ddt5_greater_than_dr
    module procedure ddt5_greater_than_rd
  end interface

  public :: operator(>=)
  interface operator(>=)
    module procedure ddt3_greater_than_equal
    module procedure ddt3_greater_than_equal_dr
    module procedure ddt4_greater_than_equal
    module procedure ddt4_greater_than_equal_dr
    module procedure ddt5_greater_than_equal
    module procedure ddt5_greater_than_equal_dr
  end interface

  public :: operator(<)
  interface operator(<)
    module procedure ddt3_less_than
    module procedure ddt3_less_than_dr
    module procedure ddt3_less_than_rd
    module procedure ddt4_less_than
    module procedure ddt4_less_than_dr
    module procedure ddt4_less_than_rd
    module procedure ddt5_less_than
    module procedure ddt5_less_than_dr
    module procedure ddt5_less_than_rd
  end interface

  public :: operator(<=)
  interface operator(<=)
    module procedure ddt3_less_than_equal
    module procedure ddt3_less_than_equal_dr
    module procedure ddt4_less_than_equal
    module procedure ddt4_less_than_equal_dr
    module procedure ddt5_less_than_equal
    module procedure ddt5_less_than_equal_dr
  end interface

  public :: operator(==)
  interface operator(==)
    module procedure ddt3_equal
    module procedure ddt3_equal_dr
    module procedure ddt4_equal
    module procedure ddt4_equal_dr
    module procedure ddt5_equal
    module procedure ddt5_equal_dr
  end interface

  public :: log
  interface log
    module procedure ddt3_log
    module procedure ddt4_log
    module procedure ddt5_log
  end interface

  public :: ddt_exp
  interface ddt_exp
    module procedure ddt3_exp
    module procedure ddt4_exp
    module procedure ddt5_exp
  end interface

  public :: ddt_tanh
  interface ddt_tanh
    module procedure ddt3_tanh
    module procedure ddt4_tanh
    module procedure ddt5_tanh
  end interface

  public :: sin
  interface sin
    module procedure ddt3_sin
    module procedure ddt4_sin
    module procedure ddt5_sin
  end interface

  public :: cos
  interface cos
    module procedure ddt3_cos
    module procedure ddt4_cos
    module procedure ddt5_cos
  end interface

  public :: sqrt
  interface sqrt
    module procedure ddt3_sqrt
    module procedure ddt4_sqrt
    module procedure ddt5_sqrt
  end interface

  public :: ddt_max
  interface ddt_max
    module procedure ddt3_max
    module procedure ddt3_max_rd
    module procedure ddt3_max_dr
    module procedure ddt4_max
    module procedure ddt4_max_rd
    module procedure ddt4_max_dr
    module procedure ddt5_max
    module procedure ddt5_max_rd
    module procedure ddt5_max_dr
  end interface

  public :: ddt_min
  interface ddt_min
    module procedure ddt3_min
    module procedure ddt3_min_rd
    module procedure ddt3_min_dr
    module procedure ddt4_min
    module procedure ddt4_min_rd
    module procedure ddt4_min_dr
    module procedure ddt5_min
    module procedure ddt5_min_rd
    module procedure ddt5_min_dr
  end interface

  public :: ddt_abs
  interface ddt_abs
    module procedure ddt3_abs
    module procedure ddt4_abs
    module procedure ddt5_abs
  end interface

  public :: ddt_sign
  interface ddt_sign
    module procedure ddt3_sign_rd
    module procedure ddt4_sign_rd
    module procedure ddt5_sign_rd
  end interface

  real(dp), parameter :: ddt_sqrt_floor = epsilon(1.0_dp)

contains

!================================= ddt3_new ==================================80
  pure elemental function ddt3_new(f, d1, d2, d3)
    real(dp), intent(in), optional  :: f, d1, d2, d3
    type(ddt3)                      :: ddt3_new

    real(dp), parameter :: zero = 0.0_dp

    continue

    if (present(f)) then
      ddt3_new%f = f
    else
      ddt3_new%f = zero
    endif

    if (present(d1)) then
      ddt3_new%d(1) = d1
    else
      ddt3_new%d(1) = zero
    endif

    if (present(d2)) then
      ddt3_new%d(2) = d2
    else
      ddt3_new%d(2) = zero
    endif

    if (present(d3)) then
      ddt3_new%d(3) = d3
    else
      ddt3_new%d(3) = zero
    endif

  end function ddt3_new

!================================= ddt4_new ==================================80
  pure elemental function ddt4_new(f, d1, d2, d3, d4)
    real(dp), intent(in), optional :: f, d1, d2, d3, d4
    type(ddt4)                     :: ddt4_new

    real(dp), parameter :: zero = 0.0_dp

    continue

    if (present(f)) then
      ddt4_new%f = f
    else
      ddt4_new%f = zero
    endif

    if (present(d1)) then
      ddt4_new%d(1) = d1
    else
      ddt4_new%d(1) = zero
    endif

    if (present(d2)) then
      ddt4_new%d(2) = d2
    else
      ddt4_new%d(2) = zero
    endif

    if (present(d3)) then
      ddt4_new%d(3) = d3
    else
      ddt4_new%d(3) = zero
    endif

    if (present(d4)) then
      ddt4_new%d(4) = d4
    else
      ddt4_new%d(4) = zero
    endif

  end function ddt4_new

!================================= ddt5_new ==================================80
  pure elemental function ddt5_new(f, d1, d2, d3, d4, d5)
    real(dp), intent(in), optional :: f, d1, d2, d3, d4, d5
    type(ddt5)                     :: ddt5_new

    real(dp), parameter :: zero = 0.0_dp

    continue

    if (present(f)) then
      ddt5_new%f = f
    else
      ddt5_new%f = zero
    endif

    if (present(d1)) then
      ddt5_new%d(1) = d1
    else
      ddt5_new%d(1) = zero
    endif

    if (present(d2)) then
      ddt5_new%d(2) = d2
    else
      ddt5_new%d(2) = zero
    endif

    if (present(d3)) then
      ddt5_new%d(3) = d3
    else
      ddt5_new%d(3) = zero
    endif

    if (present(d4)) then
      ddt5_new%d(4) = d4
    else
      ddt5_new%d(4) = zero
    endif

    if (present(d5)) then
      ddt5_new%d(5) = d5
    else
      ddt5_new%d(5) = zero
    endif

  end function ddt5_new

!================================= ddt*_identity =============================80
  pure function ddt3_identity(f)
    real(dp),   dimension(3), intent(in) :: f
    type(ddt3), dimension(3)             :: ddt3_identity

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt3_identity = f
    ddt3_identity(1)%d(1) = one
    ddt3_identity(2)%d(2) = one
    ddt3_identity(3)%d(3) = one

  end function ddt3_identity

  pure function ddt4_identity(f)
    real(dp),   dimension(4), intent(in) :: f
    type(ddt4), dimension(4)             :: ddt4_identity

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt4_identity = f
    ddt4_identity(1)%d(1) = one
    ddt4_identity(2)%d(2) = one
    ddt4_identity(3)%d(3) = one
    ddt4_identity(4)%d(4) = one

  end function ddt4_identity

  pure function ddt5_identity(f)
    real(dp),   dimension(5), intent(in) :: f
    type(ddt5), dimension(5)             :: ddt5_identity

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt5_identity = f
    ddt5_identity(1)%d(1) = one
    ddt5_identity(2)%d(2) = one
    ddt5_identity(3)%d(3) = one
    ddt5_identity(4)%d(4) = one
    ddt5_identity(5)%d(5) = one

  end function ddt5_identity

!================================= ddt*_jacobian =============================80
  pure function ddt3_jacobian(ddt)

    type(ddt3), dimension(3), intent(in) :: ddt
    real(dp),   dimension(3,3)           :: ddt3_jacobian

    integer :: row

    continue

    do row = 1, 3
      ddt3_jacobian(row,:) = ddt(row)%d
    end do

  end function ddt3_jacobian

  pure function ddt4_jacobian(ddt)

    type(ddt4), dimension(4), intent(in) :: ddt
    real(dp),   dimension(4,4)           :: ddt4_jacobian

    integer :: row

    continue

    do row = 1, 4
      ddt4_jacobian(row,:) = ddt(row)%d
    end do

  end function ddt4_jacobian

  pure function ddt5_jacobian(ddt)

    type(ddt5), dimension(5), intent(in) :: ddt
    real(dp),   dimension(5,5)           :: ddt5_jacobian

    integer :: row

    continue

    do row = 1, 5
      ddt5_jacobian(row,:) = ddt(row)%d
    end do

  end function ddt5_jacobian

!================================= ddt3_assign ===============================80
  pure elemental subroutine ddt3_assign_dr(d, r)
    type(ddt3), intent(out) :: d
    real(dp),   intent(in)  :: r

    real(dp), parameter :: zero = 0.0_dp

    continue

    d%f = r
    d%d = zero

  end subroutine ddt3_assign_dr

  pure elemental subroutine ddt3_assign_rd(r, d)
    real(dp),   intent(out) :: r
    type(ddt3), intent(in)  :: d

    continue

    r = d%f

  end subroutine ddt3_assign_rd

!================================= ddt4_assign ===============================80
  pure elemental subroutine ddt4_assign_dr(d, r)
    type(ddt4), intent(out) :: d
    real(dp),   intent(in)  :: r

    real(dp), parameter :: zero = 0.0_dp

    continue

    d%f = r
    d%d = zero

  end subroutine ddt4_assign_dr

  pure elemental subroutine ddt4_assign_rd(r, d)
    real(dp),   intent(out)  :: r
    type(ddt4), intent(in)   :: d

    continue

    r = d%f

  end subroutine ddt4_assign_rd

!================================= ddt5_assign ===============================80
  pure elemental subroutine ddt5_assign_dr(d, r)
    type(ddt5), intent(out) :: d
    real(dp),   intent(in)  :: r

    real(dp), parameter :: zero = 0.0_dp

    continue

    d%f = r
    d%d = zero

  end subroutine ddt5_assign_dr

  pure elemental subroutine ddt5_assign_rd(r, d)
    real(dp),   intent(out)  :: r
    type(ddt5), intent(in)   :: d

    continue

    r = d%f

  end subroutine ddt5_assign_rd

!================================= ddt3_plus =================================80
  pure elemental function ddt3_plus(a, b)
    type(ddt3), intent(in) :: a, b
    type(ddt3)             :: ddt3_plus

    continue

    ddt3_plus%f = a%f + b%f
    ddt3_plus%d = a%d + b%d

  end function ddt3_plus

  pure elemental function ddt3_plus_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_plus_dr

    continue

    ddt3_plus_dr%f = d%f + r
    ddt3_plus_dr%d = d%d

  end function ddt3_plus_dr

  pure elemental function ddt3_plus_rd(r, d)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_plus_rd

    continue

    ddt3_plus_rd%f = d%f + r
    ddt3_plus_rd%d = d%d

  end function ddt3_plus_rd

!================================= ddt4_plus =================================80
  pure elemental function ddt4_plus(a, b)
    type(ddt4), intent(in) :: a, b
    type(ddt4)             :: ddt4_plus

    continue

    ddt4_plus%f = a%f + b%f
    ddt4_plus%d = a%d + b%d

  end function ddt4_plus

  pure elemental function ddt4_plus_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_plus_dr

    continue

    ddt4_plus_dr%f = d%f + r
    ddt4_plus_dr%d = d%d

  end function ddt4_plus_dr

  pure elemental function ddt4_plus_rd(r, d)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_plus_rd

    continue

    ddt4_plus_rd%f = d%f + r
    ddt4_plus_rd%d = d%d

  end function ddt4_plus_rd

!================================= ddt5_plus =================================80
  pure elemental function ddt5_plus(a, b)
    type(ddt5), intent(in) :: a, b
    type(ddt5)             :: ddt5_plus

    continue

    ddt5_plus%f = a%f + b%f
    ddt5_plus%d = a%d + b%d

  end function ddt5_plus

  pure elemental function ddt5_plus_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_plus_dr

    continue

    ddt5_plus_dr%f = d%f + r
    ddt5_plus_dr%d = d%d

  end function ddt5_plus_dr

  pure elemental function ddt5_plus_rd(r, d)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_plus_rd

    continue

    ddt5_plus_rd%f = d%f + r
    ddt5_plus_rd%d = d%d

  end function ddt5_plus_rd

!================================= ddt3_minus ================================80
  pure elemental function ddt3_minus(a, b)
    type(ddt3), intent(in) :: a, b
    type(ddt3)             :: ddt3_minus

    continue

    ddt3_minus%f = a%f - b%f
    ddt3_minus%d = a%d - b%d

  end function ddt3_minus

  pure elemental function ddt3_minus_rd(r, d)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_minus_rd

    continue

    ddt3_minus_rd%f = r - d%f
    ddt3_minus_rd%d =   - d%d

  end function ddt3_minus_rd

  pure elemental function ddt3_minus_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_minus_dr

    continue

    ddt3_minus_dr%f = d%f - r
    ddt3_minus_dr%d = d%d

  end function ddt3_minus_dr

  pure elemental function ddt3_minus_unary(d)
    type(ddt3), intent(in) :: d
    type(ddt3)             :: ddt3_minus_unary

    continue

    ddt3_minus_unary%f = -d%f
    ddt3_minus_unary%d = -d%d

  end function ddt3_minus_unary

!================================= ddt4_minus ================================80
  pure elemental function ddt4_minus(a, b)
    type(ddt4), intent(in) :: a, b
    type(ddt4)             :: ddt4_minus

    continue

    ddt4_minus%f = a%f - b%f
    ddt4_minus%d = a%d - b%d

  end function ddt4_minus

  pure elemental function ddt4_minus_rd(r, d)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_minus_rd

    continue

    ddt4_minus_rd%f = r - d%f
    ddt4_minus_rd%d =   - d%d

  end function ddt4_minus_rd

  pure elemental function ddt4_minus_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_minus_dr

    continue

    ddt4_minus_dr%f = d%f - r
    ddt4_minus_dr%d = d%d

  end function ddt4_minus_dr

  pure elemental function ddt4_minus_unary(d)
    type(ddt4), intent(in) :: d
    type(ddt4)             :: ddt4_minus_unary

    continue

    ddt4_minus_unary%f = -d%f
    ddt4_minus_unary%d = -d%d

  end function ddt4_minus_unary

!================================= ddt5_minus ================================80
  pure elemental function ddt5_minus(a, b)
    type(ddt5), intent(in) :: a, b
    type(ddt5)             :: ddt5_minus

    continue

    ddt5_minus%f = a%f - b%f
    ddt5_minus%d = a%d - b%d

  end function ddt5_minus

  pure elemental function ddt5_minus_rd(r, d)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_minus_rd

    continue

    ddt5_minus_rd%f = r - d%f
    ddt5_minus_rd%d =   - d%d

  end function ddt5_minus_rd

  pure elemental function ddt5_minus_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_minus_dr

    continue

    ddt5_minus_dr%f = d%f - r
    ddt5_minus_dr%d = d%d

  end function ddt5_minus_dr

  pure elemental function ddt5_minus_unary(d)
    type(ddt5), intent(in) :: d
    type(ddt5)             :: ddt5_minus_unary

    continue

    ddt5_minus_unary%f = -d%f
    ddt5_minus_unary%d = -d%d

  end function ddt5_minus_unary

!================================= ddt3_multiply =============================80
  pure elemental function ddt3_multiply(a, b)
    type(ddt3), intent(in) :: a, b
    type(ddt3)             :: ddt3_multiply

    continue

    ddt3_multiply%f = a%f * b%f
    ddt3_multiply%d = a%f*b%d + b%f*a%d

  end function ddt3_multiply

  pure elemental function ddt3_multiply_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_multiply_dr

    continue

    ddt3_multiply_dr%f = r*d%f
    ddt3_multiply_dr%d = r*d%d

  end function ddt3_multiply_dr

  pure elemental function ddt3_multiply_rd(r, d)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_multiply_rd

    continue

    ddt3_multiply_rd%f = r*d%f
    ddt3_multiply_rd%d = r*d%d

  end function ddt3_multiply_rd

!================================= ddt4_multiply =============================80
  pure elemental function ddt4_multiply(a, b)
    type(ddt4), intent(in) :: a, b
    type(ddt4)             :: ddt4_multiply

    continue

    ddt4_multiply%f = a%f * b%f
    ddt4_multiply%d = a%f*b%d + b%f*a%d

  end function ddt4_multiply

  pure elemental function ddt4_multiply_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_multiply_dr

    continue

    ddt4_multiply_dr%f = r*d%f
    ddt4_multiply_dr%d = r*d%d

  end function ddt4_multiply_dr

  pure elemental function ddt4_multiply_rd(r, d)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_multiply_rd

    continue

    ddt4_multiply_rd%f = r*d%f
    ddt4_multiply_rd%d = r*d%d

  end function ddt4_multiply_rd

!================================= ddt5_multiply =============================80
  pure elemental function ddt5_multiply(a, b)
    type(ddt5), intent(in) :: a, b
    type(ddt5)             :: ddt5_multiply

    continue

    ddt5_multiply%f = a%f * b%f
    ddt5_multiply%d = a%f*b%d + b%f*a%d

  end function ddt5_multiply

  pure elemental function ddt5_multiply_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_multiply_dr

    continue

    ddt5_multiply_dr%f = r*d%f
    ddt5_multiply_dr%d = r*d%d

  end function ddt5_multiply_dr

  pure elemental function ddt5_multiply_rd(r, d)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_multiply_rd

    continue

    ddt5_multiply_rd%f = r*d%f
    ddt5_multiply_rd%d = r*d%d

  end function ddt5_multiply_rd

!================================= ddt3_divide ===============================80
  pure elemental function ddt3_divide(a, b)
    type(ddt3), intent(in) :: a, b
    type(ddt3)             :: ddt3_divide

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt3_divide%f = a%f / b%f

    denom = one / b%f / b%f

    ddt3_divide%d = ( a%d*b%f - b%d*a%f ) * denom

  end function ddt3_divide

  pure elemental function ddt3_divide_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_divide_dr

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt3_divide_dr%f = d%f / r

    denom = one / r / r

    ddt3_divide_dr%d = ( d%d*r ) * denom

  end function ddt3_divide_dr

  pure elemental function ddt3_divide_rd(r, d)
    real(dp),   intent(in) :: r
    type(ddt3), intent(in) :: d
    type(ddt3)             :: ddt3_divide_rd

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt3_divide_rd%f = r / d%f

    denom = one / d%f / d%f

    ddt3_divide_rd%d = ( - d%d*r ) * denom

  end function ddt3_divide_rd

!================================= ddt4_divide ===============================80
  pure elemental function ddt4_divide(a, b)
    type(ddt4), intent(in) :: a, b
    type(ddt4)             :: ddt4_divide

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt4_divide%f = a%f / b%f

    denom = one / b%f / b%f

    ddt4_divide%d = ( a%d*b%f - b%d*a%f ) * denom

  end function ddt4_divide

  pure elemental function ddt4_divide_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_divide_dr

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt4_divide_dr%f = d%f / r

    denom = one / r / r

    ddt4_divide_dr%d = ( d%d*r ) * denom

  end function ddt4_divide_dr

  pure elemental function ddt4_divide_rd(r, d)
    real(dp),   intent(in) :: r
    type(ddt4), intent(in) :: d
    type(ddt4)             :: ddt4_divide_rd

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt4_divide_rd%f = r / d%f

    denom = one / d%f / d%f

    ddt4_divide_rd%d = ( - d%d*r ) * denom

  end function ddt4_divide_rd

!================================= ddt5_divide ===============================80
  pure elemental function ddt5_divide(a, b)
    type(ddt5), intent(in) :: a, b
    type(ddt5)             :: ddt5_divide

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt5_divide%f = a%f / b%f

    denom = one / b%f / b%f

    ddt5_divide%d = ( a%d*b%f - b%d*a%f ) * denom

  end function ddt5_divide

  pure elemental function ddt5_divide_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_divide_dr

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt5_divide_dr%f = d%f / r

    denom = one / r / r

    ddt5_divide_dr%d = ( d%d*r ) * denom

  end function ddt5_divide_dr

  pure elemental function ddt5_divide_rd(r, d)
    real(dp),   intent(in) :: r
    type(ddt5), intent(in) :: d
    type(ddt5)             :: ddt5_divide_rd

    real(dp)            :: denom
    real(dp), parameter :: one = 1.0_dp

    continue

    ddt5_divide_rd%f = r / d%f

    denom = one / d%f / d%f

    ddt5_divide_rd%d = ( - d%d*r ) * denom

  end function ddt5_divide_rd

!================================= ddt3_power ================================80
  pure elemental function ddt3_power(a, b)
    type(ddt3), intent(in) :: a, b
    type(ddt3)             :: ddt3_power

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt3_power%f = a%f**b%f
    ddt3_power%d = b%f * a%f**(b%f-one) * a%d + log(a%f) * a%f**b%f * b%d

  end function ddt3_power

  pure elemental function ddt3_power_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt3)             :: ddt3_power_dr

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt3_power_dr%f = d%f**r
    ddt3_power_dr%d = r * d%f**(r-one) * d%d

  end function ddt3_power_dr

  pure elemental function ddt3_power_di(d, i)
    type(ddt3), intent(in) :: d
    integer,    intent(in) :: i
    type(ddt3)             :: ddt3_power_di

    continue

    ddt3_power_di%f = d%f**i
    ddt3_power_di%d = real(i,dp) * d%f**(i-1) * d%d

  end function ddt3_power_di

!================================= ddt4_power ================================80
  pure elemental function ddt4_power(a, b)
    type(ddt4), intent(in) :: a, b
    type(ddt4)             :: ddt4_power

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt4_power%f = a%f**b%f
    ddt4_power%d = b%f * a%f**(b%f-one) * a%d + log(a%f) * a%f**b%f * b%d

  end function ddt4_power

  pure elemental function ddt4_power_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt4)             :: ddt4_power_dr

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt4_power_dr%f = d%f**r
    ddt4_power_dr%d = r * d%f**(r-one) * d%d

  end function ddt4_power_dr

  pure elemental function ddt4_power_di(d, i)
    type(ddt4), intent(in) :: d
    integer,    intent(in) :: i
    type(ddt4)             :: ddt4_power_di

    continue

    ddt4_power_di%f = d%f**i
    ddt4_power_di%d = real(i,dp) * d%f**(i-1) * d%d

  end function ddt4_power_di

!================================= ddt5_power ================================80
  pure elemental function ddt5_power(a, b)
    type(ddt5), intent(in) :: a, b
    type(ddt5)             :: ddt5_power

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt5_power%f = a%f**b%f
    ddt5_power%d = b%f * a%f**(b%f-one) * a%d + log(a%f) * a%f**b%f * b%d

  end function ddt5_power

  pure elemental function ddt5_power_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    type(ddt5)             :: ddt5_power_dr

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt5_power_dr%f = d%f**r
    ddt5_power_dr%d = r * d%f**(r-one) * d%d

  end function ddt5_power_dr

  pure elemental function ddt5_power_di(d, i)
    type(ddt5), intent(in) :: d
    integer,    intent(in) :: i
    type(ddt5)             :: ddt5_power_di

    continue

    ddt5_power_di%f = d%f**i
    ddt5_power_di%d = real(i,dp) * d%f**(i-1) * d%d

  end function ddt5_power_di

!================================= ddt3_greater_than =========================80
  pure elemental function ddt3_greater_than(a, b)
    type(ddt3), intent(in) :: a, b
    logical                :: ddt3_greater_than

    continue

    ddt3_greater_than = ( a%f > b%f )

  end function ddt3_greater_than

  pure elemental function ddt3_greater_than_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt3_greater_than_dr

    continue

    ddt3_greater_than_dr = ( d%f > r )

  end function ddt3_greater_than_dr

  pure elemental function ddt3_greater_than_rd(r, d)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt3_greater_than_rd

    continue

    ddt3_greater_than_rd = ( r > d%f )

  end function ddt3_greater_than_rd

!================================= ddt4_greater_than =========================80
  pure elemental function ddt4_greater_than(a, b)
    type(ddt4), intent(in) :: a, b
    logical                :: ddt4_greater_than

    continue

    ddt4_greater_than = ( a%f > b%f )

  end function ddt4_greater_than

  pure elemental function ddt4_greater_than_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt4_greater_than_dr

    continue

    ddt4_greater_than_dr = ( d%f > r )

  end function ddt4_greater_than_dr

  pure elemental function ddt4_greater_than_rd(r, d)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt4_greater_than_rd

    continue

    ddt4_greater_than_rd = ( r > d%f )

  end function ddt4_greater_than_rd

!================================= ddt5_greater_than =========================80
  pure elemental function ddt5_greater_than(a, b)
    type(ddt5), intent(in) :: a, b
    logical                :: ddt5_greater_than

    continue

    ddt5_greater_than = ( a%f > b%f )

  end function ddt5_greater_than

  pure elemental function ddt5_greater_than_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt5_greater_than_dr

    continue

    ddt5_greater_than_dr = ( d%f > r )

  end function ddt5_greater_than_dr

  pure elemental function ddt5_greater_than_rd(r, d)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt5_greater_than_rd

    continue

    ddt5_greater_than_rd = ( r > d%f )

  end function ddt5_greater_than_rd

!================================= ddt3_greater_than_equal ===================80
  pure elemental function ddt3_greater_than_equal(a, b)
    type(ddt3), intent(in) :: a, b
    logical                :: ddt3_greater_than_equal

    continue

    ddt3_greater_than_equal = ( a%f >= b%f )

  end function ddt3_greater_than_equal

  pure elemental function ddt3_greater_than_equal_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt3_greater_than_equal_dr

    continue

    ddt3_greater_than_equal_dr = ( d%f >= r )

  end function ddt3_greater_than_equal_dr

!================================= ddt4_greater_than_equal ===================80
  pure elemental function ddt4_greater_than_equal(a, b)
    type(ddt4), intent(in) :: a, b
    logical                :: ddt4_greater_than_equal

    continue

    ddt4_greater_than_equal = ( a%f >= b%f )

  end function ddt4_greater_than_equal

  pure elemental function ddt4_greater_than_equal_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt4_greater_than_equal_dr

    continue

    ddt4_greater_than_equal_dr = ( d%f >= r )

  end function ddt4_greater_than_equal_dr

!================================= ddt5_greater_than_equal ===================80
  pure elemental function ddt5_greater_than_equal(a, b)
    type(ddt5), intent(in) :: a, b
    logical                :: ddt5_greater_than_equal

    continue

    ddt5_greater_than_equal = ( a%f >= b%f )

  end function ddt5_greater_than_equal

  pure elemental function ddt5_greater_than_equal_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt5_greater_than_equal_dr

    continue

    ddt5_greater_than_equal_dr = ( d%f >= r )

  end function ddt5_greater_than_equal_dr

!================================= ddt3_less_than ============================80
  pure elemental function ddt3_less_than(a, b)
    type(ddt3), intent(in) :: a, b
    logical                :: ddt3_less_than

    continue

    ddt3_less_than = ( a%f < b%f )

  end function ddt3_less_than

  pure elemental function ddt3_less_than_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt3_less_than_dr

    continue

    ddt3_less_than_dr = ( d%f < r )

  end function ddt3_less_than_dr

  pure elemental function ddt3_less_than_rd(r, d)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt3_less_than_rd

    continue

    ddt3_less_than_rd = ( r < d%f )

  end function ddt3_less_than_rd

!================================= ddt4_less_than ============================80
  pure elemental function ddt4_less_than(a, b)
    type(ddt4), intent(in) :: a, b
    logical                :: ddt4_less_than

    continue

    ddt4_less_than = ( a%f < b%f )

  end function ddt4_less_than

  pure elemental function ddt4_less_than_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt4_less_than_dr

    continue

    ddt4_less_than_dr = ( d%f < r )

  end function ddt4_less_than_dr

  pure elemental function ddt4_less_than_rd(r, d)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt4_less_than_rd

    continue

    ddt4_less_than_rd = ( r < d%f )

  end function ddt4_less_than_rd

!================================= ddt5_less_than ============================80
  pure elemental function ddt5_less_than(a, b)
    type(ddt5), intent(in) :: a, b
    logical                :: ddt5_less_than

    continue

    ddt5_less_than = ( a%f < b%f )

  end function ddt5_less_than

  pure elemental function ddt5_less_than_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt5_less_than_dr

    continue

    ddt5_less_than_dr = ( d%f < r )

  end function ddt5_less_than_dr

  pure elemental function ddt5_less_than_rd(r, d)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt5_less_than_rd

    continue

    ddt5_less_than_rd = ( r < d%f )

  end function ddt5_less_than_rd

!================================= ddt3_less_than_equal ======================80
  pure elemental function ddt3_less_than_equal(a, b)
    type(ddt3), intent(in) :: a, b
    logical                :: ddt3_less_than_equal

    continue

    ddt3_less_than_equal = ( a%f <= b%f )

  end function ddt3_less_than_equal

  pure elemental function ddt3_less_than_equal_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt3_less_than_equal_dr

    continue

    ddt3_less_than_equal_dr = ( d%f <= r )

  end function ddt3_less_than_equal_dr

!================================= ddt4_less_than_equal ======================80
  pure elemental function ddt4_less_than_equal(a, b)
    type(ddt4), intent(in) :: a, b
    logical                :: ddt4_less_than_equal

    continue

    ddt4_less_than_equal = ( a%f <= b%f )

  end function ddt4_less_than_equal

  pure elemental function ddt4_less_than_equal_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt4_less_than_equal_dr

    continue

    ddt4_less_than_equal_dr = ( d%f <= r )

  end function ddt4_less_than_equal_dr

!================================= ddt5_less_than_equal ======================80
  pure elemental function ddt5_less_than_equal(a, b)
    type(ddt5), intent(in) :: a, b
    logical                :: ddt5_less_than_equal

    continue

    ddt5_less_than_equal = ( a%f <= b%f )

  end function ddt5_less_than_equal

  pure elemental function ddt5_less_than_equal_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt5_less_than_equal_dr

    continue

    ddt5_less_than_equal_dr = ( d%f <= r )

  end function ddt5_less_than_equal_dr

!================================= ddt3_equal ================================80
  pure elemental function ddt3_equal(a, b)
    type(ddt3), intent(in) :: a, b
    logical                :: ddt3_equal

    continue

    ddt3_equal = ( a%f == b%f )

  end function ddt3_equal

  pure elemental function ddt3_equal_dr(d, r)
    type(ddt3), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt3_equal_dr

    continue

    ddt3_equal_dr = ( d%f == r )

  end function ddt3_equal_dr

!================================= ddt4_equal ================================80
  pure elemental function ddt4_equal(a, b)
    type(ddt4), intent(in) :: a, b
    logical                :: ddt4_equal

    continue

    ddt4_equal = ( a%f == b%f )

  end function ddt4_equal

  pure elemental function ddt4_equal_dr(d, r)
    type(ddt4), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt4_equal_dr

    continue

    ddt4_equal_dr = ( d%f == r )

  end function ddt4_equal_dr

!================================= ddt5_equal ================================80
  pure elemental function ddt5_equal(a, b)
    type(ddt5), intent(in) :: a, b
    logical                :: ddt5_equal

    continue

    ddt5_equal = ( a%f == b%f )

  end function ddt5_equal

  pure elemental function ddt5_equal_dr(d, r)
    type(ddt5), intent(in) :: d
    real(dp),   intent(in) :: r
    logical                :: ddt5_equal_dr

    continue

    ddt5_equal_dr = ( d%f == r )

  end function ddt5_equal_dr

!================================= ddt3_exp ==================================80
  pure elemental function ddt3_exp(a)
    type(ddt3), intent(in) :: a
    type(ddt3)             :: ddt3_exp

    continue

    ddt3_exp%f = exp(a%f)
    ddt3_exp%d = ddt3_exp%f * a%d

  end function ddt3_exp

!================================= ddt4_exp ==================================80
  pure elemental function ddt4_exp(a)
    type(ddt4), intent(in) :: a
    type(ddt4)             :: ddt4_exp

    continue

    ddt4_exp%f = exp(a%f)
    ddt4_exp%d = ddt4_exp%f * a%d

  end function ddt4_exp

!================================= ddt5_exp ==================================80
  pure elemental function ddt5_exp(a)
    type(ddt5), intent(in) :: a
    type(ddt5)             :: ddt5_exp

    continue

    ddt5_exp%f = exp(a%f)
    ddt5_exp%d = ddt5_exp%f * a%d

  end function ddt5_exp

!================================= ddt3_log ==================================80
  pure elemental function ddt3_log(a)
    type(ddt3), intent(in) :: a
    type(ddt3)             :: ddt3_log

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt3_log%f = log(a%f)
    ddt3_log%d = ( one / a%f ) * a%d

  end function ddt3_log

!================================= ddt4_log ==================================80
  pure elemental function ddt4_log(a)
    type(ddt4), intent(in) :: a
    type(ddt4)             :: ddt4_log

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt4_log%f = log(a%f)
    ddt4_log%d = ( one / a%f ) * a%d

  end function ddt4_log

!================================= ddt5_log ==================================80
  pure elemental function ddt5_log(a)
    type(ddt5), intent(in) :: a
    type(ddt5)             :: ddt5_log

    real(dp), parameter :: one = 1.0_dp

    continue

    ddt5_log%f = log(a%f)
    ddt5_log%d = ( one / a%f ) * a%d

  end function ddt5_log

!================================= ddt3_tanh =================================80
  pure elemental function ddt3_tanh(a)
    type(ddt3), intent(in) :: a
    type(ddt3)             :: ddt3_tanh

    real(dp)            :: tanh_f
    real(dp), parameter :: one = 1.0_dp

    continue

    tanh_f = tanh(a%f)
    ddt3_tanh%f = tanh_f
    ddt3_tanh%d = ( one - tanh_f**2 ) * a%d

  end function ddt3_tanh

!================================= ddt4_tanh =================================80
  pure elemental function ddt4_tanh(a)
    type(ddt4), intent(in) :: a
    type(ddt4)             :: ddt4_tanh

    real(dp)            :: tanh_f
    real(dp), parameter :: one = 1.0_dp

    continue

    tanh_f = tanh(a%f)
    ddt4_tanh%f = tanh_f
    ddt4_tanh%d = ( one - tanh_f**2 ) * a%d

  end function ddt4_tanh

!================================= ddt5_tanh =================================80
  pure elemental function ddt5_tanh(a)
    type(ddt5), intent(in) :: a
    type(ddt5)             :: ddt5_tanh

    real(dp)            :: tanh_f
    real(dp), parameter :: one = 1.0_dp

    continue

    tanh_f = tanh(a%f)
    ddt5_tanh%f = tanh_f
    ddt5_tanh%d = ( one - tanh_f**2 ) * a%d

  end function ddt5_tanh

!================================= ddt3_cos ==================================80
  pure elemental function ddt3_cos(a)
    type(ddt3), intent(in) :: a
    type(ddt3)             :: ddt3_cos

    continue

    ddt3_cos%f =  cos(a%f)
    ddt3_cos%d = -sin(a%f) * a%d

  end function ddt3_cos

!================================= ddt4_cos ==================================80
  pure elemental function ddt4_cos(a)
    type(ddt4), intent(in) :: a
    type(ddt4)             :: ddt4_cos

    continue

    ddt4_cos%f =  cos(a%f)
    ddt4_cos%d = -sin(a%f) * a%d

  end function ddt4_cos

!================================= ddt5_cos ==================================80
  pure elemental function ddt5_cos(a)
    type(ddt5), intent(in) :: a
    type(ddt5)             :: ddt5_cos

    continue

    ddt5_cos%f =  cos(a%f)
    ddt5_cos%d = -sin(a%f) * a%d

  end function ddt5_cos

!================================= ddt3_sin ==================================80
  pure elemental function ddt3_sin(a)
    type(ddt3), intent(in) :: a
    type(ddt3)             :: ddt3_sin

    continue

    ddt3_sin%f = sin(a%f)
    ddt3_sin%d = cos(a%f) * a%d

  end function ddt3_sin

!================================= ddt4_sin ==================================80
  pure elemental function ddt4_sin(a)
    type(ddt4), intent(in) :: a
    type(ddt4)             :: ddt4_sin

    continue

    ddt4_sin%f = sin(a%f)
    ddt4_sin%d = cos(a%f) * a%d

  end function ddt4_sin

!================================= ddt5_sin ==================================80
  pure elemental function ddt5_sin(a)
    type(ddt5), intent(in) :: a
    type(ddt5)             :: ddt5_sin

    continue

    ddt5_sin%f = sin(a%f)
    ddt5_sin%d = cos(a%f) * a%d

  end function ddt5_sin

!================================= ddt3_sqrt =================================80
  pure elemental function ddt3_sqrt(a)
    type(ddt3), intent(in) :: a
    type(ddt3)             :: ddt3_sqrt

    real(dp)            :: sqrt_f
    real(dp), parameter :: half = 0.5_dp

    continue

    sqrt_f = sqrt(a%f)

    ddt3_sqrt%f = sqrt_f

    prevent_inf_derivative : if ( sqrt_f > ddt_sqrt_floor ) then
      ddt3_sqrt%d = ( half / sqrt_f ) * a%d
    else
      ddt3_sqrt%d = ( half / ddt_sqrt_floor ) * a%d
    end if prevent_inf_derivative

  end function ddt3_sqrt

!================================= ddt4_sqrt =================================80
  pure elemental function ddt4_sqrt(a)
    type(ddt4), intent(in) :: a
    type(ddt4)             :: ddt4_sqrt

    real(dp)            :: sqrt_f
    real(dp), parameter :: half = 0.5_dp

    continue

    sqrt_f = sqrt(a%f)

    ddt4_sqrt%f = sqrt_f

    prevent_inf_derivative : if ( sqrt_f > ddt_sqrt_floor ) then
      ddt4_sqrt%d = ( half / sqrt_f ) * a%d
    else
      ddt4_sqrt%d = ( half / ddt_sqrt_floor ) * a%d
    end if prevent_inf_derivative

  end function ddt4_sqrt

!================================= ddt5_sqrt =================================80
  pure elemental function ddt5_sqrt(a)
    type(ddt5), intent(in) :: a
    type(ddt5)             :: ddt5_sqrt

    real(dp)            :: sqrt_f
    real(dp), parameter :: half = 0.5_dp

    continue

    sqrt_f = sqrt(a%f)

    ddt5_sqrt%f = sqrt_f

    prevent_inf_derivative : if ( sqrt_f > ddt_sqrt_floor ) then
      ddt5_sqrt%d = ( half / sqrt_f ) * a%d
    else
      ddt5_sqrt%d = ( half / ddt_sqrt_floor ) * a%d
    end if prevent_inf_derivative

  end function ddt5_sqrt

!================================= ddt3_max ==================================80
  pure elemental function ddt3_max(a,b)
    type(ddt3), intent(in) :: a, b
    type(ddt3)             :: ddt3_max

    continue

    if ( a > b ) then
      ddt3_max = a
    else
      ddt3_max = b
    endif

  end function ddt3_max

  pure elemental function ddt3_max_rd(r,b)
    real(dp),   intent(in) :: r
    type(ddt3), intent(in) :: b
    type(ddt3)             :: ddt3_max_rd

    continue

    if ( r > b ) then
      ddt3_max_rd = r
    else
      ddt3_max_rd = b
    endif

  end function ddt3_max_rd

  pure elemental function ddt3_max_dr(d,r)
    real(dp),   intent(in) :: r
    type(ddt3), intent(in) :: d
    type(ddt3)             :: ddt3_max_dr

    continue

    if ( r > d ) then
      ddt3_max_dr = r
    else
      ddt3_max_dr = d
    endif

  end function ddt3_max_dr

!================================= ddt4_max ==================================80
  pure elemental function ddt4_max(a,b)
    type(ddt4), intent(in) :: a, b
    type(ddt4)             :: ddt4_max

    continue

    if ( a > b ) then
      ddt4_max = a
    else
      ddt4_max = b
    endif

  end function ddt4_max

  pure elemental function ddt4_max_rd(r,b)
    real(dp),   intent(in) :: r
    type(ddt4), intent(in) :: b
    type(ddt4)             :: ddt4_max_rd

    continue

    if ( r > b ) then
      ddt4_max_rd = r
    else
      ddt4_max_rd = b
    endif

  end function ddt4_max_rd

  pure elemental function ddt4_max_dr(d,r)
    real(dp),   intent(in) :: r
    type(ddt4), intent(in) :: d
    type(ddt4)             :: ddt4_max_dr

    continue

    if ( r > d ) then
      ddt4_max_dr = r
    else
      ddt4_max_dr = d
    endif

  end function ddt4_max_dr

!================================= ddt5_max ==================================80
  pure elemental function ddt5_max(a,b)
    type(ddt5), intent(in) :: a, b
    type(ddt5)             :: ddt5_max

    continue

    if ( a > b ) then
      ddt5_max = a
    else
      ddt5_max = b
    endif

  end function ddt5_max

  pure elemental function ddt5_max_rd(r,b)
    real(dp),   intent(in) :: r
    type(ddt5), intent(in) :: b
    type(ddt5)             :: ddt5_max_rd

    continue

    if ( r > b ) then
      ddt5_max_rd = r
    else
      ddt5_max_rd = b
    endif

  end function ddt5_max_rd

  pure elemental function ddt5_max_dr(d,r)
    real(dp),   intent(in) :: r
    type(ddt5), intent(in) :: d
    type(ddt5)             :: ddt5_max_dr

    continue

    if ( r > d ) then
      ddt5_max_dr = r
    else
      ddt5_max_dr = d
    endif

  end function ddt5_max_dr

!================================= ddt3_min ==================================80
  pure elemental function ddt3_min(a,b)
    type(ddt3), intent(in) :: a, b
    type(ddt3)             :: ddt3_min

    continue

    if ( a < b ) then
      ddt3_min = a
    else
      ddt3_min = b
    endif

  end function ddt3_min

  pure elemental function ddt3_min_rd(r,b)
    real(dp),   intent(in) :: r
    type(ddt3), intent(in) :: b
    type(ddt3)             :: ddt3_min_rd

    continue

    if ( r < b ) then
      ddt3_min_rd = r
    else
      ddt3_min_rd = b
    endif

  end function ddt3_min_rd

  pure elemental function ddt3_min_dr(d,r)
    real(dp),   intent(in) :: r
    type(ddt3), intent(in) :: d
    type(ddt3)             :: ddt3_min_dr

    continue

    if ( r < d ) then
      ddt3_min_dr = r
    else
      ddt3_min_dr = d
    endif

  end function ddt3_min_dr

!================================= ddt4_min ==================================80
  pure elemental function ddt4_min(a,b)
    type(ddt4), intent(in) :: a, b
    type(ddt4)             :: ddt4_min

    continue

    if ( a < b ) then
      ddt4_min = a
    else
      ddt4_min = b
    endif

  end function ddt4_min

  pure elemental function ddt4_min_rd(r,b)
    real(dp),   intent(in) :: r
    type(ddt4), intent(in) :: b
    type(ddt4)             :: ddt4_min_rd

    continue

    if ( r < b ) then
      ddt4_min_rd = r
    else
      ddt4_min_rd = b
    endif

  end function ddt4_min_rd

  pure elemental function ddt4_min_dr(d,r)
    real(dp),   intent(in) :: r
    type(ddt4), intent(in) :: d
    type(ddt4)             :: ddt4_min_dr

    continue

    if ( r < d ) then
      ddt4_min_dr = r
    else
      ddt4_min_dr = d
    endif

  end function ddt4_min_dr

!================================= ddt5_min ==================================80
  pure elemental function ddt5_min(a,b)
    type(ddt5), intent(in) :: a, b
    type(ddt5)             :: ddt5_min

    continue

    if ( a < b ) then
      ddt5_min = a
    else
      ddt5_min = b
    endif

  end function ddt5_min

  pure elemental function ddt5_min_rd(r,b)
    real(dp),   intent(in) :: r
    type(ddt5), intent(in) :: b
    type(ddt5)             :: ddt5_min_rd

    continue

    if ( r < b ) then
      ddt5_min_rd = r
    else
      ddt5_min_rd = b
    endif

  end function ddt5_min_rd

  pure elemental function ddt5_min_dr(d,r)
    real(dp),   intent(in) :: r
    type(ddt5), intent(in) :: d
    type(ddt5)             :: ddt5_min_dr

    continue

    if ( r < d ) then
      ddt5_min_dr = r
    else
      ddt5_min_dr = d
    endif

  end function ddt5_min_dr

!================================= ddt3_abs ==================================80
  pure elemental function ddt3_abs(a)
    type(ddt3), intent(in) :: a
    type(ddt3)             :: ddt3_abs

    real(dp), parameter :: zero = 0.0_dp

    continue

    if ( a < zero ) then
      ddt3_abs = -a
    else
      ddt3_abs = a
    endif

  end function ddt3_abs

!================================= ddt4_abs ==================================80
  pure elemental function ddt4_abs(a)
    type(ddt4), intent(in) :: a
    type(ddt4)             :: ddt4_abs

    real(dp), parameter :: zero = 0.0_dp

    continue

    if ( a < zero ) then
      ddt4_abs = -a
    else
      ddt4_abs = a
    endif

  end function ddt4_abs

!================================= ddt5_abs ==================================80
  pure elemental function ddt5_abs(a)
    type(ddt5), intent(in) :: a
    type(ddt5)             :: ddt5_abs

    real(dp), parameter :: zero = 0.0_dp

    continue

    if ( a < zero ) then
      ddt5_abs = -a
    else
      ddt5_abs = a
    endif

  end function ddt5_abs

!================================= ddt3_sign_rd ==============================80
  pure elemental function ddt3_sign_rd(a,b)
    real(dp),   intent(in) :: a
    type(ddt3), intent(in) :: b
    type(ddt3)             :: ddt3_sign_rd

    real(dp), parameter :: zero = 0.0_dp

    continue

    if ( b >= zero) then
      ddt3_sign_rd =  abs(a)
    else
      ddt3_sign_rd = -abs(a)
    end if

  end function ddt3_sign_rd

!================================= ddt4_sign_rd ==============================80
  pure elemental function ddt4_sign_rd(a,b)
    real(dp),   intent(in) :: a
    type(ddt4), intent(in) :: b
    type(ddt4)             :: ddt4_sign_rd

    real(dp), parameter :: zero = 0.0_dp

    continue

    if ( b >= zero) then
      ddt4_sign_rd =  abs(a)
    else
      ddt4_sign_rd = -abs(a)
    end if

  end function ddt4_sign_rd

!================================= ddt5_sign_rd ==============================80
  pure elemental function ddt5_sign_rd(a,b)
    real(dp),   intent(in) :: a
    type(ddt5), intent(in) :: b
    type(ddt5)             :: ddt5_sign_rd

    real(dp), parameter :: zero = 0.0_dp

    continue

    if ( b >= zero) then
      ddt5_sign_rd =  abs(a)
    else
      ddt5_sign_rd = -abs(a)
    end if

  end function ddt5_sign_rd

end module ddt
