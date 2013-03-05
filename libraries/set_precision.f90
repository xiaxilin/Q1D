!============================== set_precision ================================80
!
! Provides IEEE 754 compliant real kinds
!
!=============================================================================80
module set_precision

  implicit none

  private

  public :: sngl, dbl, quad, dp

  integer, parameter :: sngl = selected_real_kind( 6, 37)
  integer, parameter :: dbl  = selected_real_kind(15, 307)
  integer, parameter :: quad = selected_real_kind(33, 4931)
  integer, parameter :: dp   = dbl

! one byte integer (greater than 10e2)
  integer, parameter, public :: system_i1=selected_int_kind(2)
  integer, parameter, public :: i1=system_i1

! two byte integer (greater than 10e3)
  integer, parameter, public :: system_i2=selected_int_kind(3)
  integer, parameter, public :: i2=system_i2


end module set_precision
