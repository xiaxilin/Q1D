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

end module set_precision
