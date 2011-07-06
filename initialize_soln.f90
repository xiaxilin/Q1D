module initialize_soln

  use set_precision, only : dp

  implicit none

  private

  public :: allocate_soln
  public :: deallocate_soln
  public :: initial_soln

  public :: prim_cc
  public :: cons_cc

  public :: restart
  public :: mref       ! Initial mach number in nozzle
  public :: to         ! Inflow stag. temp
  public :: po         ! Inflow stag. pressure
  public :: pback      ! Outflow backpressure, negative for extrapolation

! Set defaults
  logical :: restart = .false.

  real(dp) :: mref  = 1.5_dp
  real(dp) :: to    = 600.0_dp
  real(dp) :: po    = 300000.0_dp
  real(dp) :: pback = -1.0_dp

  real(dp), allocatable, dimension(:,:) :: prim_cc   ! primitive vars
  real(dp), allocatable, dimension(:,:) :: cons_cc   ! conserved vars

contains

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine allocate_soln(neq, cells)

    implicit none

    integer, intent(in) :: neq, cells

    continue

    allocate( prim_cc(neq, cells+2), cons_cc(neq, cells+2) )

  end subroutine allocate_soln

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine deallocate_soln()

    implicit none

    continue

    deallocate( prim_cc, cons_cc )

  end subroutine deallocate_soln

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine initial_soln(cells)

    use set_constants,   only : half, one
    use fluid_constants, only : r, gamma, gm1, gxgm1

    implicit none

    integer, intent(in) :: cells

    integer  :: cell, soln_unit
    real(dp) :: p, psi, t

    continue

    if (restart) then

      soln_unit = find_available_unit()

      open(soln_unit, file='q1d.sln', status='old')

      do cell = 1, cells+2
        read(soln_unit,*) prim_cc(1, cell), prim_cc(2, cell), prim_cc(3,cell)
        cons_cc(:,cell) = primitive_to_conserved_1D( prim_cc(:,cell) )
      end do

      close(soln_unit)

    else

      psi = one + half*gm1*mref*mref
      t   = to/psi
      p = po/(psi**gxgm1)

      do cell = 1, cells+2
        prim_cc(1,cell) = p/(r*t)
        prim_cc(2,cell) = mref*sqrt(gamma*r*t)
        prim_cc(3,cell) = p

        cons_cc(:,cell) = primitive_to_conserved_1D( prim_cc(:,cell) )
      end do

    endif

  end subroutine initial_soln

  include 'find_available_unit.f90'
  include 'primitive_to_conserved_1D.f90'

end module initialize_soln
