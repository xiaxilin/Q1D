!FIXME: change initial conditions to be smarter...
!ie, find throat if it exists, set pback in outflow cell if > 0, etc

module initialize_soln

  use set_precision, only : dp

  implicit none

  private

  public :: allocate_soln
  public :: deallocate_soln
  public :: initial_soln

  public :: prim_cc
  public :: cons_cc

  public :: L1_init
  public :: L2_init
  public :: Linf_init

  public :: restart
  public :: mref       ! Initial mach number in nozzle
  public :: to         ! Inflow stag. temp
  public :: po         ! Inflow stag. pressure
  public :: pback      ! Outflow backpressure, negative for extrapolation

! Set defaults
  logical :: restart

  real(dp) :: mref
  real(dp) :: to
  real(dp) :: po
  real(dp) :: pback

  real(dp), dimension(3) :: L1_init, L2_init, Linf_init

  real(dp), allocatable, dimension(:,:) :: prim_cc   ! primitive vars
  real(dp), allocatable, dimension(:,:) :: cons_cc   ! conserved vars

contains

!=============================== allocate_soln ===============================80
!
! Allocates the soln variables
!
!=============================================================================80
  subroutine allocate_soln(neq, cells)

    implicit none

    integer, intent(in) :: neq, cells

    continue

    allocate( prim_cc(neq, cells+2), cons_cc(neq, cells+2) )
!    allocate( prim_cc(neq, cells+2) )

  end subroutine allocate_soln

!============================== deallocate_soln ==============================80
!
! Deallocates the soln variables in the correct order
!
!=============================================================================80
  subroutine deallocate_soln()

    implicit none

    continue

    deallocate( cons_cc, prim_cc )
!    deallocate( prim_cc )

  end subroutine deallocate_soln

!================================ initial_soln ===============================80
!
! Sets the initial soln
!
!=============================================================================80
  subroutine initial_soln(cells)

    use set_constants,   only : half, one
    use fluid_constants, only : r, gamma, gm1, gxgm1

    implicit none

    integer, intent(in) :: cells

    integer  :: cell, restart_unit
    real(dp) :: psi, t, p, m

    continue

    L1_init(:)   = one
    L2_init(:)   = one
    Linf_init(:) = one

    if (restart) then

      restart_unit = find_available_unit()

      open(restart_unit, file='q1d.rst', status='old')

      read(restart_unit,*) L1_init(1), L1_init(2), L1_init(3)
      read(restart_unit,*) L2_init(1), L2_init(2), L2_init(3)
      read(restart_unit,*) Linf_init(1), Linf_init(2), Linf_init(3)

      do cell = 1, cells+2
        read(restart_unit,*) prim_cc(1, cell), prim_cc(2, cell), prim_cc(3,cell)
        cons_cc(:,cell) = primitive_to_conserved_1D( prim_cc(:,cell) )
      end do

      close(restart_unit)

    else

! Constant initial conditions
!      psi = one + half*gm1*mref*mref
!      t   = to/psi
!      p   = po/(psi**gxgm1)

!      do cell = 1, cells+2
!        prim_cc(1,cell) = p/(r*t)
!        prim_cc(2,cell) = mref*sqrt(gamma*r*t)
!        prim_cc(3,cell) = p
!      end do

! Inflow velocity set, everything else stagnant
!      psi = one + half*gm1*mref*mref
!      t   = to/psi
!      p   = po/(psi**gxgm1)

!      prim_cc(1,1) = p/(r*t)
!      prim_cc(2,1) = mref*sqrt(gamma*r*t)
!      prim_cc(3,1) = p

!      do cell = 2, cells+2
!        prim_cc(1,cell) = po/(r*to)
!        prim_cc(2,cell) = zero
!        prim_cc(3,cell) = po
!      end do

! Linear ramp from mref to Mach = 1 at the throat
      do cell = 1, cells+2
        m = mref + (one-mref)*real(cell,dp)/real(cells/2+1,dp)
        psi = one + half*gm1*m**2
        t   = to/psi
        p   = po/(psi**gxgm1)

        prim_cc(1,cell) = p/(r*t)
        prim_cc(2,cell) = m*sqrt(gamma*r*t)
        prim_cc(3,cell) = p
      end do

      do cell = 1, cells+2
        cons_cc(:,cell) = primitive_to_conserved_1D( prim_cc(:,cell) )
      end do

    endif

  end subroutine initial_soln

  include 'find_available_unit.f90'
  include 'primitive_to_conserved_1D.f90'

end module initialize_soln
