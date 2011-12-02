module namelist

  use solvers,         only : iterations, firstorder, itercheck, iter_out,     &
                              rkorder, cfl, limiter, muscl, kappa, toler,      &
                              flux_type, k2, k4
  use initialize_soln, only : restart, mref, to, po, pback
  use fluid_constants, only : gamma, r

  implicit none

  private

  public :: read_nml

! Define namelist inputs

  namelist /code_control/ iterations, firstorder, itercheck, iter_out,         &
                          rkorder, cfl, limiter, muscl, kappa, toler

  namelist /flux/ flux_type, k2, k4

  namelist /conditions/ restart, mref, to, po, pback

  namelist /gas_properties/ gamma, r

contains

!================================= read_nml ==================================80
!
! Reads the q1d.nml file
!
!=============================================================================80
  subroutine read_nml

    use set_precision, only : dp

    integer :: nml_unit

    continue

    nml_unit = find_available_unit()

    open(nml_unit, file='q1d.nml', status='old')

! set defaults and read &code_control

    iterations = 100000
    !sets the maximum number of iterations

    firstorder = 10000
    !if running 2nd order, run this many 1st order iters to set up flow field

    itercheck  = 1000
    !check for convergence every itercheck iteration

    iter_out   = -1
    !output a solution every iter_out iteration, -1 if only at the end of run

    rkorder    = 1
    !sets the m-step Runge-Kutta order, 1 for Euler explicit

    muscl      = .false.
    !turns MUSCL scheme on or off

    cfl        = 1.0_dp
    !sets the CFL number

    kappa      = -1.0_dp
    !sets the type of upwinding, -1 for fully upwind

    toler      = 1.0e-13_dp
    !sets the resisdual convergence tolerance

    limiter    = 'none'
    !sets the type of variable limiter, only needed with MUSCLE = .true.

    rewind(nml_unit)
    read(nml_unit, nml=code_control)

! set defaults and read &flux

    flux_type  = 'jst'
    !sets the flux type

    k2 = 0.5_dp
    !sets the adaptive dissipation for central difference flux

    k4 = 0.03125_dp
    !sets the adaptive dissipation for central difference flux

    rewind(nml_unit)
    read(nml_unit, nml=flux)

! set defaults and read &conditions

    restart = .false.
    !start code from restart file

    mref  = 1.5_dp
    !reference inflow mach number if starting from scratch

    to    = 600.0_dp
    !reference total temperature, Kelvins

    po    = 300000.0_dp
    !reference total pressure, Pa

    pback = -1.0_dp
    !back pressure, Pa.  set less than 0 for supersonic extrapolation

    rewind(nml_unit)
    read(nml_unit, nml=conditions)

! set defaults and read &gas_properties

    r = 287.0_dp
    !gas constant,

    gamma = 1.4_dp
    !ratio of specific heats, 7/5 for diatomic gas, 5/3 for monatomic

    rewind(nml_unit)
    read(nml_unit, nml=gas_properties)

    close(nml_unit)

  end subroutine read_nml

  include 'find_available_unit.f90'

end module namelist
