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
    firstorder = 10000
    itercheck  = 1000
    iter_out   = -1
    rkorder    = 1
    muscl      = .false.
    cfl        = 1.0_dp
    kappa      = -1.0_dp
    toler      = 1.0e-13_dp
    limiter    = 'none'

    rewind(nml_unit)
    read(nml_unit, nml=code_control)

! set defaults and read &flux

    flux_type  = 'jst'
    k2 = 0.5_dp
    k4 = 0.03125_dp

    rewind(nml_unit)
    read(nml_unit, nml=flux)

! set defaults and read &conditions

    restart = .false.
    mref  = 1.5_dp
    to    = 600.0_dp
    po    = 300000.0_dp
    pback = -1.0_dp

    rewind(nml_unit)
    read(nml_unit, nml=conditions)

! set defaults and read &gas_properties

    r = 287.0_dp
    gamma = 1.4_dp

    rewind(nml_unit)
    read(nml_unit, nml=gas_properties)

    close(nml_unit)

  end subroutine read_nml

  include 'find_available_unit.f90'

end module namelist
