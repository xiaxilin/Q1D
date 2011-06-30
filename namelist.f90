module namelist

  use set_precision,   only : dp
  use solvers,         only : iterations, firstorder, itercheck, rkorder, cfl, &
                              limiter, muscl, kappa, toler, flux_type, k2, k4
  use initialize_soln, only : mref, to, po, pback
  use fluid_constants, only : gamma, r

  implicit none

  private

  public :: read_nml

! Define namelist inputs

  integer           :: iterations, firstorder, itercheck, rkorder
  real(dp)          :: cfl, kappa, toler
  character(len=10) :: limiter
  logical           :: muscl
  namelist /code_control/ iterations, firstorder, itercheck, rkorder, cfl, &
                          limiter, muscl, kappa, toler

  character(len=10) :: flux_type
  real(dp)          :: k2, k4
  namelist /flux/ flux_type, k2, k4

  real(dp) :: mref, to, po, pback
  namelist /conditions/ mref, to, po, pback

  real(dp) :: gamma, r
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

    rewind(nml_unit)
    read(nml_unit, nml=code_control)

! set defaults and read &flux

    rewind(nml_unit)
    read(nml_unit, nml=flux)

! set defaults and read &conditions

    rewind(nml_unit)
    read(nml_unit, nml=conditions)

! set defaults and read &gas_properties

    rewind(nml_unit)
    read(nml_unit, nml=gas_properties)

    close(nml_unit)

  end subroutine read_nml

  include 'find_available_unit.f90'

end module namelist
