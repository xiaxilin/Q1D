module namelist

  use set_precision, only : dp

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
  real(dp)          :: k2, k4    ! only for flux_type = 'jst'
  namelist /flux/ flux_type, k2, k4

  real(dp) :: mref, to, po pback
  namelist /conditions/ mreft, to, po, pback

  real(dp) :: gamma, r
  namelist /gas_properties/ gamma, r

contains

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
    rkorder    = 1
    cfl        = 1.0
    limiter    = 'none'
    muscl      = .false.
    kappa      = -1.0
    toler      = 1.0e-13

    rewind(nml_unit)
    read(nml_unit, nml=code_control)

! set defaults and read &flux

    flux_type  = 'jst'
    k2         = 0.5
    k4         = 0.03125

    rewind(nml_unit)
    read(nml_unit, nml=flux)

! set defaults and read &conditions

    mref  = 1.5
    to    = 600.0
    po    = 300000.0
    pback = -1.0

    rewind(nml_unit)
    read(nml_unit, nml=conditions)

! set defaults and read &gas_properties

    gamma = 1.4
    r     = 287.0

    rewind(nml_unit)
    read(nml_unit, nml=gas_properties)

    close(nml_unit)

  end subroutine read_nml

end module namelist
