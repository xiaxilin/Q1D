module namelist

  use residual,        only : muscl, kappa, limiter, firstorder,               &
                              flux_type, k2, k4, inflow_gc, inflow_face,       &
                              outflow_gc, outflow_face
  use lhs,             only : lhs_order
  use solvers,         only : iterations, itercheck, iter_out, solver,         &
                              rkorder, cfl, toler, cfl_end, cfl_ramp
  use initialize_soln, only : restart, mref, to, po, pback
  use fluid_constants, only : gamma, r

  implicit none

  private

  public :: read_nml

! Define namelist inputs

  namelist /code_control/ solver, lhs_order, iterations, firstorder, itercheck,&
                          iter_out, rkorder, cfl, limiter, muscl, kappa, toler,&
                          cfl_end, cfl_ramp

  namelist /bc/ inflow_gc, inflow_face, outflow_gc, outflow_face

  namelist /flux/ flux_type, k2, k4

  namelist /conditions/ restart, mref, to, po, pback

  namelist /gas_properties/ gamma, r

!  namelist /ref_properties/ gamma, r, l_ref, rho_ref, a_ref

contains

!================================== read_nml =================================80
!
! Reads the q1d.nml file
! FIXME: rethink where the variables live so the input file can be cleaner
!
!=============================================================================80
  subroutine read_nml

    use set_precision, only : dp

    integer :: nml_unit, err, err_tot

    continue

    nml_unit = find_available_unit()

    open(nml_unit, file='q1d.nml', status='old')

! set defaults and read &code_control

    solver = 'explicit'
    !sets solver type between explicit or implicit

    lhs_order = 1
    !sets first (1) or second (2) order lhs

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

    cfl_end    = 1.0_dp
    !sets the max CFL to ramp to

    cfl_ramp   = 1
    !number of iterations over which to ramp the CFL

    kappa      = -1.0_dp
    !sets the type of upwinding, -1 for fully upwind

    toler      = 1.0e-13_dp
    !sets the resisdual convergence tolerance

    limiter    = 'none'
    !sets the type of variable limiter, only needed with MUSCLE = .true.

    rewind(nml_unit)
    read(nml_unit, nml=code_control, iostat=err)

    if ( err < 0 ) then
      write(*,*) "NML CODE_CONTROL not found, using defaults."
    else if ( err > 0 ) then
      write(*,*) "WARNING: Error in CODE_CONTROL inputs!"
      err_tot = err_tot + err
    end if

! set defaults and read&bc

    inflow_gc    = 0

    inflow_face  = 0

    outflow_gc   = 0

    outflow_face = 0

    rewind(nml_unit)
    read(nml_unit, nml=bc, iostat=err)

    if ( err < 0 ) then
      write(*,*) "NML BC not found, using defaults."
    else if ( err > 0 ) then
      write(*,*) "WARNING: Error in BC inputs!"
      err_tot = err_tot + err
    end if

! set defaults and read &flux

    flux_type  = 'jst'
    !sets the flux type

    k2 = 0.5_dp
    !sets the adaptive dissipation for central difference flux

    k4 = 0.03125_dp
    !sets the adaptive dissipation for central difference flux

    rewind(nml_unit)
    read(nml_unit, nml=flux, iostat=err)

    if ( err < 0 ) then
      write(*,*) "NML FLUX not found, using defaults."
    else if ( err > 0 ) then
      write(*,*) "WARNING: Error in FLUX inputs!"
      err_tot = err_tot + err
    end if

! set defaults and read &conditions

    restart = .false.
    !start code from restart file

    mref  = 1.5_dp
    !reference inflow mach number if starting from scratch

    to    = 600.0_dp
    !reference total temperature, assumed SI, Kelvins

    po    = 300000.0_dp
    !reference total pressure, assumed SI, Pa

!    so    = ?

    pback = -1.0_dp
    !back pressure, Pa.  set less than 0 for supersonic extrapolation

    rewind(nml_unit)
    read(nml_unit, nml=conditions, iostat=err)

    if ( err < 0 ) then
      write(*,*) "NML CONDITIONS not found, using defaults."
    else if ( err > 0 ) then
      write(*,*) "WARNING: Error in CONDITIONS inputs!"
      err_tot = err_tot + err
    end if

! set defaults and read &gas_properties

    r = 287.058_dp
    !gas constant,

    gamma = 1.4_dp
    !ratio of specific heats, 7/5 for diatomic gas, 5/3 for monatomic

!    l_ref = 1.0_dp
    !reference length

!    rho_ref = 1.225_dp
    !reference density, assumed SI, STP

!    t_ref = 288.15_dp
    !reference temperature, assumed SI, STP

!    a_ref = sqrt(gamma*r*t_ref)
    !reference speed of sound, assumed SI, STP

! Derived reference quantities
!    p_ref = rho_ref*a_ref**2
    !reference pressure... p_ref = rho*a**2 = gamma*p, at STP p_ref = 141855 Pa

    rewind(nml_unit)
    read(nml_unit, nml=gas_properties, iostat=err)

    if ( err < 0 ) then
      write(*,*) "NML GAS_PROPERTIES not found, using defaults."
    else if ( err > 0 ) then
      write(*,*) "WARNING: Error in GAS_PROPERTIES inputs!"
      err_tot = err_tot + err
    end if

    close(nml_unit)

  end subroutine read_nml

  include 'find_available_unit.f90'

end module namelist
