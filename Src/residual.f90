module residual

  use set_precision, only : dp

  implicit none

  private

  public :: create_residual
  public :: muscl_extrapolation

  public :: muscl      ! MUSCL extrapolation, .true. or .false.
  public :: kappa      ! form of MUSCL
  public :: limiter    ! type of variable limiting
  public :: firstorder ! if 2nd order, how many 1st order iters for stability?

  public :: flux_type  ! 'central', 'jst', 'sw', 'vanleer', 'roe', 'ausm'
  public :: k2         ! JST damping coefficient, only for flux_type = 'jst'
  public :: k4         ! JST damping coefficient, only for flux_type = 'jst'

  logical       :: muscl
  real(dp)      :: kappa
  character(10) :: limiter
  integer       :: firstorder

  character(10) :: flux_type
  real(dp)      :: k2
  real(dp)      :: k4

contains

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine create_residual( cells, faces, iteration, prim_cc, cons_cc,       &
                              area_f, dadx_cc, dx, residual )

    use set_precision, only : dp
    use set_constants, only : zero

    implicit none

    integer,                        intent(in)  :: cells, faces, iteration
    real(dp), dimension(3,cells+2), intent(in)  :: prim_cc
    real(dp), dimension(3,cells+2), intent(in)  :: cons_cc
    real(dp), dimension(faces),     intent(in)  :: area_f
    real(dp), dimension(cells+2),   intent(in)  :: dadx_cc
    real(dp), dimension(cells+2),   intent(in)  :: dx
    real(dp), dimension(3,cells+2), intent(out) :: residual

    integer :: cell
    real(dp), dimension(3,faces)   :: F
    real(dp), dimension(3,cells+2) :: S

    continue

    residual = zero

    call create_fluxes( cells, faces, iteration, prim_cc, cons_cc, F )
    call create_source( cells, prim_cc(3,:), dadx_cc, S )

    do cell = 2, cells+1
      residual(:,cell) = area_f(cell)*F(:,cell) - area_f(cell-1)*F(:,cell-1)   &
                       - S(:,cell)*dx(cell)
    end do

  end subroutine create_residual

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine create_fluxes( cells, faces, iteration, prim_cc, cons_cc, flux )

    use set_precision, only : dp

    implicit none

    integer,                         intent(in)  :: cells, faces, iteration
    real(dp), dimension(3, cells+2), intent(in)  :: prim_cc
    real(dp), dimension(3, cells+2), intent(in)  :: cons_cc
    real(dp), dimension(3, faces),   intent(out) :: flux

    real(dp), dimension(3, faces) :: prim_left, prim_right
    integer :: i

    continue

! call muscl extrapolation only where appropriate
    if (trim(flux_type) /= 'jst' .or. trim(flux_type) /= 'central') then
      call muscl_extrapolation(cells, faces, iteration, &
                               prim_cc, prim_left, prim_right)

      do i = 1, faces
        prim_left(:,i)  = floor_primitive_vars(prim_left(:,i))
        prim_right(:,i) = floor_primitive_vars(prim_right(:,i))
      end do

    end if

! create flux vectors
    select case(trim(flux_type))
    case ('central')
      do i = 1, faces
        flux(:,i) = flux_central( prim_cc(:,i), prim_cc(:,i+1) )
      end do
    case ('jst')
      do i = 1, faces
        flux(:,i) = flux_central( prim_cc(:,i), prim_cc(:,i+1) )
      end do
      call add_jst_damping( cells, faces, prim_cc, cons_cc, flux )
    case('vanleer')
      do i = 1, faces
        flux(:,i) = flux_vanleer( prim_left(:,i), prim_right(:,i) )
      end do
    case('sw')
      do i = 1, faces
        flux(:,i) = flux_sw( prim_left(:,i), prim_right(:,i) )
      end do
    case('ausm')
      do i = 1, faces
        flux(:,i) = flux_ausm( prim_left(:,i), prim_right(:,i) )
      end do
    case('ausmplus')
      do i = 1, faces
        flux(:,i) = flux_ausm_plus( prim_left(:,i), prim_right(:,i) )
      end do
    case('roe')
      do i = 1, faces
        flux(:,i) = flux_roe( prim_left(:,i), prim_right(:,i) )
      end do
    case('hllc')

    end select

  end subroutine create_fluxes

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine create_source( cells, pressure, dadx_cc, source )

    use set_precision, only : dp
    use set_constants, only : zero

    implicit none

    integer,                         intent(in)  :: cells
    real(dp), dimension(cells+2),    intent(in)  :: pressure, dadx_cc
    real(dp), dimension(3, cells+2), intent(out) :: source

    integer :: i

    continue

    source = zero

    do i = 2, cells+1
      source(2,i) = pressure(i)*dadx_cc(i)
    end do

  end subroutine create_source

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine add_jst_damping( cells, faces, prim_cc, cons_cc, central_flux )

    use set_precision, only : dp
    use set_constants, only : zero, half, two, three

    implicit none

    integer,                         intent(in)    :: cells, faces
    real(dp), dimension(3, cells+2), intent(in)    :: prim_cc
    real(dp), dimension(3, cells+2), intent(in)    :: cons_cc
    real(dp), dimension(3, faces),   intent(inout) :: central_flux

    integer                      :: i
    real(dp), dimension(cells+2) :: nu, a
    real(dp), dimension(3)       :: dissipation
    real(dp)                     :: lambda, epstwo, epsfour

    continue

! calculate pressure switch
    do i = 2, cells+1
      nu(i) = abs((prim_cc(3,i-1) - two*prim_cc(3,i) + prim_cc(3,i+1))         &
            /     (prim_cc(3,i-1) + two*prim_cc(3,i) + prim_cc(3,i+1)))
    end do
    nu(1)       = two*nu(2) - nu(3)
    nu(cells+2) = zero

    do i = 1, cells+2
      a(i) = speed_of_sound(prim_cc(3,i), prim_cc(1,i))
    end do

! calculate smoothing terms and dissipation vector, add to flux
    i = 1
    epstwo  = k2*max(nu(i), nu(i+1), nu(i+2))
    epsfour = max(zero, k4-epstwo)

    lambda  = half*(abs(prim_cc(2,i+1))+a(i+1) + abs(prim_cc(2,i))+a(i))

    dissipation(:) = -lambda* (epstwo*(cons_cc(:,i+1) - cons_cc(:,i))          &
       - epsfour*(cons_cc(:,i+2) - three*cons_cc(:,i+1) + two*cons_cc(:,i)))

    central_flux(:,i) = central_flux(:,i) + dissipation(:)

    do i = 2, faces-1
      epstwo  = k2*max(nu(i-1), nu(i), nu(i+1), nu(i+2))
      epsfour = max(zero, k4-epstwo)

      lambda  = half*(abs(prim_cc(2,i+1))+a(i+1) + abs(prim_cc(2,i))+a(i))

      dissipation(:) = -lambda*(epstwo*(cons_cc(:,i+1)-cons_cc(:,i))           &
        - epsfour*(cons_cc(:,i+2) - three*cons_cc(:,i+1) + three*cons_cc(:,i)  &
        - cons_cc(:,i-1)))

      central_flux(:,i) = central_flux(:,i) + dissipation(:)
    end do

    i = faces
    epstwo  = k2*max(nu(i-1), nu(i), nu(i+1))
    epsfour = max(zero, k4-epstwo)

    lambda  = half*(abs(prim_cc(2,i+1))+a(i+1) + abs(prim_cc(2,i))+a(i))

    dissipation(:) = -lambda*(epstwo*(cons_cc(:,i+1) - cons_cc(:,i))           &
      - epsfour*(-two*cons_cc(:,i+1) + three*cons_cc(:,i) - cons_cc(:,i-1)))

    central_flux(:,i) = central_flux(:,i) + dissipation(:)

  end subroutine add_jst_damping

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine muscl_extrapolation( cells, faces, iteration, vars_cc,           &
                                  vars_left, vars_right, limL, limR )

    use set_precision, only : dp
    use set_constants, only : zero, fourth, half, one, onep5, two, small

    implicit none

    integer,                        intent(in)  :: cells, faces, iteration
    real(dp), dimension(3,cells+2), intent(in)  :: vars_cc
    real(dp), dimension(3,faces),   intent(out) :: vars_left, vars_right
    real(dp), optional, dimension(3,cells+2), intent(out) :: limL, limR

    integer                       :: i
    real(dp), dimension(3, faces) :: r_L, r_R, psi_L, psi_R

    real(dp), parameter :: small_factor = 0.0000001_dp

    continue

    second_order : if ( iteration <= firstorder .or. .not. muscl ) then
! skip excess computations if only first order
      do i = 1, faces
        vars_left(:,i)  = vars_cc(:,i)
        vars_right(:,i) = vars_cc(:,i+1)
      end do
    else

! calculate left side variations, r>=0
      do i = 1, faces-1
        r_L(:,i)  = max( zero, (vars_cc(:,i+2) - vars_cc (:,i+1))              &
                             / (vars_cc(:,i+1) - vars_cc(:,i) + small_factor) )
      end do
      r_L(:,faces) = zero

! calculate right side variations, r>=0

      r_R(:,1)     = zero
      do i = 2, faces
        r_R(:,i) = max( zero, (vars_cc(:,i) - vars_cc(:,i-1))                  &
                            / (vars_cc(:,i+1) - vars_cc(:,i) + small_factor) )
      end do

! apply appropriate limiter
      select case(trim(limiter))
      case('minmod')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(3, r_L(:,i), one)
          psi_R(:,i) = limiter_sweby(3, r_R(:,i), one)
        end do
      case('superbee')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(3, r_L(:,i), two)
          psi_R(:,i) = limiter_sweby(3, r_R(:,i), two)
        end do
      case('sweby')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(3, r_L(:,i), onep5)
          psi_R(:,i) = limiter_sweby(3, r_R(:,i), onep5)
        end do
      case('ospre')
        do i = 1, faces
          psi_L(:,i) = limiter_ospre(3, r_L(:,i))
          psi_R(:,i) = limiter_ospre(3, r_R(:,i))
        end do
      case('vanleer')
        do i = 1, faces
          psi_L(:,i) = limiter_vanleer(3, r_L(:,i))
          psi_R(:,i) = limiter_vanleer(3, r_R(:,i))
        end do
      case('vanalbada')
        do i = 1, faces
          psi_L(:,i) = limiter_vanalbada(3, r_L(:,i))
          psi_R(:,i) = limiter_vanalbada(3, r_R(:,i))
        end do
      case default
        psi_L(:,:) = one
        psi_R(:,:) = one
      end select

! perform MUSCL extrapolation... note that there is no limiting at in/outflow
! the commented out portions of code are old versions
      i = 1
      vars_left(:,i) = vars_cc(:,i)
!     vars_left(:,i) = vars_cc(:,i) + fourth                                   &
!               * ((one+kappa)*(vars_cc(:,i+1) - vars_cc(:,i)))

      do i = 2, faces
        vars_left(:,i)  = vars_cc(:,i) + fourth                                &
                * ((one+kappa)*psi_R(:,i)   * (vars_cc(:,i+1) - vars_cc(:,i))  &
                +  (one-kappa)*psi_L(:,i-1) * (vars_cc(:,i)   - vars_cc(:,i-1)))
      end do

      do i = 1, faces-1
        vars_right(:,i) = vars_cc(:,i+1) - fourth                              &
                * ((one-kappa)*psi_R(:,i+1) * (vars_cc(:,i+2) - vars_cc(:,i+1))&
                +  (one+kappa)*psi_L(:,i)   * (vars_cc(:,i+1) - vars_cc(:,i)))
      end do

      i = faces
      vars_right(:,i) = vars_cc(:,i+1)
!     vars_right(:,i) = vars_cc(:,i+1) - fourth                                &
!               * ((one+kappa)*(vars_cc(:,i+1) - vars_cc(:,i)))

! If extrapolated density or pressure/energy is < 0 then
! cancel the extrapolation and go first order

      do i = 1, faces
        if (vars_left(1,i) <= zero .or. vars_left(3,i) <= zero) then
          vars_left(:,i) = vars_cc(:,i)
        end if
        if (vars_right(1,i) <= zero .or. vars_right(3,i) <= zero) then
          vars_right(:,i) = vars_cc(:,i+1)
        end if
      end do

      if ( present(limL) .and. present(limR) ) then
        limL = psi_L
        limR = psi_R
      end if

    end if second_order

  end subroutine muscl_extrapolation

  include 'conserved_to_primitive_1D.f90'
  include 'floor_primitive_vars.f90'
  include 'primitive_to_conserved_1D.f90'
  include 'speed_of_sound.f90'

  include 'flux_central.f90'
  include 'flux_vanleer.f90'
  include 'flux_sw.f90'
  include 'flux_ausm.f90'
  include 'flux_ausm_plus.f90'
  include 'flux_roe.f90'
!  include 'flux_*.f90'

  include 'limiter_sweby.f90'
  include 'limiter_ospre.f90'
  include 'limiter_vanleer.f90'
  include 'limiter_vanalbada.f90'
!  include 'limiter_*.f90'

end module residual
