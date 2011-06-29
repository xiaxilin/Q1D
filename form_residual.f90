module form_residual

  implicit none

  private

  public :: create_residual

contains

!=============================================================================80
!
! 
!
!=============================================================================80

  subroutine create_residual( cells, faces, iteration, prim_cc,                &
                              area_f, area_cc, dadx_cc, residual )

    use set_precision, only : dp
    use set_constants, only : zero

    implicit none

    integer,                        intent(in)  :: cells, faces, iterations
    real(dp), dimension(3,cells+2), intent(in)  :: prim_cc
    real(dp), dimension(faces),     intent(in)  :: area_f
    real(dp), dimension(cells+2),   intent(in)  :: area_cc
    real(dp), dimension(cells+2),   intent(in)  :: dadx_cc
    real(dp), dimension(3,cells+2), intent(out) :: residual

    integer :: i
    real(dp), dimension(3,faces) :: F, S

    continue

    call create_fluxes( faces, iteration, prim_cc, F )
    call create_source( cells, prim_cc(3,:), dadx_cc, S )

    residual = zero
    do i = 2, cells+1
      residual(:,i) = ( dx*S(:,i) - area_f(i)*F(:,i) + area_f(i-1)*F(:i-1) )   &
                    / (dx*area_cc(i))
    end do

  end subroutine create_residual

!=============================================================================80
!
! 
!
!=============================================================================80

  subroutine create_fluxes( faces, iteration, prim_cc, flux )

    use set_precision, only : dp
! FIXME... figure out where these will come from
    use XXXXXXX, only : flux_type

    implicit none

    integer,                     intent(in)  :: faces, iteration
    real, dimension(3, faces+1), intent(in)  :: prim_cc
    real, dimension(3, faces),   intent(out) :: flux

    real, allocatable, dimension(:,:) :: prim_left, prim_right
    integer :: i

    continue

! call muscl extrapolation only where appropriate
    if (trim(flux_type) /= 'jst' .or. trim(flux_type) /= 'central') then
      allocate( prim_left(3, faces), prim_right(3, faces) )
      call muscl_extrapolation(iteration, prim_cc, prim_left, prim_right)
      do i = 1, faces
        prim_left(:,i)  = floor_primitive_vars(prim_left(:,i))
        prim_right(:,i) = floor_primitive_vars(prim_right(:,i))
      end do
    end if

    select case(trim(flux_type))
    case ('central')
      do i = 1, faces
        flux(:,i) = flux_central( prim_cc(:,i), prim_cc(:,i+1) )
      end do
    case ('jst')
      do i = 1, faces
        flux(:,i) = flux_central( prim_cc(:,i), prim_cc(:,i+1) )
      end do
      call add_jst_damping( prim_cc, flux )
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
    real(dp), dimension(1, cells+2), intent(in)  :: pressure, dadx_cc
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
! FIXME: broken, doesn't have conserved variables
!
!=============================================================================80

  subroutine add_jst_damping( cells, faces, prim_cc, central_flux )

    use set_precision, only : dp
    use set_constants, only : zero, half, two, three

    implicit none

    integer,                         intent(in)    :: cells, faces
    real(dp), dimension(3, cells+2), intent(in)    :: prim_cc
    real(dp), dimension(3, faces),   intent(inout) :: central_flux

    integer                      :: i
    real(dp), dimension(cells+2) :: nu
    real(dp), dimension(3)       :: dissipation, cons_left, cons_right
    real(dp)                     :: lambda, epstwo, epsfour

    continue

! calculate pressure switch
    do i = 2, cells+1
      nu(i) = abs((prim_cc(3,i-1) - two*prim_cc(3,i) + prim_cc(3,i+1))         &
            /     (prim_cc(3,i-1) + two*prim_cc(3,i) + prim_cc(3,i+1)))
    end do
    nu(1)       = two*nu(2) - nu(3)
    nu(cells+2) = zero

! calculate smoothing terms and dissipation vector, add to flux
    i = 1
    epstwo  = ktwo*max(nu(i), nu(i+1), nu(i+2))
    epsfour = max(zero, kfour-epstwo)

    lambda  = half*(abs(prim_cc(2,i+1))+a(i+1) + abs(prim_cc(2,i))+a(i))

    dissipation(:) = -lambda* (epstwo*(Q(:,i+1) - Q(:,i))                     &
                   - epsfour*(Q(:,i+2) - three*Q(:,i+1) + two*Q(:,i)))

    central_flux(:,i) = central_flux(:,i) + dissipation(:)
                
    do i = 2, faces-1
      epstwo  = ktwo*max(nu(i-1), nu(i), nu(i+1), nu(i+2))
      epsfour = max(zero, kfour-epstwo)

      lambda  = half*(abs(prim_cc(2,i+1))+a(i+1) + abs(prim_cc(2,i))+a(i))

      dissipation(:) = -lambda*(epstwo*(Q(:,i+1)-Q(:,i))                      &
        - epsfour*(Q(:,i+2) - three*Q(:,i+1) + three*Q(:,i) - Q(:,i-1)))

      central_flux(:,i) = central_flux(:,i) + dissipation(:)
    end do

    i = faces
    epstwo  = ktwo*max(nu(i-1), nu(i), nu(i+1))
    epsfour = max(zero, kfour-epstwo)

    lambda  = half*(abs(prim_cc(2,i+1))+a(i+1) + abs(prim_cc(2,i))+a(i))

    dissipation(:) = -lambda*(epstwo*(Q(:,i+1) - Q(:,i))                      &
                   - epsfour*(-two*Q(:,i+1) + three*Q(:,i) - Q(:,i-1)))

    central_flux(:,i) = central_flux(:,i) + dissipation(:)

  end subroutine add_jst_damping

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine muscl_extrapolation( cells, faces, iterations, vars_cc,           &
                                  vars_left, vars_right )

    use set_precision, only : dp
    use set_constants, only : zero, fourth, one, small
! FIXME: figure out where these will come from
    use XXXXXX,        only :  firstorder, muscl, limiter, kappa

    implicit none

    integer,                        intent(in)  :: cells, faces, iterations
    real(dp), dimension(3,cells+2), intent(in)  :: vars_cc
    real(dp), dimension(3,faces),   intent(out) :: vars_left, vars_right

    integer                       :: i
    real(dp), dimension(3, faces) :: r_L, r_R, psi_L, psi_R

    real(dp), parameter :: small_factor = 0.0000001_dp

    continue

    if ( iterations <= firstorder .or. .not. muscl ) then
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
      r_L(:,faces) = one

! calculate right side variations, r>=0
      r_R(:,1)     = one
      do i = 2, faces
        r_R(:,i) = max( zero, (vars_cc(:,i) - vars_cc(:,i-1))                  &
                            / (vars_cc(:,i+1) - vars_cc(:,i) + small_factor) )
      end do

! apply appropriate limiter

      select case(trim(limiter))
      case('minmod')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(r_L(i), 1)
          psi_R(:,i) = limiter_sweby(r_R(i), 1)        
        end do
      case('superbee')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(r_L(i), 2)
          psi_R(:,i) = limiter_sweby(r_R(i), 2)  
        end do
      case('sweby')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(r_L(i), 1.5)
          psi_R(:,i) = limiter_sweby(r_R(i), 1.5)
        end do
      case('ospre')
        do i = 1, faces
          psi_L(:,i) = limiter_ospre(r_L(i))
          psi_R(:,i) = limiter_ospre(r_R(i))
        end do
      case('vanleer')
        do i = 1, faces
          psi_L(:,i) = limiter_vanleer(r_L(i))
          psi_R(:,i) = limiter_vanleer(r_R(i))
        end do
      case('vanalbada')
        do i = 1, faces
          psi_L(:,i) = limiter_vanalbada(r_L(i))
          psi_R(:,i) = limiter_vanalbada(r_R(i))
        end do
      case( default )
        psi_L(:,:) = one
        psi_R(:,:) = one
      end select

! perform MUSCL extrapolation
      i = 1
      vars_left(:,i) = vars_cc(:,i) + fourth                                   &
                * ((one+kappa)*(vars_cc(:,i+1) - vars_cc(:,i)))

      do i = 2, faces

        vars_left(:,i)    = vars_cc(:,i) + fourth                              &
                * ((one+kappa)*psi_R(:,i)   * (vars_cc(:,i+1) - vars_cc(:,i))  &
                +  (one-kappa)*psi_L(:,i-1) * (vars_cc(:,i)   - vars_cc(:,i-1)))

        vars_right(:,i-1) = vars_cc(:,i) - fourth                              &
                * ((one-kappa)*psi_R(:,i)   * (vars_cc(:,i+1) - vars_cc(:,i))  &
                +  (one+kappa)*psi_L(:,i-1) * (vars_cc(:,i)   - vars_cc(:,i-1)))

      end do

      i = faces
      vars_right(:,i) = vars_cc(:,i) - fourth                                  &
                * ((one+kappa)*(vars_cc(:,i) - vars_cc(:,i-1)))

    end if

  end subroutine muscl_extrapolation

! begin include statements
  include 'floor_primitive_vars.f90'
  include 'flux_central.f90'
  include 'flux_vanleer.f90'
  include 'flux_sw.f90'
  include 'flux_ausm.f90'
  include 'flux_roe.f90'
!  include 'flux_*.f90'

  include 'limiter_sweby.f90'
  include 'limiter_ospre.f90'
  include 'limiter_vanleer.f90'
  include 'limiter_vanalbada.f90'
!  include 'limiter_*.f90'

end module form_residual
