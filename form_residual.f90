!=============================================================================80
!
! 
!
!=============================================================================80

subroutine form_residual( cells, faces, iteration, dt, prim_cc,)

  use set_precision, only : dp
  use set_constants, only : one

! FIXME... figure out where these will come from
  use XXXXXXX, only : flux_type



! FIXME... needed by muscl_extrapolation
  use XXXXXX, only : firstorder, kappa, muscl, limiter

  implicit none

  integer, intent(in) :: cells, faces, iterations
  real(dp), 


  integer :: i

  continue

! FIXME: formulate as:
!  call create_fluxes
!  call create_source

!  do i = 2, cells+1
!    residual(:,i) = (one/(dx*area_cc(i))) 
!                  * ( dx*S(:,i) - (F(:,i)*area_f(i)) + (F(:,i-1)*area_f(i-1)) )
!  end do

end subroutine form_residual

!=============================================================================80
!
! 
!
!=============================================================================80

subroutine create_fluxes( faces, iteration, prim_cc, flux )

  use set_precision, only : dp

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
    call muscl_extrapolation(iteration, prim_left, prim_right)
    do i = 1, faces
      call floor_primitive_vars(prim_left(:,i))
      call floor_primitive_vars(prim_right(:,i))
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

! FIXME: includes!

 include 'flux_*.f90'
 include 'flux_*.f90'
 include 'flux_*.f90'
 include 'flux_*.f90'
 include 'flux_*.f90'
 include 'flux_*.f90'
