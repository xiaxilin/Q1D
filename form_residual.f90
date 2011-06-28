!=============================================================================80
!
! 
!
!=============================================================================80

subroutine form_residual(iteration, dt, prim_cc,)

  use set_precision, only : dp

! FIXME... figure out where these will come from
  use XXXXXXX, only : flux_type, firstorder, kappa, muscl, limiter

  implicit none




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

  
