module jacobians

  implicit none

  private

  public :: jac_vanleer_1D

contains

!========================= van_leer_jacobian =================================80
!
! This subroutine takes the left and right *conserved* variables at a face
! and returns the left and right flux jacobians
!
!=============================================================================80
  subroutine jac_vanleer_1D( qL, qR, jac_l, jac_r )

    use set_precision,   only : dp
    use set_constants,   only : zero, fourth, half, one, two, four
    use fluid_constants, only : gamma, gm1, xg, xg2m1

    implicit none

    real(dp), dimension(3),   intent(in)  :: ql, qr
    real(dp), dimension(3,3), intent(out) :: jac_l, jac_r

    real(dp) :: rho, rhoinv, u, p, a, m, fa, fb

    real(dp), dimension(3)   :: drho_dq, du_dq, dp_dq, dq3_dq
    real(dp), dimension(3)   :: da_dq, dm_dq, dfa_dq, dfb_dq
    real(dp), dimension(3,3) :: mat

    continue

! handle left side first

! calculate primitive vars from conserved
    rho = qL(1)
    rhoinv = one/rho
    u = qL(2)*rhoinv
    p = gm1*( qL(3) - half * rho * u**2 )
    a = speed_of_sound(p, rho)
    m = u/a

! linearization of primitive variables
    drho_dq(1) = one
    drho_dq(2) = zero
    drho_dq(3) = zero

    du_dq(1) = -u*rhoinv
    du_dq(2) = rhoinv
    du_dq(3) = zero
        
    dp_dq(1) = half*gm1*u**2
    dp_dq(2) = -gm1*u
    dp_dq(3) = gm1

    dq3_dq(1) = zero
    dq3_dq(2) = zero
    dq3_dq(3) = one

! speed of sound linearization
    da_dq(:) = (half*gamma/a) * ( rhoinv*dp_dq(:) - (p*rhoinv**2)*drho_dq(:) )

! mach number linearization

    dm_dq(:) = du_dq(:)/a - (u/a**2)*da_dq(:)

! linearization of FVS mass term
    fa =  fourth*rho*a*(m+one)**2

    dfa_dq(:) =  fourth*( drho_dq(:)*a*(m+one)**2                              &
              + rho*da_dq(:)*(m+one)**2 + rho*a*two*(m+one)*dm_dq(:) ) 

! linearization of FVS energy term
    fb = a*(gm1*m + two)

    dfb_dq(:) = da_dq(:)*(gm1*m + two) + dm_dq(:)*gm1*a

!    Flux = (/ fa, fa*fb*xg, half*fa*fb*fb*xg2m1 /)

! now, formulate jacobian

    mat = zero

    if( abs(m)<1.0 ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))

    elseif( m>=1.0 ) then
      mat(1,:) = drho_dq(:)*u + rho*du_dq(:)
      mat(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
      mat(3,:) = ( dq3_dq(:) + dp_dq(:) )*u + (ql(3) + p)*du_dq(:)

    endif

    jac_l = mat

! handle right side second

! calculate primitive vars from conserved
    rho = qR(1)
    rhoinv = one/rho
    u = qR(2)*rhoinv
    p = gm1*( qR(3) - half * rho * u**2 )
    a = speed_of_sound(p,rho)
    m = u/a

! linearization of right primitive variables
    drho_dq(1) = one
    drho_dq(2) = zero
    drho_dq(3) = zero
       
    du_dq(1) = -u*rhoinv
    du_dq(2) = rhoinv
    du_dq(2) = zero
        
    dp_dq(1) =  half*gm1*u**2
    dp_dq(2) = -gm1*u
    dp_dq(3) =  gm1

    dq3_dq(1) = zero
    dq3_dq(2) = zero
    dq3_dq(3) = one

! speed of sound linearization
    da_dq(:) = (half*gamma/a) * ( rhoinv*dp_dq(:) - (p*rhoinv**2)*drho_dq(:) )

! mach number linearization
    dm_dq(:) = du_dq(:)/a - (u/a**2)*da_dq(:)

! linearization of FVS mass term
    fa = -fourth*rho*a*(m-one)**2

    dfa_dq(:) = -fourth*( drho_dq(:)*a*(m-one)**2                              &
              + rho*da_dq(:)*(m-one)**2 + rho*a*two*(m-one)*dm_dq(:) )

! linearization of FVS energy term

    fb = a*(gm1*m - two)

    dfb_dq(:) = da_dq(:)*(gm1*m - two) + dm_dq(:)*gm1*a

! now, form jacobian

    mat = zero

    if( abs(m)<1.0 ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))

    elseif( m<=-1.0 ) then
      mat(1,:) = drho_dq(:)*u   + rho*du_dq(:)
      mat(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
      mat(3,:) = ( dq3_dq(:) + dp_dq(:) )*u + (qR(3) + p)*du_dq(:)

    endif

    jac_r = mat

  end subroutine jac_vanleer_1D

!  subroutine jac_source_1D
!  end subroutine jac_source_1D

  include 'speed_of_sound.f90'

end module jacobians
