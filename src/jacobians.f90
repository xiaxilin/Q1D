module jacobians

  implicit none

  private

  public :: jac_vanleer_1D, jac_vanleer_q_1D
  public :: jac_source_1D, jac_source_q_1D
  public :: dconserved_dprimitive

contains

!=============================== jac_vanleer_1D ==============================80
!
! This subroutine takes the left and right *conserved* variables at a face
! and returns the left and right flux jacobians wrt conserved variables
! This Jacobian has been checked and found to be correct.
!
!=============================================================================80
  subroutine jac_vanleer_1D( qL, qR, jac_l, jac_r )

    use set_precision,   only : dp
    use set_constants,   only : zero, fourth, half, one, two
    use fluid_constants, only : gamma, gm1, xg, xg2m1

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

! linearization of primitive variables wrt conserved
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

!    Flux = [ fa, fa*fb*xg, half*fa*fb*fb*xg2m1 ]

! now, formulate jacobian

    mat = zero

    if( abs(m)<1.0_dp ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))

    elseif( m>=1.0_dp ) then
      mat(1,:) = drho_dq(:)*u + rho*du_dq(:)
      mat(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
      mat(3,:) = ( dq3_dq(:) + dp_dq(:) )*u + (qL(3) + p)*du_dq(:)

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

! linearization of right primitive variables wrt conserved
    drho_dq(1) = one
    drho_dq(2) = zero
    drho_dq(3) = zero

    du_dq(1) = -u*rhoinv
    du_dq(2) = rhoinv
    du_dq(3) = zero

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

    if( abs(m)<1.0_dp ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))

    elseif( m<=-1.0_dp ) then
      mat(1,:) = drho_dq(:)*u   + rho*du_dq(:)
      mat(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
      mat(3,:) = ( dq3_dq(:) + dp_dq(:) )*u + (qR(3) + p)*du_dq(:)

    endif

    jac_r = mat

  end subroutine jac_vanleer_1D

!=============================== jac_central_1D ==============================80
!
! This subroutine takes the left and right *conserved* variables at a face
! and returns the left and right flux jacobians wrt conserved variables
!
!=============================================================================80
  subroutine jac_central_1D(qL, qR, jac_l, jac_r)

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two
    use fluid_constants, only : gm1

    real(dp), dimension(3),   intent(in)  :: ql, qr
    real(dp), dimension(3,3), intent(out) :: jac_l, jac_r

    real(dp) :: rho, rhoinv, u, p

    real(dp), dimension(3) :: drho_dq, du_dq, dp_dq, dq3_dq

    continue

! Common linearizations
    drho_dq(1) = one
    drho_dq(2) = zero
    drho_dq(3) = zero

    dq3_dq(1) = zero
    dq3_dq(2) = zero
    dq3_dq(3) = one

! Do the left side...
! calculate primitive vars from conserved
    rho = qL(1)
    rhoinv = one/rho
    u = qL(2)*rhoinv
    p = gm1*( qL(3) - half * rho * u**2 )

! linearization of primitive variables wrt conserved
    du_dq(1) = -u*rhoinv
    du_dq(2) = rhoinv
    du_dq(3) = zero

    dp_dq(1) = half*gm1*u**2
    dp_dq(2) = -gm1*u
    dp_dq(3) = gm1

    jac_L(1,:) = drho_dq(:)*u + rho*du_dq(:)
    jac_L(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
    jac_L(3,:) = ( dq3_dq(:) + dp_dq(:) )*u + (qL(3) + p)*du_dq(:)

    jac_L = half*jac_L

! Do the right side...
! calculate primitive vars from conserved
    rho = qR(1)
    rhoinv = one/rho
    u = qR(2)*rhoinv
    p = gm1*( qR(3) - half * rho * u**2 )

! linearization of right primitive variables wrt conserved
    du_dq(1) = -u*rhoinv
    du_dq(2) = rhoinv
    du_dq(3) = zero

    dp_dq(1) =  half*gm1*u**2
    dp_dq(2) = -gm1*u
    dp_dq(3) =  gm1

    jac_R(1,:) = drho_dq(:)*u + rho*du_dq(:)
    jac_R(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
    jac_R(3,:) = ( dq3_dq(:) + dp_dq(:) )*u + (qR(3) + p)*du_dq(:)

    jac_R = half*jac_R

  end subroutine jac_central_1D

!================================ jac_ausm_1D ================================80
!
! This subroutine takes the left and right *conserved* variables at a face
! and returns the left and right flux jacobians wrt conserved variables
!
!=============================================================================80
  subroutine jac_ausm_1D( qL, qR, jac_l, jac_r )

    use set_precision,   only : dp
    use set_constants,   only : zero, fourth, half, one, two
    use fluid_constants, only : gamma, gm1, xg, xg2m1

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

! linearization of primitive variables wrt conserved
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

!    Flux = [ fa, fa*fb*xg, half*fa*fb*fb*xg2m1 ]

! now, formulate jacobian

    mat = zero

    if( abs(m)<=1.0_dp ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))

    elseif( m>1.0_dp ) then
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

! linearization of right primitive variables wrt conserved
    drho_dq(1) = one
    drho_dq(2) = zero
    drho_dq(3) = zero

    du_dq(1) = -u*rhoinv
    du_dq(2) = rhoinv
    du_dq(3) = zero

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

    if( abs(m)<=1.0_dp ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))

    elseif( m<-1.0_dp ) then
      mat(1,:) = drho_dq(:)*u   + rho*du_dq(:)
      mat(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
      mat(3,:) = ( dq3_dq(:) + dp_dq(:) )*u + (qR(3) + p)*du_dq(:)

    endif

    jac_r = mat

  end subroutine jac_ausm_1D

!=============================== jac_source_1D ===============================80
!
! This subroutine returns the source Jacobian wrt conserved variables
!
!=============================================================================80
  subroutine jac_source_1D(vel, dadx_cc, cell_jac, source_jac)

    use set_precision,   only : dp
    use set_constants,   only : zero, half
    use fluid_constants, only : gm1

    real(dp),                 intent(in)  :: vel, dadx_cc, cell_jac
    real(dp), dimension(3,3), intent(out) :: source_jac

    real(dp) :: sidewall_area

    continue

    sidewall_area = dadx_cc*cell_jac

    source_jac(1,1) = zero
    source_jac(2,1) = half*gm1*vel**2*sidewall_area
    source_jac(3,1) = zero
    source_jac(1,2) = zero
    source_jac(2,2) = -gm1*vel*sidewall_area
    source_jac(3,2) = zero
    source_jac(1,3) = zero
    source_jac(2,3) = gm1*sidewall_area
    source_jac(3,3) = zero

  end subroutine jac_source_1D

!============================== jac_vanleer_q_1D =============================80
!
! This subroutine takes the left and right *primitive* variables at a face
! and returns the left and right flux jacobians wrt primitive variables
!
!=============================================================================80
  subroutine jac_vanleer_q_1D( qL, qR, jac_l, jac_r )

    use set_precision,   only : dp
    use set_constants,   only : zero, fourth, half, one, two
    use fluid_constants, only : gamma, gm1, xg, xgm1, xg2m1

    real(dp), dimension(3),   intent(in)  :: ql, qr
    real(dp), dimension(3,3), intent(out) :: jac_l, jac_r

    real(dp) :: rho, rhoinv, u, p, rhoet, a, m, fa, fb

    real(dp), dimension(3)   :: drho_dq, du_dq, dp_dq, drhoet_dq
    real(dp), dimension(3)   :: da_dq, dm_dq, dfa_dq, dfb_dq
    real(dp), dimension(3,3) :: mat

    continue

! handle left side first

! calculate primitive vars from conserved
    rho = qL(1)
    rhoinv = one/rho
    u = qL(2)
    p = qL(3)
    rhoet = p*xgm1 + half * rho * u**2
    a = speed_of_sound(p, rho)
    m = u/a

! linearization of primitive variables wrt conserved
    drho_dq(1) = one
    drho_dq(2) = zero
    drho_dq(3) = zero

    du_dq(1) = zero
    du_dq(2) = one
    du_dq(3) = zero

    dp_dq(1) = zero
    dp_dq(2) = zero
    dp_dq(3) = one

    drhoet_dq(1) = half*u**2
    drhoet_dq(2) = rho*u
    drhoet_dq(3) = xgm1

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

!    Flux = [ fa, fa*fb*xg, half*fa*fb*fb*xg2m1 ]

! now, formulate jacobian

    mat = zero

    if( abs(m) < 1.0_dp ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))
    elseif( m >= 1.0_dp ) then
      mat(1,:) = drho_dq(:)*u + rho*du_dq(:)
      mat(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
      mat(3,:) = ( drhoet_dq(:) + dp_dq(:) )*u + (rhoet + p)*du_dq(:)
    endif

    jac_l = mat

! handle right side second

! calculate primitive vars from conserved
    rho = qR(1)
    rhoinv = one/rho
    u = qR(2)
    p = qR(3)
    rhoet = p*xgm1 + half * rho * u**2
    a = speed_of_sound(p,rho)
    m = u/a

! linearization of right primitive variables wrt conserved
    drho_dq(1) = one
    drho_dq(2) = zero
    drho_dq(3) = zero

    du_dq(1) = zero
    du_dq(2) = one
    du_dq(3) = zero

    dp_dq(1) = zero
    dp_dq(2) = zero
    dp_dq(3) = one

    drhoet_dq(1) = half*u**2
    drhoet_dq(2) = rho*u
    drhoet_dq(3) = xgm1

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

    dfb_dq(:) = (gm1*m - two)*da_dq(:) + gm1*a*dm_dq(:)

! now, form jacobian

    mat = zero

    if( abs(m) < 1.0_dp ) then
      mat(1,:) = dfa_dq(:)
      mat(2,:) = (dfa_dq(:)*fb + dfb_dq(:)*fa)*xg
      mat(3,:) = half*xg2m1*(dfa_dq(:)*fb*fb + two*fa*fb*dfb_dq(:))
    elseif( m <= -1.0_dp ) then
      mat(1,:) = drho_dq(:)*u   + rho*du_dq(:)
      mat(2,:) = drho_dq(:)*u*u + two*rho*u*du_dq(:) + dp_dq(:)
      mat(3,:) = ( drhoet_dq(:) + dp_dq(:) )*u + (rhoet + p)*du_dq(:)
    endif

    jac_r = mat

  end subroutine jac_vanleer_q_1D

!============================== jac_source_q_1D ==============================80
!
! This subroutine returns the source Jacobian wrt conserved variables
!
!=============================================================================80
  subroutine jac_source_q_1D( dadx_cc, cell_jac, source_jac )

    use set_precision, only : dp
    use set_constants, only : zero

    real(dp),                 intent(in)  :: dadx_cc, cell_jac
    real(dp), dimension(3,3), intent(out) :: source_jac

    continue

    source_jac = zero
    source_jac(2,3) = dadx_cc*cell_jac

!    source_jac(1,1) = zero
!    source_jac(2,1) = zero
!    source_jac(3,1) = zero
!    source_jac(1,2) = zero
!    source_jac(2,2) = zero
!    source_jac(3,2) = zero
!    source_jac(1,3) = zero
!    source_jac(2,3) = dadx_cc*cell_jac
!    source_jac(3,3) = zero

  end subroutine jac_source_q_1D

!============================== jac_source_q_1D ==============================80
!
! This subroutine returns the source Jacobian wrt conserved variables
!
!=============================================================================80
  subroutine dconserved_dprimitive( q, dcdp )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one
    use fluid_constants, only : xgm1

    real(dp), dimension(3),   intent(in)  :: q
    real(dp), dimension(3,3), intent(out) :: dcdp

    continue

    dcdp(1,1) = one
    dcdp(2,1) = q(2)
    dcdp(3,1) = half*q(2)*q(2)
    dcdp(1,2) = zero
    dcdp(2,2) = q(1)
    dcdp(3,2) = q(1)*q(2)
    dcdp(1,3) = zero
    dcdp(2,3) = zero
    dcdp(3,3) = xgm1

  end subroutine dconserved_dprimitive

  include 'speed_of_sound.f90'

end module jacobians
