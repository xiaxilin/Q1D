! Will hold implicit BC's for Q1D Nozzle
! Subsonic inflow with Po and To specified, delta Po and delta To = 0
! Supersonic inflow with any variables specified, they get converted to cons.
! Subsonic outflow with p specified... ie delta p = 0
! Supersonic outflow... linear extrapolation from interior

module bc

  implicit none

  private

  public :: subsonic_inflow_explicit
  public :: outflow_explicit
  public :: subsonic_inflow
  public :: set_outflow
  public :: subsonic_inflow_prim
  public :: set_outflow_prim

contains

!========================== subsonic_inflow_explicit =========================80
!
! Uses primitive variables
!
!=============================================================================80
  subroutine subsonic_inflow_explicit( cc_in, cc_1, cc_2, cc_3 )

    use set_precision,   only : dp
    use set_constants,   only : third, half, one, onep5, two, three, four, six
    use fluid_constants, only : r, gamma, gm1, xgm1, gxgm1!, gm1xgp1, gp1
    use initialize_soln, only : po, to
    use residual,        only : inflow_gc

    real(dp), dimension(3), intent(in)  :: cc_1, cc_2, cc_3
    real(dp), dimension(3), intent(out) :: cc_in

    real(dp) :: vel_max, u_extrap, psi

! testing
!    real(dp) :: rho0, h0, s0, p, u, temp, a2, ao2

    continue

! set max, physically possible velocity
    vel_max = sqrt(two*gamma*r*to*xgm1)-one

! extrapolate velocity and limit
    select case( inflow_gc )
    case( 1 ) ! extrapolate to ghost cell (zero grad)
      u_extrap = third*( four*cc_1(2) - cc_2(2) )
    case( 2 ) ! extrapolate to ghost cell (zero curv)
      u_extrap = two*cc_1(2) - cc_2(2)
    case( 3 )
      u_extrap = ( three*cc_1(2) - three*cc_2(2) + cc_3(2) )
    case ( -1 ) ! extrapolate to face
      u_extrap = onep5*cc_1(2) - half*cc_2(2)
    case ( -2 ) ! zero grad at face
      u_extrap = ( 7._dp*cc_1(2) - cc_2(2) ) / six
    case ( -3 ) !extrap to face
      u_extrap = ( 11._dp*cc_1(2) - 7._dp*cc_2(2) + two*cc_3(2) ) / six
    end select

    cc_in(2) = max(-vel_max, min(u_extrap, vel_max))

! now calculate inflow
    psi = to/(to-(gm1*cc_in(2)**2/(two*gamma*r)))

    cc_in(1) = po/(r*to*psi**xgm1)
    cc_in(3) = po/psi**gxgm1

! to test another version...
!    ao2 = gamma*r*to
!    a2 = ao2 - gm1*half*cc_in(2)**2
!    temp = ao2/a2

!    rho0 = po/(r*to)
!    cc_in(1) = rho0*temp**gm1
!    cc_in(3) = po*temp**(gm1/gamma)

!    rho0 = po/(r*to)
!    h0 = gxgm1*po/rho0
!    s0 = r*xgm1*log(po/rho0**gamma)

!    cc_in(2) = min(max(two*cc_1(1) - cc_2(1),0.0001_dp),rho0)

!    p = exp(s0*gm1/r)*cc_in(1)**gamma
!    u = sqrt(two*(h0 - gxgm1*p/cc_in(1)))

!    cc_in(2) = u
!    cc_in(3) = p

! floor variables
    cc_in = floor_primitive_vars(cc_in)

  end subroutine subsonic_inflow_explicit

!========================= subsonic_inflow_r_explicit ========================80
!
! Uses primitive variables
!
!=============================================================================80
  subroutine subsonic_inflow_r_explicit( cc_in, cc_1 )

    use set_precision,   only : dp
    use set_constants,   only : one, half, two
    use fluid_constants, only : r, gm1, xgm1, gxgm1, gm1xgp1, gp1
    use initialize_soln, only : po, to

    real(dp), dimension(3), intent(in)  :: cc_1
    real(dp), dimension(3), intent(out) :: cc_in

    real(dp) :: a_1, a_in, a_bound, r_minus, t_in

    continue

    a_1 = speed_of_sound(cc_1(3), cc_1(1))

    a_in = sqrt(half*gm1*cc_1(2)**2+a_1**2)

    r_minus = -cc_1(2)-two*xgm1*a_1

    a_bound = -r_minus*gm1xgp1                                                 &
            * (one-sqrt((gp1*a_in**2)/(gm1*r_minus**2)-half*gm1))

    t_in = to*(a_bound**2/a_in**2)

    cc_in(3) = po*(t_in/to)**gxgm1

    cc_in(1) = cc_in(3)/(r*t_in)

    cc_in(2) = sqrt(2012.0_dp*(to-t_in))

! floor variables
    cc_in = floor_primitive_vars(cc_in)

  end subroutine subsonic_inflow_r_explicit

!============================= outflow_explicit ==============================80
!
! Uses primitive variables
!
!=============================================================================80
  subroutine outflow_explicit( cc_out, cc_1, cc_2, cc_3 )

    use set_precision,   only : dp
    use set_constants,   only : third, half, onep5, two, three, four, six
    use initialize_soln, only : pback
    use residual,        only : outflow_gc

    real(dp), dimension(3), intent(in)  :: cc_1, cc_2, cc_3
    real(dp), dimension(3), intent(out) :: cc_out

    continue

! extrapolate all variables
    select case( outflow_gc )
    case( 1 ) ! extrapolate to ghost cell (zero grad)
      cc_out = third*( four*cc_1 - cc_2 )
    case( 2 ) ! extrapolate to ghost cell (zero curv)
      cc_out = two*cc_1 - cc_2
    case( 3 )
      cc_out = ( three*cc_1 - three*cc_2 + cc_3 )
    case ( -1 ) ! extrapolate to face
      cc_out = onep5*cc_1 - half*cc_2
    case ( -2 ) ! zero grad at face
      cc_out = ( 7._dp*cc_1 - cc_2 ) / six
    case ( -3 ) !extrap to face
      cc_out = ( 11._dp*cc_1 - 7._dp*cc_2 + two*cc_3 ) / six
    end select

! set back pressure if appropriate
    if ( pback >= 0.0_dp ) cc_out(3) = pback

! floor variables
    cc_out = floor_primitive_vars(cc_out)

  end subroutine outflow_explicit

!============================== subsonic_inflow ==============================80
!
! Uses conserved variables with Po and To specified to calculate subsonic inflow
! FIXME: Check directions on extrapolations, could account for adjoint wiggles
!
!=============================================================================80
  subroutine subsonic_inflow( iter, cc_in, cc_1, cc_2, DD, DU1, DU2, RHS )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two, three, four
    use fluid_constants, only : gamma, gm1, gxgm1, xgm1, R, cv
    use initialize_soln, only : po, to
    use residual,        only : firstorder
    use lhs,             only : lhs_order

    integer,                  intent(in)  :: iter
    real(dp), dimension(3),   intent(in)  :: cc_in, cc_1, cc_2
    real(dp), dimension(3,3), intent(out) :: DD, DU1, DU2
    real(dp), dimension(3),   intent(out) :: RHS

    real(dp) :: rho, u, rhoet, p, T, m, factor
    real(dp) :: dTdrho, dTdrhou, dTdrhoet
    real(dp) :: dPdrho, dPdrhou, dPdrhoet
    real(dp) :: dM2drho, dM2drhou, dM2drhoet

    continue

    DD  = zero
    DU1 = zero
    DU2 = zero

! Extrapolate velocity from interior
    DD(2,1) = -cc_in(2)/cc_in(1)**2
    DD(2,2) = one/cc_in(1)

    DU1(2,1) = -cc_1(2)/cc_1(1)**2
    DU1(2,2) = one/cc_1(1)

    DU1 = -DU1

    DU2(2,1) = -cc_2(2)/cc_2(1)**2
    DU2(2,2) = one/cc_2(1)

! Now need to account for delta Po = delta To = 0 at inflow... DD matrix

    rho   = cc_in(1)
    u     = cc_in(2)/cc_in(1)
    rhoet = cc_in(3)
    p     = gm1*(rhoet-half*rho*u**2)
    T     = p/(rho*R)

    m = u/speed_of_sound(p, rho)

    factor = one + half*gm1*m**2

    dTdrho   = -(rhoet/rho**2 - u**2/rho)/cv
    dTdrhou  = -u/(rho*cv)
    dTdrhoet = one/(rho*cv)

    dM2drho   = -m**2*(two/rho + dTdrho/T)
    dM2drhou  = two*u/(gamma*p) - dTdrhou*m**2/T
    dM2drhoet = -dTdrhoet*m**2/T

! add To equations...
    DD(1,1) = factor*dTdrho   + half*gm1*T*dM2drho
    DD(1,2) = factor*dTdrhou  + half*gm1*T*dM2drhou
    DD(1,3) = factor*dTdrhoet + half*gm1*T*dM2drhoet

    dPdrho   = half*gm1*u**2
    dPdrhou  = -gm1*u
    dPdrhoet = gm1

! add Po equations...
    DD(3,1) = factor**gxgm1*dPdrho   + half*gamma*factor**xgm1*P*dM2drho
    DD(3,2) = factor**gxgm1*dPdrhou  + half*gamma*factor**xgm1*P*dM2drhou
    DD(3,3) = factor**gxgm1*dPdrhoet + half*gamma*factor**xgm1*P*dM2drhoet

    RHS(1) = to - T*factor
    RHS(2) = -cc_in(2)/cc_in(1) + cc_1(2)/cc_1(1)! - cc_2(2)/cc_2(1)
    RHS(3) = po - p*factor**gxgm1

    if ( iter >= firstorder .and. lhs_order /= 1 ) then
      DD(2,:) = three*DD(2,:)
      DU1 = four*DU1

      RHS(2) = -three*cc_in(2)/cc_in(1) + four*cc_1(2)/cc_1(1) - cc_2(2)/cc_2(1)
    end if

  end subroutine subsonic_inflow

!============================ subsonic_inflow_prim ===========================80
!
! Uses primitive variables with Po and To specified to calculate subsonic inflow
! FIXME: Check directions on extrapolations, could account for adjoint wiggles
!
!=============================================================================80
  subroutine subsonic_inflow_prim( iter, cc_in, cc_1, cc_2, DD, DU1, DU2, RHS )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two, three, four
    use fluid_constants, only : gamma, gm1, gxgm1, xgm1, R
    use initialize_soln, only : po, to
    use residual,        only : firstorder
    use lhs,             only : lhs_order

    integer,                  intent(in)  :: iter
    real(dp), dimension(3),   intent(in)  :: cc_in, cc_1, cc_2
    real(dp), dimension(3,3), intent(out) :: DD, DU1, DU2
    real(dp), dimension(3),   intent(out) :: RHS

    real(dp) :: rho, u, p, T, m, factor

    real(dp), dimension(3) :: dTdq, dM2dq, dpdq

    continue

    DD  = zero
    DU1 = zero
    DU2 = zero

! Extrapolate velocity from interior
    DD(2,2)  = one

    DU1(2,2) = -one

    DU2(2,2) = one

! Now need to account for delta Po = delta To = 0 at inflow... DD matrix

    rho = cc_in(1)
    u   = cc_in(2)
    p   = cc_in(3)
    T   = p/(rho*R)

    m = u/speed_of_sound(p, rho) ! u/sqrt(gamma*P/rho)

    factor = one + half*gm1*m**2

    dTdq(1) = -p/(R*rho**2)
    dTdq(2) = zero
    dTdq(3) = one/(rho*R)

! M**2 = (rho*u**2) / (gamma*P)
    dM2dq(1) = u**2/(gamma*P)
    dM2dq(2) = two*rho*u/(gamma*P)
    dM2dq(3) = -rho*u**2 / (gamma*P**2)

! add To equations...
    DD(1,:) = factor*dTdq + half*gm1*T*dM2dq

    dpdq(1) = zero
    dpdq(2) = zero
    dpdq(3) = one

! add Po equations...
    DD(3,:) = factor**gxgm1*dPdq + half*gamma*factor**xgm1*p*dM2dq

    RHS(1) = to - T*factor
    RHS(2) = -cc_in(2) + cc_1(2)! - cc_2(2)
    RHS(3) = po - p*factor**gxgm1

    if ( iter >= firstorder .and. lhs_order /= 1 ) then
      DD(2,2)  = three
      DU1(2,2) = -four

      RHS(2) = -three*cc_in(2) + four*cc_1(2) - cc_2(2)
    end if

  end subroutine subsonic_inflow_prim

!================================ set_outflow ================================80
!
! Creates matrices for variable extrapolation
!
!=============================================================================80
  subroutine set_outflow( iter, cc_out, cc_1, cc_2, DD, DL1, DL2, RHS )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, three, four
    use fluid_constants, only : gm1
    use initialize_soln, only : pback
    use residual,        only : firstorder
    use lhs,             only : lhs_order

    integer,                  intent(in)  :: iter
    real(dp), dimension(3),   intent(in)  :: cc_out, cc_1, cc_2
    real(dp), dimension(3,3), intent(out) :: DD, DL1, DL2
    real(dp), dimension(3),   intent(out) :: RHS

    real(dp) :: u_out, u_1, u_2
    real(dp) :: p_out, p_1, p_2

    continue

    u_out = cc_out(2) / cc_out(1)
    u_1   = cc_1(2)   / cc_1(1)
    u_2   = cc_2(2)   / cc_2(1)

    p_out = gm1*(cc_out(3) - half*cc_out(1)*u_out**2)
    p_1   = gm1*(cc_1(3)   - half*cc_1(1)*u_1**2)
    p_2   = gm1*(cc_2(3)   - half*cc_2(1)*u_2**2)

! Extrapolate velocity from interior
    DD(1,1) = one
    DD(2,1) = -u_out/cc_out(1)
    DD(3,1) = half*gm1*u_out**2
    DD(1,2) = zero
    DD(2,2) = one/cc_out(1)
    DD(3,2) = -gm1*u_out
    DD(1,3) = zero
    DD(2,3) = zero
    DD(3,3) = gm1

    DL1(1,1) = one
    DL1(2,1) = -u_1/cc_1(1)
    DL1(3,1) = half*gm1*u_1**2
    DL1(1,2) = zero
    DL1(2,2) = one/cc_1(1)
    DL1(3,2) = -gm1*u_1
    DL1(1,3) = zero
    DL1(2,3) = zero
    DL1(3,3) = gm1

    DL1 = -DL1

    DL2(1,1) = one
    DL2(2,1) = -u_2/cc_2(1)
    DL2(3,1) = half*gm1*u_2**2
    DL2(1,2) = zero
    DL2(2,2) = one/cc_2(1)
    DL2(3,2) = -gm1*u_2
    DL2(1,3) = zero
    DL2(2,3) = zero
    DL2(3,3) = gm1

    RHS(1) = -cc_out(1) + cc_1(1)! - cc_2(1)
    RHS(2) = -u_out     + u_1!     - u_2
    RHS(3) = -p_out     + p_1!     - p_2

    if ( iter >= firstorder .and. lhs_order /= 1 ) then
      DD  = three*DD
      DL1 = four*DL1

      RHS(1) = -three*cc_out(1) + four*cc_1(1) - cc_2(1)
      RHS(2) = -three*u_out     + four*u_1     - u_2
      RHS(3) = -three*p_out     + four*p_1     - p_2
    end if

    if (pback > zero) then
      if ( iter >= firstorder .and. lhs_order /= 1 ) DD(3,:) = DD(3,:)/three
      DL1(3,:) = zero
      DL2(3,:) = zero
      RHS(3)   = pback - p_out
    end if

  end subroutine set_outflow

!============================== set_outflow_prim =============================80
!
! Creates matrices for variable extrapolation wrt primitive variables
!
!=============================================================================80
  subroutine set_outflow_prim( iter, cc_out, cc_1, cc_2, DD, DL1, DL2, RHS )

    use set_precision,   only : dp
    use set_constants,   only : zero, one, three, four
    use initialize_soln, only : pback
    use residual,        only : firstorder
    use lhs,             only : lhs_order

    integer,                  intent(in)  :: iter
    real(dp), dimension(3),   intent(in)  :: cc_out, cc_1, cc_2
    real(dp), dimension(3,3), intent(out) :: DD, DL1, DL2
    real(dp), dimension(3),   intent(out) :: RHS

    real(dp), dimension(3,3) :: ident3x3

    continue

    ident3x3 = reshape( [one, zero, zero, zero, one, zero, zero, zero, one],   &
                        [3,3] )
    DD  = ident3x3
    DL1 = -ident3x3
    DL2 = ident3x3

    RHS = -cc_out + cc_1! - cc_2

    if ( iter >= firstorder .and. lhs_order /= 1 ) then
      DD  = three*DD
      DL1 = four*DL1

      RHS = -three*cc_out + four*cc_1 - cc_2
    end if

    if ( pback > zero ) then
      if ( iter >= firstorder .and. lhs_order /= 1 ) DD(3,3) = one
      DL1(3,3) = zero
      DL2(3,3) = zero

      RHS(3)   = pback - cc_out(3)
    end if

  end subroutine set_outflow_prim

  include 'speed_of_sound.f90'
  include 'floor_primitive_vars.f90'

end module bc
