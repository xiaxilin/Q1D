!========================= van_leer_jacobian ==================================80
!
! This subroutine takes a face normal vector, left and right conserved variables
! and returns the left and right flux jacobians
!
!==============================================================================80
subroutine van_leer_jacobian( normal_vec, ql, qr, jac_l, jac_r )

  use set_precision, only : dp
  use constants,     only : zero, fourth, one, two, four
  use flow_vars,     only : gamma

  implicit none

  real(dp), dimension(3),   intent(in)  :: normal_vec
  real(dp), dimension(5),   intent(in)  :: ql, qr
  real(dp), dimension(5,5), intent(out) :: jac_l, jac_r

  real(dp) :: gm1
  real(dp) :: nx, ny, nz, length
  real(dp) :: rho, rhoinv, u, v, w, p, h, a, unorm , m

  real(dp), dimension(5) :: drho_dq, du_dq, dv_dq, dw_dq, dp_dq, dq5_dq
  real(dp), dimension(5) :: da_dq, dunorm_dq, dm_dq, dfa_dq, dfb_dq
  real(dp), dimension(5,5) :: mat

  continue

  gm1 = gamma - one

! rescale normal vector to get proper normal velocity later
  nx = normal_vec(1)
  ny = normal_vec(2)
  nz = normal_vec(3)
        
  length = sqrt(nx**2 + ny**2 + nz**2)

  nx = nx/area
  ny = ny/area
  nz = nz/area

! handle left side first

! calculate primitive vars from conserved
  rho = ql(1)
  rhoinv = one/rho
  u = ql(2)*rhoinv
  v = ql(3)*rhoinv
  w = ql(4)*rhoinv
  p = gm1*( ql(5) - half * rho * (u**2 + v**2 + w**2) )
  h = (ql(5) + p)*rhoinv
  a = sqrt(gamma*p*rhoinv)

! normal velocity
  unorm = u*nx + v*ny + w*nz

! linearization of primitive variables
  drho_dq = zero
  du_dq   = zero
  dv_dq   = zero
  dw_dq   = zero
  dp_dq   = zero
  dq5_dq  = zero

  drho_dq(1) = one
       
  du_dq(1) = -u*rhoinv
  du_dq(2) = rhoinv
        
  dv_dq(1) = -v*rhoinv
  dv_dq(3) = rhoinv

  dw_dq(1) = -w*rhoinv
  dw_dq(4) = rhoinv

  dp_dq(1) = half*gm1*(u**2 + v**2 + w**2)
  dp_dq(2) = -gm1*u
  dp_dq(3) = -gm1*v
  dp_dq(4) = -gm1*w
  dp_dq(5) = gm1

  dq5_dq(5) = one

  da_dq(:) = (half*gamma/a) * ( rhoinv*dp_dq(:) - (p*rhoinv**2)*drho_dq(:) )

! face normal velocity and mach number linearization
  dunorm_dq(:) = du_dq(:)*nx + dv_dq(:)*ny + dw_dq(:)*nz

  m = unorm/a

  dm_dq(:) = dunorm_dq(:)/a - (unorm/a**2)*da_dq(:)

! linearization of FVS mass term
  fa =  fourth*rho*a*(m+one)**2

  dfa_dq(:) =  fourth*( drho_dq(:)*a*(m+one)**2 
            + rho*da_dq(:)*(m+one)**2 + rho*a*two*(m+one)*dm_dq(:) ) 

! linearization of FVS energy term

  fb = -gm1*unorm**2 + two*a**2 + two*gm1*unorm*a
  fb = fb/(gamma**2 - one) + half*( u**2 + v**2 +w**2 )
  fb = fb*fa

  dfb_dq(:) = -gm1*two*unorm*dunorm_dq(:) + four*a*da_dq(:)                     &
            + two*gm1*(dunorm_dq(:)*a + unorm*da_dq(:))
  dfb_dq(:) = dfb_dq(:)/(gamma**2 - one) + u*du_dq(:) + v*dv_dq(:) + w*dw_dq(:)
  dfb_dq(:) = dfb_dq(:)*fa + fb*dfa_dq(:)/fa

! now, formulate jacobian

  mat = zero

  if( abs(m)<1.0 ) then

    mat(1,:) = dfa_dq(:)
    mat(2,:) = dfa_dq(:)*( u + nx*(-unorm + two*a)/gamma )                      &
             + fa*( du_dq(:) + nx*( -dunorm_dq(:) + two*da_dq(:) )/gamma )
    mat(3,:) = dfa_dq(:)*( v + ny*(-unorm + two*a)/gamma )                      &
             + fa*( dv_dq(:) + ny*( -dunorm_dq(:) + two*da_dq(:) )/gamma )
    mat(4,:) = dfa_dq(:)*( w + nz*(-unorm + two*a)/gamma )                      &
             + fa*( dw_dq(:) + nz*( -dunorm_dq(:) + two*da_dq(:) )/gamma )
    mat(5,:) = dfb_dq(:)

  elseif( m>=1.0 ) then

    mat(1,:) = drho_dq(:)*unorm   + rho*dunorm_dq(:)
    mat(2,:) = drho_dq(:)*unorm*u + rho*dunorm_dq(:)*u + rho*unorm*du_dq(:)
    mat(3,:) = drho_dq(:)*unorm*v + rho*dunorm_dq(:)*v + rho*unorm*dv_dq(:)
    mat(4,:) = drho_dq(:)*unorm*w + rho*dunorm_dq(:)*w + rho*unorm*dw_dq(:)

    mat(2,:) = mat(2,:) + nx*dp_dq(:)
    mat(3,:) = mat(3,:) + ny*dp_dq(:)
    mat(4,:) = mat(4,:) + nz*dp_dq(:)

    mat(5,:) = ( dq5_dq(:) + dp_dq(:) )*unorm + (ql(5) + p)*dunorm_dq(:)

  endif

  jac_l = mat

! handle right side second

! calculate primitive vars from conserved
  rho = qr(1)
  rhoinv = one/rho
  u = qr(2)*rhoinv
  v = qr(3)*rhoinv
  w = qr(4)*rhoinv
  p = gm1*( qr(5) - half * rho * (u**2 + v**2 + w**2) )
  h = (qr(5) + p)*rhoinv
  a = sqrt(gamma*p*rhoinv)

! normal velocity
  unorm = u*nx + v*ny + w*nz

! linearization of right primitive variables
  drho_dq = zero
  du_dq   = zero
  dv_dq   = zero
  dw_dq   = zero
  dp_dq   = zero
  dq5_dq  = zero

  drho_dq(1) = one
       
  du_dq(1) = -u*rhoinv
  du_dq(2) = rhoinv
        
  dv_dq(1) = -v*rhoinv
  dv_dq(3) = rhoinv

  dw_dq(1) = -w*rhoinv
  dw_dq(4) = rhoinv

  dp_dq(1) =  half*gm1*(u**2 + v**2 + w**2)
  dp_dq(2) = -gm1*u
  dp_dq(3) = -gm1*v
  dp_dq(4) = -gm1*w
  dp_dq(5) =  gm1

  dq5_dq(5) = one

  da_dq(:) = (half*gamma/a) * ( rhoinv*dp_dq(:) - (p*rhoinv**2)*drho_dq(:) )

! face normal velocity and mach number linearization
  dunorm_dq(:) = du_dq(:)*nx + dv_dq(:)*ny + dw_dq(:)*nz

  m = unorm/a

  dm_dq(:) = dunorm_dq(:)/a - (unorm/a**2)*da_dq(:)

! linearization of FVS mass term
  fa = -fourth*rho*a*(m-one)**2

  dfa_dq(:) = -fourth*( drho_dq(:)*a*(m-one)**2 
            + rho*da_dq(:)*(m-one)**2 + rho*a*two*(m-one)*dm_dq(:) )

! linearization of FVS energy term

  fb = -gm1*unorm**2 + two*a**2 - two*gm1*unorm*a
  fb = fb/(gamma**2 - one) + half*( u**2 + v**2 + w**2 )
  fb = fb*fa

  dfb_dq(:) = -gm1*two*unorm*dunorm_dq(:) + four*a*da_dq(:)                     &
            - two*gm1*(dunorm_dq(:)*a + unorm*da_dq(:)) 
  dfb_dq(:) = dfb_dq(:)/(gamma**2 - one) + u*du_dq(:) + v*dv_dq(:) + w*dw_dq(:)
  dfb_dq(:) = dfb_dq(:)*fa + fb*dfa_dq(:)/fa

! now, form jacobian

  mat = zero

  if( abs(m)<1.0 ) then

    mat(1,:) = dfa_dq(:)
    mat(2,:) = dfa_dq(:) * ( u + nx*(-unorm - two*a)/gamma )  &
             + fa * ( du_dq(:) + nx*( -dunorm_dq(:) - two*da_dq(:) )/gamma )
    mat(3,:) = dfa_dq(:) * ( v + ny*(-unorm - two*a)/gamma )  &
             + fa * ( dv_dq(:) + ny*( -dunorm_dq(:) - two*da_dq(:) )/gamma )
    mat(4,:) = dfa_dq(:) * ( w + nz*(-unorm - two*a)/gamma )  &
             + fa * ( dw_dq(:) + nz*( -dunorm_dq(:) - two*da_dq(:) )/gamma )
    mat(5,:) = dfb_dq(:)

  elseif( m<=-1.0 ) then

    mat(1,:) = drho_dq(:)*unorm   + rho*dunorm_dq(:)
    mat(2,:) = drho_dq(:)*unorm*u + rho*dunorm_dq(:)*u + rho*unorm*du_dq(:)
    mat(3,:) = drho_dq(:)*unorm*v + rho*dunorm_dq(:)*v + rho*unorm*dv_dq(:)
    mat(4,:) = drho_dq(:)*unorm*w + rho*dunorm_dq(:)*w + rho*unorm*dw_dq(:)

    mat(2,:) = mat(2,:) + nx*dp_dq(:)
    mat(3,:) = mat(3,:) + ny*dp_dq(:)
    mat(4,:) = mat(4,:) + nz*dp_dq(:)

    mat(5,:) = ( dq5_dq(:) + dp_dq(:) )*unorm + (qr(5) + p)*dunorm_dq(:)

  endif

  jac_r = mat

end subroutine van_leer_jacobian
