subroutine roe_jacobian( normal_vec, q_l, q_r, jac_l, jac_r )

  use set_precision, only : dp
  use constants,     only : zero, half, one
  use flow_vars,     only : gamma

  implicit none

  real*8  :: nxyz(3)
  real*8  :: ql(5),qr(5)
  real*8  :: lmat(5,5),rmat(5,5)
  real*8  :: gam
  integer :: imode

  real*8  :: gm1
  real*8  :: nxd,nyd,nzd,area,nx,ny,nz
  real*8  :: rol,ul,vl,wl,pl,hl
  real*8  :: ror,ur,vr,wr,pr,hr
  real*8  :: uconl,uconr
  real*8  :: ubar,vbar,wbar,hbar,uconbar,cbar,robar
  real*8  :: dpress,dro,du,dv,dw
  real*8  :: eig1,eig2,eig3

  real*8  :: fact,A,B,term1,term2,del1,del2

  integer :: k,i,j
  real*8  :: dro_dql(5),dro_dqr(5)
  real*8  :: du_dql(5),du_dqr(5)
  real*8  :: dv_dql(5),dv_dqr(5)
  real*8  :: dw_dql(5),dw_dqr(5)
  real*8  :: dp_dql(5),dp_dqr(5)
  real*8  :: ducon_dql(5),ducon_dqr(5)
  real*8  :: ddel1_dql(5),ddel1_dqr(5)
  real*8  :: ddel2_dql(5),ddel2_dqr(5)

  real*8  :: dq5_dql(5),dq5_dqr(5)
  real*8  :: dh_dql(5),dh_dqr(5)
  real*8  :: dfact_dql(5),dfact_dqr(5)
  real*8  :: dA_dql(5),dA_dqr(5)
  real*8  :: dB_dql(5),dB_dqr(5)
  real*8  :: drobar_dql(5),dubar_dql(5),dvbar_dql(5),dwbar_dql(5)
  real*8  :: drobar_dqr(5),dubar_dqr(5),dvbar_dqr(5),dwbar_dqr(5)
  real*8  :: dhbar_dql(5),duconbar_dql(5),dcbar_dql(5)
  real*8  :: dhbar_dqr(5),duconbar_dqr(5),dcbar_dqr(5)

  real*8  :: deig1_dql(5),deig2_dql(5),deig3_dql(5)
  real*8  :: deig1_dqr(5),deig2_dqr(5),deig3_dqr(5)
  real*8  :: dterm1_dql(5),dterm1_dqr(5)
  real*8  :: dterm2_dql(5),dterm2_dqr(5)
  real*8  :: imat(5,5)
  real*8  :: cl,cr,dc_dql(5),dc_dqr(5)
  real*8  :: t1a,t1b,t2a,t2b,t3a,t3b
  real*8  :: eps1,eps2,eps3

  real*8  :: lmat1(5,5),rmat1(5,5)
!------------------------------------------------------------------------------

  continue

  gm1 = gamma - one

  nxd = nxyz(1)
  nyd = nxyz(2)
  nzd = nxyz(3)

  area = sqrt(nxd*nxd + nyd*nyd + nzd*nzd)

  nx = nxd/area
  ny = nyd/area
  nz = nzd/area

! calculate primitive variables

  rhol = ql(1)
  rhoinvl = one/rhol
  ul  = ql(2)*rhoinvl
  vl  = ql(3)*rhoinvl
  wl  = ql(4)*rhoinvl
  pl  = gm1*( ql(5) - half * rhol * (ul**2 + vl**2 + wl**2) )
  hl  = (ql(5) + pl)*rhoinvl
  al  = sqrt(gamma*pl*rhoinvl)
        
  rhor = qr(1)
  rhoinvr = one/rhor
  ur  = qr(2)*rhoinvr
  vr  = qr(3)*rhoinvr
  wr  = qr(4)*rhoinvr
  pr  = gm1*( qr(5) - half * rhor * (ur**2 + vr**2 + wr**2) )
  hr  = (qr(5) + pr)*rhoinvr
  ar  = sqrt(gamma*pr*rhoinvr)

!----> face normal velocities

  uconl = ul*nx + vl*ny + wl*nz
  uconr = ur*nx + vr*ny + wr*nz

! jump in primitive vars across face

  drho   = rhor - rhol
  du     = ur - ul
  dv     = vr - vl
  dw     = wr - wl
  dpress = pr - pl

! Roe face variables

  fact = sqrt(rhor/rhol)

  A    = one /(one + fact)
  B    = fact/(one + fact)

  rhobar = rhol*fact
  ubar   = ul*A + ur*B
  vbar   = vl*A + vr*B
  wbar   = wl*A + wr*B
  hbar   = hl*A + hr*B
  abar   = sqrt(gm1*(hbar - half*(ubar**2 + vbar**2 + wbar**2)))
  uconbar = ubar*nx + vbar*ny + wbar*nz

! Roe face eigenvalues
  eig1 = abs(uconbar)
  eig2 = abs(uconbar + abar)
  eig3 = abs(uconbar - abar)

!------------------------------------------------------------------------------!
!-------> linearization of left and right primitive states <-------------------!
!------------------------------------------------------------------------------!
! left state

  dro_dql = zero
  du_dql  = zero
  dv_dql  = zero
  dw_dql  = zero
  dq5_dql = zero

  dro_dql(1) = one

  du_dql(1) = -rhoinvl*ul
  du_dql(2) =  rhoinvl

  dv_dql(1) = -rhoinvl*vl
  dv_dql(3) =  rhoinvl

  dw_dql(1) = -rhoinvl*wl
  dw_dql(4) =  rhoinvl

  dp_dql(1) =  half*gm1*( ul**2 + vl**2 + wl**2 )
  dp_dql(2) = -gm1*ul
  dp_dql(3) = -gm1*vl
  dp_dql(4) = -gm1*wl
  dp_dql(5) =  gm1

  dq5_dql(5) = one

  dh_dql(:) = -(ql(5)+pl)*drho_dql(:)*rhoinvl**2 + rhoinvl*(dq5_dql(:)+dp_dql(:))
 
  da_dql(:) = (half*gamma/al)*( rhoinvl*dp_dql(:) - (pl*rhoinvl**2)*drho_dql(:) )

  ducon_dql(1) = -rhoinvl*uconl
  ducon_dql(2) =  rhoinvl*nx
  ducon_dql(3) =  rhoinvl*ny
  ducon_dql(4) =  rhoinvl*nz
  ducon_dql(5) =  zero

! right state

  dro_dqr = zero
  du_dqr  = zero
  dv_dqr  = zero
  dw_dqr  = zero
  dq5_dqr = zero

  dro_dqr(1) = one

  du_dqr(1) = -rhoinvr*ur
  du_dqr(2) =  rhoinvr

  dv_dqr(1) = -rhoinvr*vr
  dv_dqr(3) =  rhoinvr

  dw_dqr(1) = -rhoinvr*wr
  dw_dqr(4) =  rhoinvr

  dp_dqr(1) =  half*gm1*( ur**2 + vr**2 + wr**2 )
  dp_dqr(2) = -gm1*ur
  dp_dqr(3) = -gm1*vr
  dp_dqr(4) = -gm1*wr
  dp_dqr(5) =  gm1

  dq5_dqr(5) = one

  dh_dqr(:) = -(qr(5)+pr)*dro_dqr(:)*rhoinvr**2 + rhoinvr*(dq5_dqr(:)+dp_dqr(:))

  da_dqr(:) = (half*gamma/ar)*( rhoinvr*dp_dqr(:) - (pr*rhoinv**2)*drho_dqr(:) )

  ducon_dqr(1) = -rhoinvr*uconr
  ducon_dqr(2) =  rhoinvr*nx
  ducon_dqr(3) =  rhoinvr*ny
  ducon_dqr(4) =  rhoinvr*nz
  ducon_dqr(5) =  zero

!------------------------------------------------------------------------------!
!--------------> approximate linearization section <---------------------------!
!------------------------------------------------------------------------------!

  if( imode==1 ) then
    term1 = -eig1 + half*(eig2 + eig3)
    term2 = half*(eig2 - eig3)
    del1  = term1*dpress/abar**2 + term2*rhobar*(uconr - uconl)/abar
    del2  = term1*(uconr - uconl)*rhobar + term2*dpress/abar

    ddel1_dql(:) = - term1*dp_dql(:)/abar**2 - term2*rhobar*ducon_dql(:)/abar
    ddel1_dqr(:) = + term1*dp_dqr(:)/abar**2 + term2*rhobar*ducon_dqr(:)/abar

    ddel2_dql(:) = - term1*ducon_dql(:)*rhobar - term2*dp_dql(:)/abar
    ddel2_dqr(:) = + term1*ducon_dqr(:)*rhobar + term2*dp_dqr(:)/abar

    goto 111
     
  endif

!------------------------------------------------------------------------------!
!-----------> linearization of Roe averaged state <----------------------------!
!------------------------------------------------------------------------------!

  dfact_dql(:) = (half/fact)*(-rhor*rhoinvl**2)*dro_dql(:)
  dfact_dqr(:) = (half/fact)*(rhoinvl)*dro_dqr(:)

  dA_dql(:) = -dfact_dql(:)/(one+fact)**2
  dA_dqr(:) = -dfact_dqr(:)/(one+fact)**2

  dB_dql(:) = dfact_dql(:)/(one + fact)**2
  dB_dqr(:) = dfact_dqr(:)/(one + fact)**2

  drhobar_dql(:) = drho_dql(:)*fact + rhol*dfact_dql(:)
  drhobar_dqr(:) =                    rhol*dfact_dqr(:)

  dubar_dql(:) = du_dql(:)*A + ul*dA_dql(:) + ur*dB_dql(:)
  dubar_dqr(:) = du_dqr(:)*B + ul*dA_dqr(:) + ur*dB_dqr(:)

  dvbar_dql(:) = dv_dql(:)*A + vl*dA_dql(:) + vr*dB_dql(:)
  dvbar_dqr(:) = dv_dqr(:)*B + vl*dA_dqr(:) + vr*dB_dqr(:)

  dwbar_dql(:) = dw_dql(:)*A + wl*dA_dql(:) + wr*dB_dql(:)
  dwbar_dqr(:) = dw_dqr(:)*B + wl*dA_dqr(:) + wr*dB_dqr(:)

  dhbar_dql(:) = dh_dql(:)*A + hl*dA_dql(:) + hr*dB_dql(:)
  dhbar_dqr(:) = dh_dqr(:)*B + hl*dA_dqr(:) + hr*dB_dqr(:)

  dabar_dql(:) = half*gm1/abar*( dhbar_dql(:) - ubar*dubar_dql(:)     &
                                              - vbar*dvbar_dql(:)     &
                                              - wbar*dwbar_dql(:) )

  dabar_dqr(:) = half*gm1/abar*( dhbar_dqr(:) - ubar*dubar_dqr(:)     &
                                              - vbar*dvbar_dqr(:)     &
                                              - wbar*dwbar_dqr(:) )

  duconbar_dql(:) = dubar_dql(:)*nx + dvbar_dql(:)*ny + dwbar_dql(:)*nz
  duconbar_dqr(:) = dubar_dqr(:)*nx + dvbar_dqr(:)*ny + dwbar_dqr(:)*nz

!------------------------------------------------------------------------------!
!------------------> linearization of Eigenvalues <----------------------------!
!------------------------------------------------------------------------------!

  deig1_dql(:) = sign(1.0_dp, uconbar)*duconbar_dql(:)
  deig1_dqr(:) = sign(1.0_dp, uconbar)*duconbar_dqr(:)

  deig2_dql(:) = sign(1.0_dp, uconbar + abar)*(duconbar_dql(:) + dabar_dql(:))
  deig2_dqr(:) = sign(1.0_dp, uconbar + abar)*(duconbar_dqr(:) + dabar_dqr(:))

  deig3_dql(:) = sign(1.0_dp, uconbar - abar)*(duconbar_dql(:) - dabar_dql(:))
  deig3_dqr(:) = sign(1.0_dp, uconbar - abar)*(duconbar_dqr(:) - dabar_dqr(:))

!------------------------------------------------------------------------------!
  term1 = half*(eig2 + eig3) - eig1
  term2 = half*(eig2 - eig3)

  del1  = term1*dpress/abar**2 + term2*rhobar*(uconr - uconl)/abar
  del2  = term1*(uconr - uconl)*rhobar + term2*dpress/abar

  dterm1_dql(:) = half*( deig2_dql(:) + deig3_dql(:) ) - deig1_dql(:) 
  dterm1_dqr(:) = half*( deig2_dqr(:) + deig3_dqr(:) ) - deig1_dql(:) 

  dterm2_dql(:) = half*( deig2_dql(:) - deig3_dql(:) )
  dterm2_dqr(:) = half*( deig2_dqr(:) - deig3_dqr(:) )

  ddel1_dql(:) = dterm1_dql(:)*dpress/abar**2 - term1*dp_dql(:)/abar**2         &
               - two*term1*dpress*dabar_dql(:)/abar**3
  ddel1_dql(:) = ddel1_dql(:) + dterm2_dql(:)*rhobar*( uconr-uconl )/abar       &
               + term2*drhobar_dql(:)*(uconr-uconl)/abar                        &
               - term2*rhobar*ducon_dql(:)/abar                                 &
               - dabar_dql(:)*term2*rhobar*(uconr-uconl)/abar**2

  ddel1_dqr(:) = dterm1_dqr(:)*dpress/abar**2 + term1*dp_dqr(:)/abar**2         &
               - two*term1*dpress*dabar_dqr(:)/abar**3
  ddel1_dqr(:) = ddel1_dqr(:) + dterm2_dqr(:)*rhobar*( uconr-uconl )/abar       &
               + term2*drhobar_dqr(:)*(uconr-uconl)/abar                        &
               + term2*rhobar*ducon_dqr(:)/abar                                 &
               - dabar_dqr(:)*term2*rhobar*(uconr-uconl)/abar**2

  ddel2_dql(:) = dterm1_dql(:)*(uconr-uconl)*rhobar                             &
               - term1*ducon_dql(:)*rhobar + term1*(uconr-uconl)*drhobar_dql(:)
  ddel2_dql(:) = ddel2_dql(:) + dterm2_dql(:)*dpress/abar                       &
               - term2*dp_dql(:)/abar - dabar_dql(:)*term2*dpress/abar**2

  ddel2_dqr(:) = dterm1_dqr(:)*(uconr-uconl)*rhobar                             &
               + term1*ducon_dqr(:)*rhobar + term1*(uconr-uconl)*drhobar_dqr(:)
  ddel2_dqr(:) = ddel2_dqr(:) + dterm2_dqr(:)*dpress/abar                       &
               + term2*dp_dqr(:)/abar - dabar_dqr(:)*term2*dpress/abar**2

!------------------------------------------------------------------------------!
111 continue
!------------------------------------------------------------------------------!
!-----------------------> Roe flux Jacobian <----------------------------------!
!------------------------------------------------------------------------------!


  !------------> common linearization terms

  lmat(:,:) = 0.0
  rmat(:,:) = 0.0

! create identity matrix     
  imat(:,:) = 0.0
  do k = 1, 5
    imat(k,k) = 1.0
  enddo

  lmat(:,:) = lmat(:,:) - eig1*imat(:,:)
  lmat(1,:) = lmat(1,:) + ddel1_dql(:)
  lmat(2,:) = lmat(2,:) + ddel1_dql(:)*ubar + ddel2_dql(:)*nx
  lmat(3,:) = lmat(3,:) + ddel1_dql(:)*vbar + ddel2_dql(:)*ny
  lmat(4,:) = lmat(4,:) + ddel1_dql(:)*wbar + ddel2_dql(:)*nz
  lmat(5,:) = lmat(5,:) + ddel1_dql(:)*hbar + ddel2_dql(:)*uconbar

  rmat(:,:) = rmat(:,:) + eig1*imat(:,:) 
  rmat(1,:) = rmat(1,:) + ddel1_dqr(:)
  rmat(2,:) = rmat(2,:) + ddel1_dqr(:)*ubar + ddel2_dqr(:)*nx
  rmat(3,:) = rmat(3,:) + ddel1_dqr(:)*vbar + ddel2_dqr(:)*ny
  rmat(4,:) = rmat(4,:) + ddel1_dqr(:)*wbar + ddel2_dqr(:)*nz
  rmat(5,:) = rmat(5,:) + ddel1_dqr(:)*hbar + ddel2_dqr(:)*uconbar

! additional terms for exact linearization

  if(imode/=1) then
    do j = 1, 5
      do i = 1, 5
        lmat(i,j) = lmat(i,j) + ( qr(i) - ql(i) )* deig1_dql(j)
        rmat(i,j) = rmat(i,j) + ( qr(i) - ql(i) )* deig1_dqr(j)
      end do
    end do

    lmat(2,:) = lmat(2,:) + del1*dubar_dql(:)
    rmat(2,:) = rmat(2,:) + del1*dubar_dqr(:)

    lmat(3,:) = lmat(3,:) + del1*dvbar_dql(:)
    rmat(3,:) = rmat(3,:) + del1*dvbar_dqr(:)

    lmat(4,:) = lmat(4,:) + del1*dwbar_dql(:)
    rmat(4,:) = rmat(4,:) + del1*dwbar_dqr(:)

    lmat(5,:) = lmat(5,:) + del1*dhbar_dql(:) + del2*duconbar_dql(:)
    rmat(5,:) = rmat(5,:) + del1*dhbar_dqr(:) + del2*duconbar_dqr(:)
  endif

! left state

  lmat1(1,:) = drho_dql(:)*uconl + rhol*ducon_dql(:)

  lmat1(2,:) = drho_dql(:)*uconl*ul + rhol*ducon_dql(:)*ul + rhol*uconl*du_dql(:)
  lmat1(3,:) = drho_dql(:)*uconl*vl + rhol*ducon_dql(:)*vl + rhol*uconl*dv_dql(:)
  lmat1(4,:) = drho_dql(:)*uconl*wl + rhol*ducon_dql(:)*wl + rhol*uconl*dw_dql(:)

  lmat1(3,:) = lmat1(3,:) + ny*dp_dql(:)
  lmat1(2,:) = lmat1(2,:) + nx*dp_dql(:)
  lmat1(4,:) = lmat1(4,:) + nz*dp_dql(:)

  lmat1(5,:) = ( dq5_dql(:) + dp_dql(:) )*uconl + (ql(5) + pl)*ducon_dql(:)

! right state

  rmat1(1,:) = drho_dqr(:)*uconr    + rhor*ducon_dqr(:)

  rmat1(2,:) = drho_dqr(:)*uconr*ur + rhor*ducon_dqr(:)*ur + rhor*uconr*du_dqr(:)
  rmat1(3,:) = drho_dqr(:)*uconr*vr + rhor*ducon_dqr(:)*vr + rhor*uconr*dv_dqr(:)
  rmat1(4,:) = drho_dqr(:)*uconr*wr + rhor*ducon_dqr(:)*wr + rhor*uconr*dw_dqr(:)

  rmat1(3,:) = rmat1(3,:) + ny*dp_dqr(:)
  rmat1(2,:) = rmat1(2,:) + nx*dp_dqr(:)
  rmat1(4,:) = rmat1(4,:) + nz*dp_dqr(:)

  rmat1(5,:) = ( dq5_dqr(:) + dp_dqr(:) )*uconr + (qr(5) + pr)*ducon_dqr(:)

! finalize left and right jacobians

  lmat(:,:) = half*( lmat1(:,:) - lmat(:,:) )
  rmat(:,:) = half*( rmat1(:,:) - rmat(:,:) )

end subroutine roe_jacobian
