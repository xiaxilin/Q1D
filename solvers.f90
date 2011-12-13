module solvers

  use set_precision, only : dp

  implicit none

  private

  public :: explicit_solve

  public :: iterations ! total number of iterations
  public :: firstorder ! if 2nd order, how many 1st order iters for stability?
  public :: itercheck  ! check for convergence every itercheck iters
  public :: iter_out   ! write solution every iter_out iterations
  public :: rkorder    ! multistep RK order 1, 2, 3, or 4
  public :: cfl        ! CFL limit
  public :: muscl      ! MUSCL extrapolation, .true. or .false.
  public :: kappa      ! form of MUSCL
  public :: limiter    ! type of variable limiting 
  public :: toler      ! convergence tolerance

  public :: flux_type  ! 'central', 'jst', 'sw', 'vanleer', 'roe', 'ausm'
  public :: k2         ! JST damping coefficient, only for flux_type = 'jst'
  public :: k4         ! JST damping coefficient, only for flux_type = 'jst'

! set initial values

  integer  :: iterations
  integer  :: firstorder
  integer  :: itercheck
  integer  :: iter_out
  integer  :: rkorder
  logical  :: muscl
  real(dp) :: cfl
  real(dp) :: kappa
  real(dp) :: toler
  character(len=10) :: limiter

  character(len=10) :: flux_type
  real(dp) :: k2
  real(dp) :: k4

  contains

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine explicit_solve( cells, faces, dxsi, prim_cc, cons_cc,             &
                             area_f, area_cc, dxdxsi_cc, dadx_cc, x_cc )

    use set_precision, only : dp
    use write_soln,    only : write_soln_line

    implicit none

    integer,                        intent(in)    :: cells, faces
    real(dp),                       intent(in)    :: dxsi
    real(dp), dimension(3,cells+2), intent(inout) :: prim_cc, cons_cc
    real(dp), dimension(faces),     intent(in)    :: area_f
    real(dp), dimension(cells+2),   intent(in)    :: area_cc, dxdxsi_cc, dadx_cc
    real(dp), dimension(cells+2),   intent(in)    :: x_cc

    integer  :: n, rk, cell, eq
    real(dp) :: dt_global
    real(dp), dimension(3,cells+2) :: cons_cc_0, residual
    real(dp), dimension(cells+2)   :: dt

    logical :: convergence_flag = .false.

    continue

    do n = 0, iterations
! set both local and global time step
      call set_time_step( cells, dxsi, prim_cc, dt, dt_global )

      dt(:) = dt_global

! make copy of solution for RK schemes... 
! wouldn't be necessary for pure Euler explicit
      cons_cc_0 = cons_cc

      do rk = 1, rkorder
        call create_residual( cells, faces, n, dxsi, prim_cc, cons_cc,         &
                              area_f, dadx_cc, dxdxsi_cc, residual )

! perform explicit iterations on interior cells

        do cell = 2,cells+1
          do eq = 1,3
            cons_cc(eq,cell) = cons_cc_0(eq,cell)                              &
                            + (dt(cell)/area_cc(cell)) * residual(eq,cell)     &
                            / real(1+rkorder-rk,dp)
          end do
          prim_cc(:,cell) = conserved_to_primitive_1D(cons_cc(:,cell))
          prim_cc(:,cell) = floor_primitive_vars(prim_cc(:,cell))
        end do

!        select case(inflow)
!        case('sub')
          call set_sub_sonic_inflow(prim_cc(:,1), prim_cc(:,2), prim_cc(:,3))
!          call set_sub_sonic_inflow_r(prim_cc(:,1), prim_cc(:,2) )
!        case('ss')
          ! don't need to do anything... they're already set
!        end select

! This routine handles both sub and supersonic conditions
          call set_outflow( prim_cc(:,cells+2),                                &
                            prim_cc(:,cells+1), prim_cc(:,cells) )

        do cell = 1, cells+2
          cons_cc(:,cell) = primitive_to_conserved_1D(prim_cc(:,cell))
        end do
      end do

      if (mod(n,itercheck) == 0) then
        call check_convergence(cells, n, residual, convergence_flag)
      end if

      if (iter_out >= 0 .and. mod(n,iter_out) == 0) then
        call write_soln_line(n, cells, x_cc, prim_cc, cons_cc)
      end if

!      if (iter_restar >= 0 .and. mod(n,iter_restart) == 0) then
!        call write_restart(cells, prim_cc)
!      end if

      if ( convergence_flag ) then
        write(*,*) 'Solution has converged!'
        exit
      end if
    end do

    if (.not. convergence_flag ) then
      write(*,*) 'Solution failed to converge...'
      write(*,*) 'Consider continuing from q1d.rst'
    end if

    call write_soln_line(iterations, cells, x_cc, prim_cc, cons_cc)

  end subroutine explicit_solve

  include 'conserved_to_primitive_1D.f90'
  include 'floor_primitive_vars.f90'
  include 'primitive_to_conserved_1D.f90'

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine set_time_step( cells, dxsi, prim_cc, dt, dt_global)

    use set_precision, only : dp
    use set_constants, only : large
    
    implicit none

    integer,                         intent(in)  :: cells
    real(dp),                        intent(in)  :: dxsi
    real(dp), dimension(3, cells+2), intent(in)  :: prim_cc
    real(dp), dimension(cells+2),    intent(out) :: dt
    real(dp),                        intent(out) :: dt_global

    integer  :: cell
    real(dp) :: a

    continue

    dt_global = large

    do cell = 1, cells+2
      a = speed_of_sound( prim_cc(3,cell), prim_cc(1,cell) )
      dt(cell) = cfl*dxsi / ( abs(prim_cc(2,cell)) + a )
      dt_global = min(dt_global, dt(cell))
    end do

  end subroutine set_time_step

  include 'speed_of_sound.f90'

!=============================================================================80
!
! 
!
!=============================================================================80

  subroutine create_residual( cells, faces, iteration, dxsi, prim_cc, cons_cc, &
                              area_f, dadx_cc, dxdxsi_cc, residual )

    use set_precision, only : dp
    use set_constants, only : zero

    implicit none

    integer,                        intent(in)  :: cells, faces, iteration
    real(dp),                       intent(in)  :: dxsi
    real(dp), dimension(3,cells+2), intent(in)  :: prim_cc
    real(dp), dimension(3,cells+2), intent(in)  :: cons_cc
    real(dp), dimension(faces),     intent(in)  :: area_f
    real(dp), dimension(cells+2),   intent(in)  :: dadx_cc
    real(dp), dimension(cells+2),   intent(in)  :: dxdxsi_cc
    real(dp), dimension(3,cells+2), intent(out) :: residual

    integer :: cell, eq
    real(dp), dimension(3,faces)   :: F
    real(dp), dimension(3,cells+2) :: S

    continue

    call create_fluxes( cells, faces, iteration, prim_cc, cons_cc, F )
    call create_source( cells, prim_cc(3,:), dadx_cc, S )

    residual = zero

!$OMP parallel 
  !$OMP do
    do cell = 2, cells+1
      do eq = 1,3
        residual(eq,cell) = S(eq,cell)                                         &
                          - (area_f(cell)   * F(eq,cell)                       &
                          -  area_f(cell-1) * F(eq,cell-1))                    &
                          / (dxdxsi_cc(cell)*dxsi)
      end do
    end do
  !$OMP end do
!$OMP end parallel

  end subroutine create_residual

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine check_convergence(cells, iteration, residuals, convergence_flag)

    use set_precision,   only : dp
    use set_constants,   only : zero
    use initialize_soln, only : restart, L1_init, L2_init, Linf_init

    implicit none

    integer,                        intent(in)  :: cells, iteration
    real(dp), dimension(3,cells+2), intent(in)  :: residuals
    logical,                        intent(out) :: convergence_flag

    integer :: cell
    real(dp), dimension(3) :: L1, L2, Linf

    continue

    L1 = zero
    L2 = zero
    Linf = zero

    do cell = 2, cells+1
      L1(:) = L1(:) + abs(residuals(:,cell))
      L2(:) = L2(:) + residuals(:,cell)**2
      Linf(:) = max(Linf(:), abs(residuals(:,cell)))
    end do

    L1(:) = L1(:)/real(cells,dp)
    L2(:) = sqrt(L2(:))/real(cells,dp)

    if (iteration == 0 .and. .not. restart) then
      L1_init(:)   = L1(:)
      L2_init(:)   = L2(:)
      Linf_init(:) = Linf(:)
    end if

    L1(:)   = L1(:)/L1_init(:)
    L2(:)   = L2(:)/L2_init(:)
    Linf(:) = Linf(:)/Linf_init(:)

    write(*,300) iteration, L2(1), L2(2), L2(3)
300 format(1X,i8,2(e15.6),3(e15.6),4(e15.6))

    convergence_flag = .false.
    if (L2(1) <= toler .and. L2(2) <= toler .and. L2(3) <= toler ) then
      convergence_flag = .true.
    end if

  end subroutine check_convergence

!=============================================================================80
!
! Uses primitive variables
!
!=============================================================================80

  subroutine set_sub_sonic_inflow( cc_in, cc_1, cc_2 )

    use set_precision,   only : dp
    use set_constants,   only : one, half, two
    use fluid_constants, only : r, gamma, gm1, xgm1, gxgm1,gm1xgp1, gp1
    use initialize_soln, only : po, to

    implicit none

    real(dp), dimension(3), intent(in)  :: cc_1, cc_2
    real(dp), dimension(3), intent(out) :: cc_in

    real(dp) :: vel_max, psi

! testing
!    real(dp) :: rho0, h0, s0, p, u

    continue

! set max, physically possible velocity
    vel_max = sqrt(two*gamma*r*to*xgm1)-one

! extrapolate velocity and limit
    cc_in(2) = max(-vel_max, min(two*cc_1(2)-cc_2(2), vel_max))

! now calculate inflow
    psi = to/(to-(gm1*cc_in(2)**2/(two*gamma*r)))

    cc_in(1) = po/(r*to*psi**xgm1)
    cc_in(3) = po/psi**gxgm1

! to test another version...

!    rho0 = po/(r*to)
!    h0 = gxgm1*po/rho0
!    s0 = r*xgm1*log(po/rho0**gamma)

!    cc_in(1) = min(max(two*cc_1(1) - cc_2(1),0.0001_dp),rho0)

!    p = exp(s0*gm1/r)*cc_in(1)**gamma
!    u = sqrt(two*(h0 - gxgm1*p/cc_in(1)))

!    cc_in(2) = u
!    cc_in(3) = p

! floor variables
    cc_in = floor_primitive_vars(cc_in)

  end subroutine set_sub_sonic_inflow
!=============================================================================80
!
! Uses primitive variables
!
!=============================================================================80

  subroutine set_sub_sonic_inflow_r( cc_in, cc_1 )

    use set_precision,   only : dp
    use set_constants,   only : one, half, two
    use fluid_constants, only : r, gamma, gm1, xgm1, gxgm1,gm1xgp1, gp1
    use initialize_soln, only : po, to

    implicit none

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

  end subroutine set_sub_sonic_inflow_r

!=============================================================================80
!
! 
!
!=============================================================================80

  subroutine set_outflow( cc_out, cc_1, cc_2 )

    use set_precision,   only : dp
    use set_constants,   only : two
    use initialize_soln, only : pback

    implicit none

    real(dp), dimension(3), intent(in)  :: cc_1, cc_2
    real(dp), dimension(3), intent(out) :: cc_out

    integer :: eq

    continue

! extrapolate all variables
    do eq = 1, 3
      cc_out(eq) = two*cc_1(eq) - cc_2(eq)
    end do

! set back pressure if appropriate
    if ( pback >= 0.0_dp ) cc_out(3) = pback

! floor variables
    cc_out = floor_primitive_vars(cc_out)

  end subroutine set_outflow


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

!$OMP parallel do
      do i = 1, faces
        prim_left(:,i)  = floor_primitive_vars(prim_left(:,i))
        prim_right(:,i) = floor_primitive_vars(prim_right(:,i))
      end do
!$OMP end parallel do

    end if

! create flux vectors
! could be made $OMP for each of these cases
    select case(trim(flux_type))
    case ('central')
!$OMP parallel do
      do i = 1, faces
        flux(:,i) = flux_central( prim_cc(:,i), prim_cc(:,i+1) )
      end do
!$OMP end parallel do
    case ('jst')
!$OMP parallel do
      do i = 1, faces
        flux(:,i) = flux_central( prim_cc(:,i), prim_cc(:,i+1) )
      end do
!$OMP end parallel do
      call add_jst_damping( cells, faces, prim_cc, cons_cc, flux )
    case('vanleer')
!$OMP parallel do
      do i = 1, faces
        flux(:,i) = flux_vanleer( prim_left(:,i), prim_right(:,i) )
      end do
!$OMP end parallel do
    case('sw')
!$OMP parallel do
      do i = 1, faces
        flux(:,i) = flux_sw( prim_left(:,i), prim_right(:,i) )
      end do
!$OMP end parallel do
    case('ausm')
!$OMP parallel do
      do i = 1, faces
        flux(:,i) = flux_ausm( prim_left(:,i), prim_right(:,i) )
      end do
!$OMP end parallel do
    case('roe')
!$OMP parallel do
      do i = 1, faces
        flux(:,i) = flux_roe( prim_left(:,i), prim_right(:,i) )
      end do
!$OMP end parallel do
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
                                  vars_left, vars_right )

    use set_precision, only : dp
    use set_constants, only : zero, fourth, one, small

    implicit none

    integer,                        intent(in)  :: cells, faces, iteration
    real(dp), dimension(3,cells+2), intent(in)  :: vars_cc
    real(dp), dimension(3,faces),   intent(out) :: vars_left, vars_right

    integer                       :: i
    real(dp), dimension(3, faces) :: r_L, r_R, psi_L, psi_R

    real(dp), parameter :: small_factor = 0.0000001_dp

    continue

    if ( iteration <= firstorder .or. .not. muscl ) then
! skip excess computations if only first order
!OMP parallel do
      do i = 1, faces
        vars_left(:,i)  = vars_cc(:,i)
        vars_right(:,i) = vars_cc(:,i+1)
      end do
!OMP end parallel do
    else

! calculate left side variations, r>=0
!OMP parallel do
      do i = 1, faces-1
        r_L(:,i)  = max( zero, (vars_cc(:,i+2) - vars_cc (:,i+1))              &
                             / (vars_cc(:,i+1) - vars_cc(:,i) + small_factor) )
      end do
!OMP end parallel do
      r_L(:,faces) = one

! calculate right side variations, r>=0

      r_R(:,1)     = one
!OMP parallel do
      do i = 2, faces
        r_R(:,i) = max( zero, (vars_cc(:,i) - vars_cc(:,i-1))                  &
                            / (vars_cc(:,i+1) - vars_cc(:,i) + small_factor) )
      end do
!OMP end parallel do

! apply appropriate limiter
!OMP for each of these
      select case(trim(limiter))
      case('minmod')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(3, r_L(:,i), 1.0_dp)
          psi_R(:,i) = limiter_sweby(3, r_R(:,i), 1.0_dp)        
        end do
      case('superbee')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(3, r_L(:,i), 2.0_dp)
          psi_R(:,i) = limiter_sweby(3, r_R(:,i), 2.0_dp)  
        end do
      case('sweby')
        do i = 1, faces
          psi_L(:,i) = limiter_sweby(3, r_L(:,i), 1.5_dp)
          psi_R(:,i) = limiter_sweby(3, r_R(:,i), 1.5_dp)
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

! perform MUSCL extrapolation
      i = 1
      vars_left(:,i) = vars_cc(:,i) + fourth                                   &
                * ((one+kappa)*(vars_cc(:,i+1) - vars_cc(:,i)))

!OMP parallel do
      do i = 2, faces
        vars_left(:,i)    = vars_cc(:,i) + fourth                              &
                * ((one+kappa)*psi_R(:,i)   * (vars_cc(:,i+1) - vars_cc(:,i))  &
                +  (one-kappa)*psi_L(:,i-1) * (vars_cc(:,i)   - vars_cc(:,i-1)))

        vars_right(:,i-1) = vars_cc(:,i) - fourth                              &
                * ((one-kappa)*psi_R(:,i)   * (vars_cc(:,i+1) - vars_cc(:,i))  &
                +  (one+kappa)*psi_L(:,i-1) * (vars_cc(:,i)   - vars_cc(:,i-1)))
      end do
!OMP end parallel do

      i = faces
      vars_right(:,i) = vars_cc(:,i) - fourth                                  &
                * ((one+kappa)*(vars_cc(:,i) - vars_cc(:,i-1)))

    end if

  end subroutine muscl_extrapolation
 
! begin include statements
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

end module solvers
