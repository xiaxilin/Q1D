module solvers

  use set_precision, only : dp

  implicit none

  private

  public :: explicit_solve
  public :: implicit_solve
  public :: set_time_step
  public :: check_convergence

  public :: iterations ! total number of iterations
  public :: itercheck  ! check for convergence every itercheck iters
  public :: iter_out   ! write solution every iter_out iterations
  public :: rkorder    ! multistep RK order 1, 2, 3, or 4
  public :: cfl        ! CFL limit at start
  public :: cfl_end    ! CFL limit at end if ramping
  public :: cfl_ramp   ! Number of iterations to ramp CFL over

  public :: toler      ! convergence tolerance

  public :: solver     ! 'explicit' or 'implicit' solver

  integer  :: iterations
  integer  :: itercheck
  integer  :: iter_out
  integer  :: rkorder
  integer  :: cfl_ramp

  real(dp) :: cfl
  real(dp) :: cfl_end

  real(dp) :: toler

  character(8)  :: solver

  contains

!=============================== explicit_solve ==============================80
!
! Routine to handle explicit solve
!
!=============================================================================80
  subroutine explicit_solve( cells, faces, prim_cc, cons_cc, cell_vol,         &
                             area_f, dx, dadx_cc, x_cc )

    use set_precision, only : dp
    use set_constants, only : zero
    use residual,      only : create_residual
    use bc,            only : subsonic_inflow_explicit, outflow_explicit
    use write_soln,    only : write_soln_line

    implicit none

    integer,                        intent(in)    :: cells, faces
    real(dp), dimension(3,cells+2), intent(inout) :: prim_cc, cons_cc
    real(dp), dimension(cells+2),   intent(in)    :: cell_vol
    real(dp), dimension(faces),     intent(in)    :: area_f
    real(dp), dimension(cells+2),   intent(in)    :: dx, dadx_cc, x_cc

    integer  :: n, rk, cell

    real(dp), dimension(3)         :: cons
    real(dp), dimension(4)         :: rk_const
    real(dp), dimension(cells+2)   :: dt
    real(dp), dimension(3,cells+2) :: resid

    logical :: convergence_flag = .false.

    continue

    rk_const = zero
    do rk = 1, rkorder
      rk_const(rk) = real(1+rkorder-rk,dp)
    end do

    iteration_loop : do n = 0, iterations
! Set both local and global time step
      dt = set_time_step( cells, dx, prim_cc )

! Store solution before RK loop. Wouldn't be necessary for pure Euler explicit
      do cell = 1, cells+2
        cons_cc(:,cell) = primitive_to_conserved_1D(prim_cc(:,cell))
      end do

      rk_loop : do rk = 1, rkorder
        call create_residual( cells, faces, n, prim_cc, area_f, dadx_cc, dx,   &
                              resid )

! Perform explicit iterations on interior cells
! Note that the conserved variables are converted to primitive and floored
        do cell = 2,cells+1
          cons(:) = cons_cc(:,cell) - dt(cell)*resid(:,cell)                   &
                  / ( rk_const(rkorder)*cell_vol(cell) )
          prim_cc(:,cell) = conserved_to_primitive_1D(cons(:))
          prim_cc(:,cell) = floor_primitive_vars(prim_cc(:,cell))
        end do

! Enforce BC's
!        select case(inflow)
!        case('sub')
        call subsonic_inflow_explicit(prim_cc(:,1), prim_cc(:,2), prim_cc(:,3))
!          call set_sub_sonic_inflow_r(prim_cc(:,1), prim_cc(:,2) )
!        case('ss')
          ! don't need to do anything... they're already set
!        end select

! This routine handles both sub and supersonic outflow
        call outflow_explicit( prim_cc(:,cells+2),                             &
                               prim_cc(:,cells+1), prim_cc(:,cells) )

      end do rk_loop

      do cell = 1, cells+2
        cons_cc(:,cell) = primitive_to_conserved_1D(prim_cc(:,cell))
      end do

      if (iter_out >= 0 .and. mod(n,iter_out) == 0) then
        call write_soln_line(n, cells, x_cc, prim_cc, cons_cc)
      end if

      if (mod(n,itercheck) == 0) then
        call check_convergence(cells, n, resid, convergence_flag)

        if ( convergence_flag ) then
          write(*,*) 'Solution has converged!'
          exit
        end if
      end if

!      if (iter_restar >= 0 .and. mod(n,iter_restart) == 0) then
!        call write_restart(cells, prim_cc)
!      end if

    end do iteration_loop

    if (.not. convergence_flag ) then
      write(*,*) 'Solution failed to converge...'
      write(*,*) 'Consider continuing from q1d.rst'
    end if

    call write_soln_line(iterations, cells, x_cc, prim_cc, cons_cc)

  end subroutine explicit_solve

!=============================== implicit_solve ==============================80
!
! Routine to handle implicit solve
!
!=============================================================================80
  subroutine implicit_solve( cells, faces, prim_cc, cons_cc, cell_vol,         &
                             area_f, dx, dadx_cc, x_cc )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two
    use fluid_constants, only : gm1
    use residual,        only : create_residual, firstorder
    use lhs,             only : fill_lhs, fill_full_lhs, lhs_order
    use bc,              only : subsonic_inflow, set_outflow, set_outflow_prim
    use matrix_manip,    only : triblocksolve, pentablocksolve
    use write_soln,      only : write_soln_line

    implicit none

    integer,                        intent(in)    :: cells, faces
    real(dp), dimension(3,cells+2), intent(inout) :: prim_cc, cons_cc
    real(dp), dimension(cells+2),   intent(in)    :: cell_vol
    real(dp), dimension(faces),     intent(in)    :: area_f
    real(dp), dimension(cells+2),   intent(in)    :: dx, dadx_cc, x_cc

    integer  :: n, eq, cell

    real(dp), dimension(3)           :: l2out
    real(dp), dimension(cells+2)     :: dt
    real(dp), dimension(3,cells+2)   :: RHS, delta_cons_cc!, delta_prim_cc
    real(dp), dimension(3,3,cells+2) :: L2, L, D, U, U2

    logical :: convergence_flag = .false.

    continue

    main_loop : do n = 0, iterations

! Time step calculations
      if (n <= cfl_ramp) then
        cfl = cfl + cfl_end/real(cfl_ramp,dp)
      end if

      dt = set_time_step( cells, dx, prim_cc )

! Form RHS
      call create_residual( cells, faces, n, prim_cc, area_f, dadx_cc, dx, RHS )
! Account for sign since the residual is moved to the RHS
      RHS = -RHS

! Check residuals for convergence before they are destroyed by the direct solver
      if (mod(n,itercheck) == 0) then
        call check_convergence(cells, n, RHS, convergence_flag, l2out)
        if ( convergence_flag ) then
          write(*,*) 'Solution has converged!'
          exit
        end if
      end if

! Form LHS
      if ( n < firstorder .or. lhs_order == 1 ) then
        call fill_lhs( cells, cell_vol, area_f, dadx_cc, dt, cons_cc, L, D, U )
!        call fill_lhs( cells, cell_vol, area_f, dadx_cc, dt, prim_cc, L, D, U )
      else
        call fill_full_lhs( cells, cell_vol, area_f, dadx_cc, &
                            dt, cons_cc, L2, L, D, U, U2 )
!        call fill_full_lhs( cells, cell_vol, area_f, dadx_cc, &
!                            dt, prim_cc, L2, L, D, U, U2 )
      end if

! Take care of BC's
! Inflow, modify according to bc
      call subsonic_inflow( n, cons_cc(:,1), cons_cc(:,2), cons_cc(:,3),       &
                            D(:,:,1), U(:,:,1), U2(:,:,1), RHS(:,1) )
!      call subsonic_inflow( n, prim_cc(:,1), prim_cc(:,2), prim_cc(:,3),      &
!                            D(:,:,1), U(:,:,1), U2(:,:,1), RHS(:,1) )

! Outflow
      call set_outflow( n, cons_cc(:,cells+2), cons_cc(:,cells+1),             &
                        cons_cc(:,cells),                                      &
                        D(:,:,cells+2), L(:,:,cells+2), L2(:,:,cells+2),       &
                        RHS(:,cells+2) )
!      call set_outflow_prim( n, prim_cc(:,cells+2), prim_cc(:,cells+1),       &
!                        prim_cc(:,cells),                                     &
!                        D(:,:,cells+2), L(:,:,cells+2), L2(:,:,cells+2),      &
!                        RHS(:,cells+2) )

! Solve the system of equations
      if ( n < firstorder .or. lhs_order == 1 ) then
        call triblocksolve(3, cells+2, L, D, U, RHS, delta_cons_cc)
!        call triblocksolve(3, cells+2, L, D, U, RHS, delta_prim_cc)
      else
        call pentablocksolve(3, cells+2, L2, L, D, U, U2, RHS, delta_cons_cc)
!        call pentablocksolve(3, cells+2, L2, L, D, U, U2, RHS, delta_prim_cc)
      end if

! Update the conserved variables
      cons_cc = cons_cc + delta_cons_cc
!      prim_cc = prim_cc + delta_prim_cc

! Floor variables for stability
      do cell = 1, cells+2
        prim_cc(:,cell) = conserved_to_primitive_1D(cons_cc(:,cell))
        if ( prim_cc(1,cell) <= zero .or. prim_cc(3,cell) <= zero ) then
          prim_cc(1,cell) = 0.01_dp
          prim_cc(2,cell) = 1.0_dp
          prim_cc(3,cell) = 3000.0_dp
        end if
        cons_cc(:,cell) = primitive_to_conserved_1D(prim_cc(:,cell))
      end do

      if (iter_out >= 0 .and. mod(n,iter_out) == 0) then
        call write_soln_line(n, cells, x_cc, prim_cc, cons_cc)
      end if

!      if (iter_restar >= 0 .and. mod(n,iter_restart) == 0) then
!        call write_restart(cells, prim_cc)
!      end if

    end do main_loop

    if (.not. convergence_flag ) then
      write(*,*) 'Solution failed to converge...'
      write(*,*) 'Consider continuing from q1d.rst'
    end if

    call write_soln_line(iterations, cells, x_cc, prim_cc, cons_cc)

  end subroutine implicit_solve

!=============================== set_time_step ===============================80
!
! Sets time step to be either local or global
!
!=============================================================================80
  pure function set_time_step( cells, dx, prim_cc ) result(dt)

    use set_precision, only : dp
    use set_constants, only : large

    implicit none

    integer,                         intent(in)  :: cells
    real(dp), dimension(cells+2),    intent(in)  :: dx
    real(dp), dimension(3, cells+2), intent(in)  :: prim_cc
    real(dp), dimension(cells+2)                 :: dt

    integer  :: cell
    real(dp) :: a, dt_global

    continue

    dt_global = large

    do cell = 1, cells+2
      a = speed_of_sound( prim_cc(3,cell), prim_cc(1,cell) )
      dt(cell) = cfl*dx(cell) / ( abs(prim_cc(2,cell)) + a )
      dt_global = min(dt_global, dt(cell))
    end do

! FIXME: fix this mess!
!    if ( solver == 'explicit' .or. time_accurate) then
!    if ( solver == 'explicit' ) then
    dt(:) = dt_global
!    end if

  end function set_time_step

!============================= check_convergence =============================80
!
! Checks for convergence based on residuals normalized at the first iteration
!
!=============================================================================80
  subroutine check_convergence( cells, iteration, residuals, convergence_flag, &
                                l2out)

    use set_precision,   only : dp
    use set_constants,   only : zero
    use initialize_soln, only : restart, L1_init, L2_init, Linf_init

    implicit none

    integer,                          intent(in)  :: cells, iteration
    real(dp), dimension(3,cells+2),   intent(in)  :: residuals
    logical,                          intent(out) :: convergence_flag
    real(dp), dimension(3), optional, intent(out) :: l2out

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

    if ( present(l2out) ) l2out = l2

  end subroutine check_convergence

! begin include statements

  include 'conserved_to_primitive_1D.f90'
  include 'floor_primitive_vars.f90'
  include 'primitive_to_conserved_1D.f90'
  include 'speed_of_sound.f90'

end module solvers
