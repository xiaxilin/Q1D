! The Q1D-nozzle has an exact solution from the mach/area relations
! This module will hold the routines necessary to calculate the exact soln,
! form the DE error norms, and calculate the TE.

! FIXME: change subsonic/supersonic flag so that it is calculated from the
!        computed solution as well as mach_init
! FIXME: find a_star by searching area_f so that this is more general

module solution_error

  implicit none

  private

  public :: calculate_exact_soln
  public :: estimate_te

contains

!============================ calculate_exact_soln ===========================80
!
! Calculates the exact solution for a converging-diverging quasi-1D nozzle
!
!=============================================================================80
  subroutine calculate_exact_soln( cells, x_cc, area_cc, a_star, a_e, cons_cc, &
                                   cons_cc_ex)

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one, two
    use fluid_constants, only : gamma, gm1, xgm1, xgp1, gxgm1, gp1xgm1, r
    use initialize_soln, only : to, po, pback

    integer,                                    intent(in)  :: cells
    real(dp), dimension(0:cells+1),             intent(in)  :: x_cc, area_cc
    real(dp),                                   intent(in)  :: a_star, a_e
    real(dp), dimension(3,cells+2),             intent(in)  :: cons_cc
    real(dp), optional, dimension(3,0:cells+1), intent(out) :: cons_cc_ex

    integer               :: i, i_throat, i_shock, unit
    integer, dimension(1) :: i_min      ! needed for minloc function

    real(dp) :: psi       ! ( = T_0/T )
    real(dp) :: temp      ! (Temperature, T)
    real(dp) :: mach_init ! Initial Mach number
    real(dp) :: ratio, M_e, M_1, po_e, a_star_new
    real(dp) :: derho, derhou, derhoet
    real(dp), dimension(3)         :: cons_exact
    real(dp), dimension(0:cells+1)   :: mach_exact
    real(dp), dimension(3,0:cells+1) :: soln_exact

    continue

    i_min = minloc( area_cc(1:cells) )
    i_throat = i_min(1)
    if ( i_throat == 0 .or. i_throat == 1 ) i_throat = -1

! calculate exact mach/area solution
    mach_init = 0.1_dp
! subsonic section to throat
    do i = 1, i_throat
      mach_exact(i) = mach_from_area(area_cc(i)/a_star, mach_init, 0)
      mach_init = mach_exact(i)
    end do
! sub/supersonic section after throat
    mach_init = mach_init + 0.01
    do i = i_throat+1, cells
      mach_exact(i) = mach_from_area(area_cc(i)/a_star, mach_init, 1)
      mach_init = mach_exact(i)
    end do

! once the isentropic solution has been found, check for shocked case
    nonisentropic : if (pback > zero) then
      ratio = pback*a_e/(po*a_star)

      M_e = subsonic_mach_at_exit(ratio)

      po_e = pback*(one + half*gm1*M_e**2)**(gxgm1)

      a_star_new = a_e*M_e*(two*xgp1*(one+half*gm1*M_e**2))**(-half*gp1xgm1)

      ratio = po_e/po

      M_1 = pre_shock_mach(ratio)

      shock_location : do i = cells, i_throat+1, -1
        if (mach_exact(i) <= M_1) then
          i_shock = i+1 ! Go upstream one cell to be post shock
          exit
        end if
      end do shock_location

      mach_init = (one + half*gm1*M_1**2)/(gamma*M_1**2-half*gm1)

      do i = i_shock, cells
        mach_exact(i) = mach_from_area(area_cc(i)/a_star_new, mach_init, 0)
        mach_init = mach_exact(i)
      end do

    end if nonisentropic

! once the mach/area solution has been found, calculate the primitive variables
    do i = 1, cells
      psi = one + half*gm1*mach_exact(i)**2
      temp = to / psi
      soln_exact(3,i) = po / ( psi**gxgm1 )
      if ( pback > zero ) then
        if ( i >= i_shock ) then
          soln_exact(3,i) = po_e / ( psi**gxgm1 )
        end if
      end if
      soln_exact(1,i) = soln_exact(3,i)/( r*temp )
      soln_exact(2,i) = mach_exact(i) &
                      * speed_of_sound(soln_exact(3,i), soln_exact(1,i))
    end do

    if ( .not. present(cons_cc_ex) ) then
      unit = find_available_unit()

! set up output file for exact solution
      open(unit,file='q1d_exact_soln.dat',status='replace')
      write(unit,*) 'TITLE = "Quasi-1D Nozzle: Exact Solution"'
      write(unit,*) 'variables="x(m)""Area(m^2)""rho(kg/m^3)""u(m/s)"&
                    &"P(N/m^2)""U1""U2""U3""DE1""DE2""DE3"'
      write(unit,*) 'ZONE T="Exact Isentropic Nozzle Solution"'
      write(unit,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE&
                    & DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'

      derho   = zero
      derhou  = zero
      derhoet = zero

      do i = 1, cells
        cons_exact = primitive_to_conserved_1D(soln_exact(:,i))

        if ( present(cons_cc_ex) ) then
          cons_cc_ex(:,i) = cons_exact
        end if

        derho   = derho   + abs(cons_cc(1,i)-cons_exact(1))
        derhou  = derhou  + abs(cons_cc(2,i)-cons_exact(2))
        derhoet = derhoet + abs(cons_cc(3,i)-cons_exact(3))

        write(unit,*) x_cc(i), area_cc(i),                                     &
                      soln_exact(1,i), soln_exact(2,i), soln_exact(3,i),       &
                      cons_exact(1), cons_exact(2), cons_exact(3),             &
                      cons_cc(1,i)-cons_exact(1), cons_cc(2,i)-cons_exact(2),  &
                      cons_cc(3,i)-cons_exact(3)
      end do

      close(unit)

    else
      do i = 1, cells
        cons_cc_ex(:,i) = primitive_to_conserved_1D(soln_exact(:,i))
      end do
    end if


    if ( .not. present(cons_cc_ex) ) then
      write(*,*) 'L1 DE of rho, rho*u rho*et'
      write(*,*) derho/real(cells,dp), derhou/real(cells,dp), &
        derhoet/real(cells,dp)
    end if

  end subroutine calculate_exact_soln

!=============================== mach_from_area ==============================80
!
! Helper function to calculate Mach from area
! FIXME: Consider passing mach to speed convergence
!
!=============================================================================80
  function mach_from_area(a_over_a_star, mach_init, mach_flag)

    use set_precision,   only : dp
    use set_constants,   only : half, one, two
    use fluid_constants, only : xgm1, xgp1, gm1, gp1

    real(dp), intent(in)    :: A_over_A_star !Ratio of local area to throat area
    real(dp), intent(inout) :: mach_init  ! Initial Mach number
    integer,  intent(in)    :: mach_flag  ! Mach number flag: = 0 for subsonic,
                                          !                   = 1 for supersonic
    real(dp) :: mach_from_area

    integer  :: niter            ! outer iterations for Mach-Area relation
    integer  :: niter_max = 500  ! Max outer iterations
    integer  :: kiter            ! inner iterations for Mach-Area relation
    integer  :: kiter_max = 1000 ! Max inner iterations for Mach-Area relation

    real(dp) :: omega = 1.0_dp   ! under-relaxation factor
    real(dp) :: psi              ! Temporary function
    real(dp) :: funct            ! Function we are trying to drive to zero
    real(dp) :: dfdm             ! Derivative of function w.r.t. Mach number
    real(dp) :: tolerance = 1.e-10_dp ! Tolerance for Newton iteration
    real(dp) :: mach      = 0.5_dp    ! Mach number

    continue

    do niter = 1, niter_max
      do kiter = 1, kiter_max
        psi = two*xgp1*( one + half*gm1*mach**2 )
        funct = psi**(gp1*xgm1) - A_over_A_star**2*mach**2
        dfdm = two*mach*( psi**( two*xgm1 ) - A_over_A_star**2 )
        mach = mach - omega*funct/dfdm
        mach = abs(mach)
        if(abs(funct).le.tolerance) exit  ! Newton iteration has converged
      end do

! Check to make sure appropriate root (subsonic or supersonic) has been chosen
      if( ((mach_flag == 0).and.(mach > 1.0_dp)) .or.                        &
          ((mach_flag == 1).and.(mach < 1.0_dp)) ) then
        if( A_over_A_star < (1.2_dp) ) then
          mach = one - (mach - one)
        else
          mach_init = mach_init*( 0.8_dp*(real(1-mach_flag, dp))             &
                    + 1.2_dp*real(mach_flag, dp) )
          mach = mach_init
        end if
        omega = 0.975_dp*omega
      else
        exit  ! Converged solution is correct root (subsonic or supersonic)
      end if
    end do

    mach_from_area = mach

  end function mach_from_area

!=========================== subsonic_mach_at_exit ===========================80
!
! Finds the subsonic exit Mach number based on pressure and area ratios
!
!=============================================================================80
  pure function subsonic_mach_at_exit(ratio) result(mach)

    use set_precision,   only : dp
    use set_constants,   only : half, one, two
    use fluid_constants, only : xgp1, gp1xgm1, gm1

    real(dp), intent(in) :: ratio
    real(dp)             :: mach

    integer  :: i
    real(dp) :: mach_new, func, df_dM
    real(dp),parameter :: tolerance = 1.e-10_dp ! Tolerance for Newton iteration

    mach = half

    do i = 1,1000
      func = one/sqrt(mach**2+half*gm1*mach**4) &
           - ratio/(two*xgp1)**(half*gp1xgm1)
      df_dM = -half/((mach**2+half*gm1*mach**4)*sqrt(mach**2+half*gm1*mach**4))&
            * two*mach+two*gm1*mach**3

      mach_new = mach - func/df_dM
      if (abs(mach_new - mach) < tolerance) then
        mach = mach_new
        exit
      end if

      mach = mach_new
    end do

  end function subsonic_mach_at_exit

!=============================== pre_shock_mach ==============================80
!
! Finds the pre-shock Mach number based on pressure and area ratios
!
!=============================================================================80
  pure function pre_shock_mach(ratio) result(mach)

    use set_precision,   only : dp
    use set_constants,   only : half, one, two, four
    use fluid_constants, only : gamma, xgp1, xgm1, gxgm1, gp1xgm1, gm1, gp1

    real(dp), intent(in) :: ratio
    real(dp)             :: mach

    integer  :: i
    real(dp) :: mach_new, func, df_dM

    mach = 1.5_dp

    do i = 1,1000
      func = (gp1*mach**2/(gm1*mach**2+two))**gxgm1 &
           * ((two*gamma*mach**2-gm1)/gp1)**(-xgm1) - ratio

      df_dm = gxgm1*(two*gp1*mach/(gm1*mach**2+two) &
            - two*gp1*gm1*mach**3/(gm1*mach**2+two)**2)**xgm1 &
            * ((two*gamma*mach**2-gm1)/gp1)**(-xgm1) &
            - xgm1*(gp1*mach**2/(gm1*mach**2+two))**(gxgm1) &
            * ((two*gamma*mach**2-gm1)/gp1)**(-xgm1-one) &
            * (four*gamma*mach*xgp1)

      mach_new = mach - func/df_dM

      if (abs(mach_new - mach) < 1.0e-10_dp) then
        mach = mach_new
        exit
      end if

      mach = mach_new

    end do

  end function pre_shock_mach

!================================ estimate_te ================================80
!
! This routine uses three point Gauss quadrature to reconstruct a cell
! average approximation from the exact solution
!
!=============================================================================80
  subroutine estimate_te( cells, faces, a_star, area_f, dadx_cc, dx, x_f )

    use set_precision,   only : dp
    use set_constants,   only : zero, half, one
    use fluid_constants, only : r, gamma, gm1, gxgm1
    use solvers,         only : iterations
    use residual,        only : create_residual
    use bc,              only : subsonic_inflow_explicit, outflow_explicit
    use initialize_soln, only : mref, to, po
    use SAM,             only : SAM_Set_Problem, SAM_Make_Positive_NonZero,    &
                                SAM_Write_Data, SAM_Adapt_Mesh

    integer,                        intent(in) :: cells, faces
    real(dp),                       intent(in) :: a_star
    real(dp), dimension(faces),     intent(in) :: area_f, x_f
    real(dp), dimension(0:cells+1), intent(in) :: dadx_cc, dx

    integer  :: new_grid, grid_te_unit, faces_te, cells_te, face, cell, i, j
    real(dp) :: psi, p, t

    real(dp), dimension(faces)         :: xout
    real(dp), dimension(cells)         :: ind_in, ind_out
    real(dp), dimension(0:1+3*cells)   :: x_te, area_te, soln_te
    real(dp), dimension(3,0:1+cells)   :: cons_cc, prim_cc, te
    real(dp), dimension(3,0:1+3*cells) :: cons_gq_te, cons_dummy

    continue

    cons_dummy = zero

    grid_te_unit = find_available_unit()

    open(grid_te_unit, file='q1d_te.grd', status='old')

    read(grid_te_unit,*) faces_te

    do face = 1, faces_te
      read(grid_te_unit, *) x_te(face), area_te(face)
    end do

    close(grid_te_unit)

    cells_te = 3*cells

    call calculate_exact_soln( cells_te, x_te, area_te, a_star, x_te(faces_te),&
                               cons_dummy, cons_gq_te )

    do i = 1, cells
      j = 3*i - 3

      cons_cc(:,i) = (5.0_dp/9.0_dp)*cons_gq_te(:,j-1) &
                   + (8.0_dp/9.0_dp)*cons_gq_te(:,j)   &
                   + (5.0_dp/9.0_dp)*cons_gq_te(:,j+1)
      cons_cc(:,i) = 0.5_dp*cons_cc(:,i)
      prim_cc(:,i) = conserved_to_primitive_1D( cons_cc(:,i) )
    end do

    psi = one + half*gm1*mref**2
    t   = to/psi
    p   = po/(psi**gxgm1)

    prim_cc(1,0) = p/(r*t)
    prim_cc(2,0) = mref*sqrt(gamma*r*t)
    prim_cc(3,0) = p

    cons_cc(:,0) = primitive_to_conserved_1D(prim_cc(:,1))

    call subsonic_inflow_explicit(prim_cc(:,0), prim_cc(:,1), prim_cc(:,2))
    call outflow_explicit( prim_cc(:,cells+1),                                 &
                           prim_cc(:,cells), prim_cc(:,cells-1) )
    cons_cc(:,cells+1) = primitive_to_conserved_1D(prim_cc(:,cells+1))

    call create_residual( cells, faces, iterations, prim_cc, area_f, dadx_cc,  &
                          dx, te )

    write(*,*) te(0,2)

    grid_te_unit = find_available_unit()

    open(grid_te_unit, file='q1d_te.dat', status='replace')
    write(grid_te_unit,*) 'TITLE = "Quasi-1D Nozzle: TE from exact Solution"'
    write(grid_te_unit,*) 'variables="x(m)""TE1""TE2""TE3"'
    write(grid_te_unit,*) 'ZONE T="TE From Exact Nozzle Solution"'
    write(grid_te_unit,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE )'

    do cell = 1, cells
      write(grid_te_unit, *) x_te(3*(cell-1)), te(1, cell), te(2, cell), te(3, cell)
    end do
    close(grid_te_unit)

    ind_in = te(1, 1:cells)

    ! Set problem for SAM
    call SAM_Set_Problem(faces)
    ! Truncate and shift Wt Function !
    call SAM_Make_Positive_NonZero(ind_in)

    ! Write initial data !
    call SAM_Write_Data(x_f, ind_in, fn=80, zn=0)
    ! Adapt Mesh
    call SAM_Adapt_Mesh(x_f, ind_in, xout, ind_out)

    open(newunit=new_grid, file="nodes.dat", status="replace")
    write(new_grid,*) cells
    do i = 1, faces
      write(new_grid,*) xout(i)
    end do
    close(new_grid)

  end subroutine estimate_te

  include 'speed_of_sound.f90'
  include 'primitive_to_conserved_1D.f90'
  include 'conserved_to_primitive_1D.f90'
  include 'find_available_unit.f90'

end module solution_error
