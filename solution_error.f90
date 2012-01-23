! The Q1D-nozzle has an exact solution from the mach/area relations
! This module will hold the routines necessary to calculate the exact soln,
! form the DE error norms, and calculate the TE.

module solution_error

  implicit none

  private

  public :: calculate_exact_soln

contains

  subroutine calculate_exact_soln(cells, x_cc, area_cc, a_star, cons_cc)

    use set_precision,   only : dp
    use set_constants,   only : half, one, two
    use fluid_constants, only : gm1, xgm1, gxgm1, r
    use initialize_soln, only : to, po

    implicit none

    integer,                      intent(in) :: cells
    real(dp), dimension(cells+2), intent(in) :: x_cc
    real(dp), dimension(cells+2), intent(in) :: area_cc
    real(dp),                     intent(in) :: a_star
    real(dp), dimension(3,cells+2), intent(in) :: cons_cc

    integer               :: i, i_throat
    integer, dimension(1) :: i_min      ! needed for minloc function

    real(dp) :: asnd      ! Define type for speed of sound function
    real(dp) :: psi       ! ( = T_0/T )
    real(dp) :: temp      ! (Temperature, T)
    real(dp) :: mach_init ! Initial Mach number
    real(dp), dimension(3)         :: cons_exact
    real(dp), dimension(cells+2)   :: mach_exact
    real(dp), dimension(3,cells+2) :: soln_exact

    continue

    i_min = minloc(area_cc(2:cells+1))
    print*, i_min
    i_throat = i_min(1)

! calculate exact mach/area solution
    mach_init = 0.1_dp
! subsonic section to throat
    do i = 2, i_throat
      mach_exact(i) = mach_from_area(area_cc(i)/a_star, mach_init, 0)
      mach_init = mach_exact(i)
    end do
! sub/supersonic section after throat
    mach_init = mach_init + 0.01
    do i = i_throat+1, cells+1
      mach_exact(i) = mach_from_area(area_cc(i)/a_star, mach_init, 1)
      mach_init = mach_exact(i)
    end do

! once the isentropic solution has been found, check for shocked case
! need Pb and Ae for this

! once the mach/area solution has been found, calculate the primitive variables
    do i = 2, cells+1
      psi = one + half*gm1*mach_exact(i)**2
      temp = to / psi
      soln_exact(3,i) = po / ( psi**gxgm1 )
      soln_exact(1,i) = soln_exact(3,i)/( r*temp )
      soln_exact(2,i) = mach_exact(i) &
                      * speed_of_sound(soln_exact(3,i), soln_exact(1,i))
    end do

! set up output file for exact solution
    open(57,file='Exact-Solution.tec',status='unknown')
    write(57,*) 'TITLE = "Quasi-1D Nozzle: Exact Isentropic Solution"'
    write(57,*) 'variables="x(m)""Area(m^2)""rho(kg/m^3)""u(m/s)""Press(N/m^2)"&
              & "U1""U2""U3""DE1""DE2""DE3"'
    write(57,*) 'ZONE T="Exact Isentropic Nozzle Solution"'
    write(57,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE&
               & DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'

    do i = 2, cells+1
      cons_exact = primitive_to_conserved_1D(soln_exact(:,i))
      write(57,*) x_cc(i), area_cc(i),                                         &
                  soln_exact(1,i), soln_exact(2,i), soln_exact(3,i),           &
                  cons_exact(1), cons_exact(2), cons_exact(3),                 &
                  cons_cc(1,i)-cons_exact(1), cons_cc(2,i)-cons_exact(2), &
                  cons_cc(3,i)-cons_exact(3)
    end do

    close(57)

  end subroutine calculate_exact_soln


!******************************************************************************

  function mach_from_area(a_over_a_star, mach_init, mach_flag)

    use set_precision,   only : dp
    use set_constants,   only : half, one, two
    use fluid_constants, only : xgm1, xgp1, gm1, gp1

    implicit none

    real(dp) :: Mach_From_Area
    real(dp) :: A_over_A_star    ! Ratio of local area to throat area
    real(dp) :: mach_init        ! (INPUT) Initial Mach number
    integer  :: mach_flag        ! (INPUT) Mach number flag: = 0 for subsonic,
                                 !                           = 1 for supersonic

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
      if( ((mach_flag == 0).and.(mach > 1.0)) .or.                            &
          ((mach_flag == 1).and.(mach < 1.0)) ) then
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

  include 'speed_of_sound.f90'
  include 'primitive_to_conserved_1D.f90'

end module solution_error
