! The Q1D-nozzle has an exact solution from the mach/area relations
! This module will hold the routines necessary to calculate the exact soln,
! form the DE error norms, and calculate the TE.

module solution_error



subroutine calculate_exact_soln(cells, throat_area, area, x_cc)

  use set_precision,   only : dp
  use set_constants,   only : half, one
  use fluid_constants, only : gm1, gxgm1, R, To, Po

  implicit none

  integer,                      intent(in) :: imax
  real(dp),                     intent(in) :: throat_area
  real(dp), dimension(cells+2), intent(in) :: area_cc
  real(dp), dimension(cells+2), intent(in) :: x_cc

  real(dp) :: Mach_From_Area
  real(dp) :: asnd      ! Define type for speed of sound function
  real(dp) :: psi       ! ( = T_0/T )
  real(dp) :: temp      ! (Temperature, T)
  real(dp) :: U1        ! entries of cons. var. vector
  real(dp) :: U2        ! entries of cons. var. vector
  real(dp) :: U3        ! entries of cons. var. vector
  real(dp) :: mach_init ! Initial Mach number
  real(dp), dimension(imax)   :: mach_exact
  real(dp), dimension(3,imax) :: soln_exact

  continue

! calculate exact mach/area solution
  mach_init = 0.1_dp
! subsonic section to throat
  do i = 1, i_throat
    mach_exact(i) = mach_from_area(area(i)/throat_area, mach_init, 0)
    mach_init = mach_exact(i)
  end do
! sub/supersonic section after throat
  mach_init = mach_init + 0.01
  do i = i_throat+1, imax
    mach_exact(i) = mach_from_area(area(i)/throat_area, mach_init, 1)
    mach_init = mach_exact(i)
  end do

! once the mach/area solution has been found, calculate the primitive variables
  do i = 1, imax
    psi = one + half*gm1*mach_exact(i)**2
    temp = To / psi
    v_exact(3,i) = Po / ( psi**gxgm1 ) )
    v_exact(1,i) = v_exact(3,i)/( R*temp )
    v_exact(2,i) = mach_exact(i)*speed_of_sound(v_exact(3,i),v_exact(1,i))
  end do

! set up output file for exact solution
  open(57,file='Exact-Solution.dat',status='unknown')
  write(57,*) 'TITLE = "Quasi-1D Nozzle: Exact Isentropic Solution"'
  write(57,*) 'variables="x(m)""Area(m^2)""rho(kg/m^3)""u(m/s)""Press(N/m^2)"&
& "Mach""U1""U2""U3"'
  write(57,*) 'ZONE T="Exact Isentropic Nozzle Solution"'
  write(57,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 
  do i = 2, cells+1
    U1 = v_exact(1,i)
    U2 = v_exact(1,i)*v_exact(2,i)
    U3 = v_exact(3,i)*xgm1 + half*v_exact(1,i)*v_exact(2,i)**2
    write(57,*) x_cc(i), area(i), v_exact(1,i), v_exact(2,i), v_exact(i,3),    &
      mach_exact(i,4), U1, U2, U3
  end do

  close(57)

end subroutine calculate_exact_soln


!******************************************************************************

function mach_from_area(A_over_A_star, mach_init, mach_flag)

  use set_precision,   only : dp
  use fluid_constants, only : xgm1, xgp1, gm1, gp1

  real(dp) :: Mach_From_Area  
  real(dp) :: A_over_A_star    ! Ratio of local area to throat area
  real(dp) :: mach_init        ! (INPUT) Initial Mach number
  integer  :: mach_flag        ! (INPUT) Mach number flag: = 0 for subsonic,
                               !                           = 1 for supersonic

  integer  :: niter            ! outer iterations for Mach-Area relation
  integer  :: niter_max = 500  ! Maximum outer iterations
  integer  :: kiter            ! inner iterations for Mach-Area relation
  integer  :: kiter_max = 1000 ! Maximum inner iterations for Mach-Area relation

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
    if( ((mach_flag.eq.0).and.(mach.gt.1.0)) .or.                              &
        ((mach_flag.eq.1).and.(mach.lt.1.0)) ) then
      if( A_over_A_star.le.(1.2_dp) ) then
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

end module solution_error
