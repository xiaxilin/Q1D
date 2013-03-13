! This program will create the grid for the Q1D nozzle solver

! Needs to know the # of cells & nozzle length, otherwise, reads x coords

module q1d_grid_functions

  use set_precision, only : dp

  implicit none

  private

  abstract interface

    pure function areaf(x)

      use set_precision, only : dp

      real(dp), intent(in) :: x
      real(dp)             :: areaf

    end function areaf

  end interface

  abstract interface

    pure function dadxf(x)

      use set_precision, only : dp

      real(dp), intent(in) :: x
      real(dp)             :: dadxf

    end function dadxf

  end interface

  public :: areaf, dadxf, choose_area_distribution

contains

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine choose_area_distribution(type, area, dadx)

    character(*),              intent(in)  :: type
    procedure(areaf), pointer, intent(out) :: area
    procedure(dadxf), pointer, intent(out) :: dadx

    continue

    area => null()
    dadx => null()

    select case(trim(type))
    case('sin')
      area => sinx_nozzle
      dadx => da_sinx_nozzle
    case('ext-sin')
      area => sinx_nozzle
      dadx => da_sinx_nozzle
    case('cosh')
      area => coshx_nozzle
      dadx => da_coshx_nozzle
    case('bump')
      area => bumpx_nozzle
      dadx => da_bumpx_nozzle
    case('tyrone')
      area => tyrone_nozzle
      dadx => da_tyrone_nozzle
    end select

  end subroutine choose_area_distribution

!=============================================================================80
!
!
!
!=============================================================================80
  pure function sinx_nozzle(x) result(area)

    use set_constants, only : pi

    real(dp), intent(in) :: x
    real(dp)             :: area

    continue

    if (x < -1.0_dp .or. x > 1.0) then
      area = 1.0_dp
    else
      area = 0.2_dp + 0.4_dp * ( 1.0_dp + sin( pi*( x - 0.5_dp ) ) )
    end if

  end function sinx_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function da_sinx_nozzle(x) result(dadx)

    use set_constants, only : pi

    real(dp), intent(in) :: x
    real(dp)             :: dadx

    continue

    if (x < -1.0_dp .or. x > 1.0) then
      dadx = 0.0_dp
    else
      dadx = 0.4_dp * pi * cos( pi*( x - 0.5_dp ) )
    end if

  end function da_sinx_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function coshx_nozzle(x) result(area)

    real(dp), intent(in) :: x
    real(dp)             :: area

    continue

    area = 0.8_dp + cosh(x)

  end function coshx_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function da_coshx_nozzle(x) result(dadx)

    real(dp), intent(in) :: x
    real(dp)             :: dadx

    continue

    dadx = sinh(x)

  end function da_coshx_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function bumpx_nozzle(x) result(area)

    real(dp), intent(in) :: x
    real(dp)             :: area

    continue

    area = 1.0_dp
    if ( x > -1.0_dp .and. x < 1.0_dp ) then
      area = 1.0_dp - exp(-1.0_dp/(1.0_dp-x**2))
    end if

  end function bumpx_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function da_bumpx_nozzle(x) result(dadx)

    real(dp), intent(in) :: x
    real(dp)             :: dadx

    continue

    dadx = 0.0_dp
    if ( x > -1.0_dp .and. x < 1.0_dp ) then
      dadx = -2.0_dp*x*exp(1.0_dp/(x**2-1.0_dp))
      dadx = dadx / (x**2 - 1.0_dp)
    end if

  end function da_bumpx_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function yee_nozzle(x) result(area)

    real(dp), intent(in) :: x
    real(dp)             :: area

    continue
! From Yee et al., JCP 1985
    area = 1.398_dp + 0.347_dp * tanh( 0.8_dp*x - 4.0_dp )

  end function yee_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function da_yee_nozzle(x) result(dadx)

    real(dp), intent(in) :: x
    real(dp)             :: dadx

    continue

    dadx = 0.8_dp*0.347_dp / cosh( 0.8_dp*x - 4.0_dp )**2

  end function da_yee_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function tyrone_nozzle(x) result(area)

    real(dp), intent(in) :: x
    real(dp)             :: area

    real(dp) :: a, c

    continue

    a = 2.0_dp
    c = 0.3_dp

    area = -a*exp(-x**2/(2.0_dp*c**2)) + tanh(4.0_dp*x) + 1.9_dp*a

  end function tyrone_nozzle

!=============================================================================80
!
!
!
!=============================================================================80
  pure function da_tyrone_nozzle(x) result(dadx)

    real(dp), intent(in) :: x
    real(dp)             :: dadx

    real(dp) :: a, c

    continue

    a = 2.0_dp
    c = 0.3_dp

    dadx = a*x*exp(-x**2/(2.0_dp*c**2))/(c**2) + 4.0_dp/cosh(4.0_dp*x)**2

  end function da_tyrone_nozzle

end module q1d_grid_functions

program q1d_grid

  use set_precision,      only : dp
  use set_constants,      only : set_pi
  use q1d_grid_functions, only : choose_area_distribution, areaf, dadxf

  implicit none

  integer       :: i, j,  cells, faces, available_unit, find_available_unit
  real(dp)      :: nozzlelength, dx
  character(10) :: nozzle_func

  real(dp), allocatable, dimension(:) :: x, areaface, areacent, dadx
  real(dp), allocatable, dimension(:) :: x_te, areaface_te

  procedure(areaf), pointer :: area_func
  procedure(dadxf), pointer :: dadx_func

  continue

  call set_pi()

  write(*,*) "What nozzle function?"

  nozzle_func = 'sin'

  call choose_area_distribution(nozzle_func, area_func, dadx_func)

  write(*,*) "How many cells to create?  Enter 0 to read from file."
  read(*,*) cells

  if ( cells > 0 ) then

    write(*,*) "What is the length of the nozzle"
    read(*,*) nozzlelength

    faces = cells+1
    allocate( x(faces) )

    dx = nozzlelength/real(cells,dp)

! Special extended nozzle
!    dx = (2.0_dp+nozzlelength)/real(cells,dp)

    do i = 1, faces
      x(i) = -0.5_dp*nozzlelength + dx*real(i-1,dp)
! Special extended nozzle
!      x(i) = -10.0_dp -0.5_dp*nozzlelength + dx*real(i-1,dp)
! For supersonic diffuser
!      x(i) = dx*real(i-1,dp)
    end do

  else

    available_unit = find_available_unit()
    open(available_unit, file='nodes.dat', status='old')
    read(available_unit,*) cells

    faces = cells+1
    allocate( x(faces) )

    do i = 1, faces
      read(available_unit,*) x(i)
    end do
    close(available_unit)

  end if

! Now that the face locations are set up, calculate areas and dA/dx
  allocate( areaface(faces), areacent(cells+1), dadx(cells+1) )

! Arrays for TE estimation
  allocate( x_te(2+3*(faces-1)), areaface_te(2+3*(faces-1)) )

  x_te(1) = x(1)
  do i = 2, faces
    j = i*3 - 3
    dx = x(i)-x(i-1)

    x_te(j)   = 0.5_dp*dx + x(i-1)
    x_te(j-1) = x_te(j) - 0.5_dp*dx*sqrt(3.0_dp/5.0_dp)
    x_te(j+1) = x_te(j) + 0.5_dp*dx*sqrt(3.0_dp/5.0_dp)

  end do
  x_te(2+3*(faces-1)) = x(faces)

! Calculate the face areas
  do i = 1, faces
    areaface(i) = area_func(x(i))
  end do
  do i = 1, 2+3*(faces-1)
    areaface_te(i) = area_func(x_te(i))
  end do

! Calculate the cell center areas and slopes
  do i = 2, cells+1
    areacent(i) = area_func(0.5_dp*(x(i)+x(i-1)))
    dadx(i) = dadx_func(0.5_dp*(x(i)+x(i-1)))
  end do
  areacent(1) = areacent(2)
  dadx(1) = dadx(2)

! Write data to file
  available_unit = find_available_unit()

  open(available_unit, file='q1d.grd', status='replace')
  write(available_unit,*) cells

  do i = 1, faces
    write(available_unit,*) x(i), areaface(i), areacent(i), dadx(i)
  end do

  close(available_unit)

! Write TE grid data to file
  available_unit = find_available_unit()

  open(available_unit, file='q1d_te.grd', status='replace')
  write(available_unit,*) 2+3*(faces-1)

  do i = 1, 2+3*(faces-1)
    write(available_unit,*) x_te(i), areaface_te(i)
  end do

  write(*,*) 2+3*(faces-1), x_te( 2+3*(faces-1) ), areaface_te( 2+3*(faces-1) )

  close(available_unit)

  deallocate( x, areaface, areacent, dadx, x_te, areaface_te )

end program q1d_grid

include 'find_available_unit.f90'
