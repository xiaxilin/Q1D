! This program will create the grid for the Q1D nozzle solver

! Needs to know the # of cells & nozzle length, otherwise, reads x coords

program q1d_grid

  use set_precision, only : dp
  use set_constants, only : pi, set_pi

  implicit none

  integer  :: i, j,  cells, faces, available_unit, find_available_unit
  real(dp) :: nozzlelength, dx
  real(dp), allocatable, dimension(:) :: x, areaface, areacent, dadx
  real(dp), allocatable, dimension(:) :: x_te, areaface_te

  continue

  call set_pi()

  write(*,*) "How many cells to create?  Enter 0 to read from file."
  read(*,*) cells

  if ( cells > 0 ) then

    write(*,*) "What is the length of the nozzle"
    read(*,*) nozzlelength

    faces = cells+1
    allocate( x(faces) )

    dx = nozzlelength/real(cells,dp)

    do i = 1, faces
      x(i) = -0.5_dp*nozzlelength + dx*real(i-1,dp)
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

  do i = 1, faces
    areaface(i) = 0.2_dp + 0.4_dp * (1.0_dp + sin(pi*(x(i) - 0.5_dp)))
! From Yee et al., JCP 1985
!    areaface(i) = 1.398_dp + 0.347_dp*tanh(0.8_dp*x(i)-4.0_dp)
  end do
  do i = 1, 2+3*(faces-1)
    areaface_te(i) = 0.2_dp + 0.4_dp * (1.0_dp + sin(pi*(x_te(i) - 0.5_dp)))
! From Yee et al., JCP 1985
!    areaface_te(i) = 1.398_dp + 0.347_dp*tanh(0.8_dp*x_te(i)-4.0_dp)
  end do

  do i = 2, cells+1
    areacent(i) = 0.2_dp + 0.4_dp &
                * (1.0_dp + sin(pi*( 0.5_dp*(x(i)+x(i-1)) - 0.5_dp )))
    dadx(i) = 0.4_dp * pi * cos(pi*( 0.5_dp*(x(i)+x(i-1)) - 0.5_dp ))
!    areacent(i) = 1.398_dp + 0.347_dp*tanh(0.4_dp*(x(i)+x(i-1))-4.0_dp)
!    dadx(i) = 0.8_dp*0.347_dp/cosh(0.4_dp*(x(i)+x(i-1))-4.0_dp)**2
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

  close(available_unit)

  deallocate( x, areaface, areacent, dadx )

end program q1d_grid

include 'find_available_unit.f90'
