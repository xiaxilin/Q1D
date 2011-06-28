! This program will create the grid for the Q1D nozzle solver

! Needs to know the # of cells & nozzle length, otherwise, reads x coords

program q1d_grid

  use set_precision, only : dp
  use set_constants, only : pi, set_pi

  implicit none

  integer  :: i, cells, faces, available_unit, find_available_unit
  real(dp) :: nozzlelength, dx
  real(dp), allocatable, dimension(:) :: x, areaface, areacent, dadx

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

  do i = 1, faces
    areaface(i) = 0.2_dp + 0.4_dp * (1.0_dp + sin(pi*(x(i) - 0.5_dp)))
  end do

  do i = 1, cells
    areacent(i) = 0.2_dp + 0.4_dp &
                * (1.0_dp + sin(pi*( 0.5_dp*(x(i)+x(i+1)) - 0.5_dp )))
    dadx(i) = 0.4_dp * pi * cos(pi*( 0.5_dp*(x(i)+x(i+1)) - 0.5_dp ))
  end do
  areacent(cells+1) = areacent(cells)
  dadx(cells+1) = dadx(cells)

! Write data to file
  available_unit = find_available_unit()

  open(available_unit, file='q1d.grd', status='replace')
  write(available_unit,*) cells

  do i = 1, faces
    write(available_unit,*) x(i), areaface(i), areacent(i), dadx(i)
  end do

  close(available_unit)

  deallocate( x, areaface, areacent, dadx )

end program q1d_grid

include 'find_available_unit.f90'