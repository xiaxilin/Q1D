! Holds routines to read and allocate the grid

module initialize_grid

  use set_precision, only : dp

  implicit none

  private

  public :: cells
  public :: faces
  public :: area_f
  public :: area_cc
  public :: x_f
  public :: x_cc
  public :: dadx_cc
  public :: dx_cc

  public :: read_grid

  integer :: cells
  integer :: faces

  real(dp), allocatable, dimension(:) :: area_f  ! face area
  real(dp), allocatable, dimension(:) :: area_cc ! cell center area
  real(dp), allocatable, dimension(:) :: x_f     ! face node location
  real(dp), allocatable, dimension(:) :: x_cc    ! cell center node location
  real(dp), allocatable, dimension(:) :: dadx_cc ! cell center da/dx
  real(dp), allocatable, dimension(:) :: dx_cc   ! cell dx (length)

contains

!=============================================================================80
!
!=============================================================================80

  subroutine read_grid

    use set_constants, only : half

    implicit none

    integer :: grid_unit, face

    continue

    grid_unit = find_available_unit()

    open(grid_unit, file='q1d.grd', status='old')

    read(grid_unit,*) cells

    faces = cells + 1

    call allocate_grid()

    do face = 1, faces
      read(grid_unit, *) x_f(face), area_f(face), area_cc(face), dadx_cc(face)
    end do

    close(grid_unit)

    do face = 2, faces
      dx_cc(face) = x_f(face) - x_f(face-1)
      x_cc(face)  = x_f(face-1) + half*dx_cc(face)
    end do

  end subroutine read_grid

!=============================================================================80
!
!=============================================================================80

  subroutine allocate_grid

    use set_constants, only : zero

    implicit none

    continue

    allocate( area_f(faces), x_f(faces), area_cc(cells+2), x_cc(cells+2) )
    allocate( dadx_cc(cells+2), dx_cc(cells+2))

    area_f  = zero
    area_cc = zero
    dadx_cc  = zero
    dadx_cc = zero
    x_f     = zero
    x_cc    = zero

  end subroutine allocate_grid

  include 'find_available_unit.f90'

end module initialize_grid
