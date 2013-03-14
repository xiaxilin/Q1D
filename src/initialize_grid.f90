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
  public :: dx
  public :: dadx_cc
  public :: cell_vol

  public :: read_grid
  public :: deallocate_grid

  integer  :: cells ! # of interior cells
  integer  :: faces ! # of cell faces

  real(dp), allocatable, dimension(:) :: area_f   ! face area
  real(dp), allocatable, dimension(:) :: area_cc  ! cell center area
  real(dp), allocatable, dimension(:) :: x_f      ! face node location
  real(dp), allocatable, dimension(:) :: x_cc     ! cell center node location
  real(dp), allocatable, dimension(:) :: dadx_cc  ! cell center da/dx
  real(dp), allocatable, dimension(:) :: dx       ! cell dx
  real(dp), allocatable, dimension(:) :: cell_vol ! cell area/volume

contains

!================================= read_grid =================================80
!
! Reads the grid file
! FIXME: consider removing dadx_cc and area_cc
!
!=============================================================================80
  subroutine read_grid

    use set_constants, only : half

    integer  :: grid_unit, face, cell

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

    do cell = 2, cells+1
      dx(cell)        = x_f(cell) - x_f(cell-1)
      x_cc(cell)      = x_f(cell-1) + half*dx(cell)
      cell_vol(cell)  = dx(cell)*area_cc(cell)
    end do
    dx(1)       = dx(2)
    dx(cells+2) = dx(cells+1)
    cell_vol(1)       = cell_vol(2)
    cell_vol(cells+1) = cell_vol(cells+1) !FIXME: review this!

!    do cell = 1, cells
!      dx(cell)        = x_f(cell+1) - x_f(cell)
!      x_cc(cell)      = x_f(cell) + half*dx(cell)
!      cell_vol(cell)  = dx(cell)*area_cc(cell)
!    end do
!    dx(0)       = dx(1)
!    dx(cells+1) = dx(cells)
!    cell_vol(0)       = cell_vol(1)
!    cell_vol(cells+1) = cell_vol(cell)

  end subroutine read_grid

!=============================================================================80
!
! Allocates and initializes the grid variables
!
!=============================================================================80
  subroutine allocate_grid

    use set_constants, only : zero

    continue

    allocate( area_f(faces), x_f(faces) )
    allocate( area_cc(cells+2), x_cc(cells+2) )
    allocate( dadx_cc(cells+2), dx(cells+2), cell_vol(cells+2) )
!    allocate( area_cc(0:cells+1), x_cc(0:cells+1) )
!    allocate( dadx_cc(0:cells+1), dx(0:cells+1), cell_vol(0:cells+1) )

    area_f   = zero
    x_f      = zero
    area_cc  = zero
    x_cc     = zero
    dadx_cc  = zero
    dx       = zero
    cell_vol = zero

  end subroutine allocate_grid

!============================== deallocate_grid ==============================80
!
! Deallocates the grid variables in the correct order
!
!=============================================================================80
  subroutine deallocate_grid

    continue

    deallocate( cell_vol, dx, dadx_cc, x_cc, area_cc, x_f, area_f )

  end subroutine deallocate_grid

  include 'find_available_unit.f90'

end module initialize_grid
