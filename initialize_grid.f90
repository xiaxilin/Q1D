! Holds routines to read and allocate the grid

module initialize_grid

  use set_precision, only : dp

  implicit none

  private

  public :: cells
  public :: faces
  public :: dxsi
  public :: area_f
  public :: area_cc
  public :: x_f
  public :: x_cc
  public :: dadx_cc
  public :: dxdxsi_cc

  public :: read_grid
  public :: deallocate_grid

  integer  :: cells ! # of interior cells
  integer  :: faces ! # of cell faces
  real(dp) :: dxsi  ! curvilinear stretching term

  real(dp), allocatable, dimension(:) :: area_f  ! face area
  real(dp), allocatable, dimension(:) :: area_cc ! cell center area
  real(dp), allocatable, dimension(:) :: x_f     ! face node location
  real(dp), allocatable, dimension(:) :: x_cc    ! cell center node location
  real(dp), allocatable, dimension(:) :: dadx_cc ! cell center da/dx
  real(dp), allocatable, dimension(:) :: dxdxsi_cc ! cell dx/dxsi metric

contains

!=============================================================================80
!
!=============================================================================80

  subroutine read_grid

    use set_constants, only : half, one

    implicit none

    integer  :: grid_unit, face
    real(dp) :: dx

    continue

    grid_unit = find_available_unit()

    open(grid_unit, file='q1d.grd', status='old')

    read(grid_unit,*) cells

    faces = cells + 1

    dxsi = one/real(cells, dp)

    call allocate_grid()

    do face = 1, faces
      read(grid_unit, *) x_f(face), area_f(face), area_cc(face), dadx_cc(face)
    end do

    close(grid_unit)

    do face = 2, faces
      dx = x_f(face) - x_f(face-1)
      dxdxsi_cc(face) = dx/dxsi
      x_cc(face)  = x_f(face-1) + half*dx
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
    allocate( dadx_cc(cells+2), dxdxsi_cc(cells+2))

    area_f  = zero
    x_f     = zero
    area_cc = zero
    x_cc    = zero
    dadx_cc = zero
    dxdxsi_cc = zero


  end subroutine allocate_grid

!=============================================================================80
!
!=============================================================================80

  subroutine deallocate_grid

    implicit none

    continue

    deallocate( area_f, x_f, area_cc, x_cc, dadx_cc, dxdxsi_cc)

  end subroutine deallocate_grid

  include 'find_available_unit.f90'

end module initialize_grid
