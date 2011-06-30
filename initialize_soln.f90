module initialize_soln

  implicit none

  private

  public :: prim_cc
  public :: cons_cc

  real(dp), allocatable, dimension(:,:) :: prim_cc   ! primitive vars
  real(dp), allocatable, dimension(:,:) :: cons_cc   ! conserved vars
  real(dp), allocatable, dimension(:,:) :: cons_cc_0 ! old iter conserved vars

contains

  subroutine allocate_soln(neq, cells)

    implicit none

    integer, intent(in) :: cells

    continue

    allocate( prim_cc(neq, cells+2), cons_cc(neq, cells+2) )

  end subroutine allocate_soln

  subroutine initial_soln(cells)

    use set_constants,   only : half, one
    use XXXXX,           only : restart
    use fluid_constants, only : r, gamma, gm1, gxgm1

    implicit none

    integer, intent(in) :: cells

    continue

    if (restart) then

      soln_unit = find_available_unit()

      open(soln_unit, file='q1d.sln', status='old')

      do cell = 1, cells+2
        read(soln_unit,*) prim_cc(1, cell), prim_cc(2, cell), prim_cc(3,cell)
        cons_cc(:,cell) = primitive_to_conserved_1D( prim_cc(:,cell) )
      end do

      close(soln_unit)

    else

      psi = one + half*gm1*mref*mref
      t   = to/psi
      p = po/(psi**gxgm1)

      do cell = 1, cells+2
        prim_cc(1,cell) = p/(r*t)
        prim_cc(2,cell) = mref*sqrt(gamma*r*t)
        prim_cc(3,cell) = p

        cons_cc(:,cell) = primitive_to_conserved_1D( prim_cc(:,cell) )
      end do

    endif

  end subroutine initial_soln

  include 'find_available_unit.f90'
  include 'primitive_to_conserved_1D.f90'

end module initialize_soln
