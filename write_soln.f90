module write_soln

  implicit none

  private

  public :: write_restart

contains

!=============================================================================80
!
! Writes basic restart file
!
!=============================================================================80
  subroutine write_restart(cells, prim_cc)

    use set_precision,   only : dp
    use initialize_soln, only : L1_init, L2_init, Linf_init

    implicit none

    integer,                        intent(in) :: cells
    real(dp), dimension(3,cells+2), intent(in) :: prim_cc

    integer :: restart_unit, cell

    continue

    restart_unit = find_available_unit()
    open(restart_unit, file='q1d.rst', status='replace')

    write(restart_unit,*) L1_init(1), L1_init(2), L1_init(3)
    write(restart_unit,*) L2_init(1), L2_init(2), L2_init(3)
    write(restart_unit,*) Linf_init(1), Linf_init(2), Linf_init(3)

    do cell = 1, cells+2
      write(restart_unit,*) prim_cc(1, cell), prim_cc(2, cell), prim_cc(3,cell)
    end do

    close(restart_unit)

  end subroutine write_restart

!=============================================================================80
!
! Writes solution as simple X-Y line output
!
!=============================================================================80

  subroutine write_soln_line()

  end subroutine write_soln_line

!=============================================================================80
!
! Writes solution as cell output ( fancy color plot )
!
!=============================================================================80

  subroutine write_soln_cells()

  end subroutine write_soln_cells

  include 'find_available_unit.f90'

end module write_soln
