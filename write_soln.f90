module write_soln

  implicit none

contains

!=============================================================================80
!
! Writes basic restart file
!
!=============================================================================80
  subroutine write_restart(cells, prim_cc)

    use set_precision, only : dp

    implicit none

    integer,                        intent(in) :: cells
    real(dp), dimension(3,cells+2), intent(in) :: prim_cc

    continue

    restart_unit = find_available_unit()
    open(restart_unit, file="q1d.rst", status=replace)

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
