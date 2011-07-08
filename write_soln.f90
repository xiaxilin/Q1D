module write_soln

  implicit none

  private

  public :: write_restart
  public :: write_soln_line
  public :: init_write_files

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

!

  subroutine init_write_files()

    implicit none

    integer :: line_unit

    continue

    open(line_unit, file='q1d_lines.tec', status='replace')
    close(line_unit)

  end subroutine init_write_files

!=============================================================================80
!
! Writes solution as simple X-Y line output
!
!=============================================================================80

  subroutine write_soln_line(iteration, cells, x_cc, prim_cc, cons_cc)

    use set_precision, only : dp

    implicit none

    integer,                        intent(in) :: iteration, cells
    real(dp), dimension(cells+2),   intent(in) :: x_cc
    real(dp), dimension(3,cells+2), intent(in) :: prim_cc, cons_cc

    integer :: i, var, line_unit

    continue

    line_unit = find_available_unit()

    open(line_unit, file='q1d_lines.tec', access='append', status='old')

    write(line_unit,*) 'VARIABLES = "X_cc", "rho", "U", "P", "rho*U", "rho*E"'
    write(line_unit,*) 'ZONE DATAPACKING=BLOCK, I=', cells
    write(line_unit,*) 'T="', iteration,'" '
    do i = 2, cells+1
      write(line_unit,*) x_cc(i)
    end do
    do var = 1,3
      do i = 2, cells+1
        write(line_unit,*) prim_cc(var,i)
      end do
    end do
    do var = 2,3
      do i = 2, cells+1
        write(line_unit,*) cons_cc(var,i)
      end do
    end do

    close(line_unit)

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
