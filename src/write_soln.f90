module write_soln

  implicit none

  private

  public :: write_restart
  public :: write_soln_line
  public :: init_write_files
  public :: write_entropy

contains

!=============================== write_restart ===============================80
!
! Writes basic restart file
!
!=============================================================================80
  subroutine write_restart(cells, prim_cc)

    use set_precision,   only : dp
    use initialize_soln, only : L1_init, L2_init, Linf_init

    integer,                        intent(in) :: cells
    real(dp), dimension(3,cells+2), intent(in) :: prim_cc
!    real(dp), dimension(3,0:cells+1), intent(in) :: prim_cc

    integer :: restart_unit, cell

    continue

    restart_unit = find_available_unit()
    open(restart_unit, file='q1d.rst', status='replace')

    write(restart_unit,*) L1_init(1), L1_init(2), L1_init(3)
    write(restart_unit,*) L2_init(1), L2_init(2), L2_init(3)
    write(restart_unit,*) Linf_init(1), Linf_init(2), Linf_init(3)

    do cell = 1, cells+2!0, cells+1
      write(restart_unit,*) prim_cc(1, cell), prim_cc(2, cell), prim_cc(3,cell)
    end do

    close(restart_unit)

  end subroutine write_restart

!============================== init_write_files =============================80
!
! Open the lines file so that it replaces the old
!
!=============================================================================80
  subroutine init_write_files()

    integer :: line_unit

    continue

    line_unit = find_available_unit()
    open(line_unit, file='q1d_lines.dat', status='replace')
    close(line_unit)

  end subroutine init_write_files

!============================== write_soln_line ==============================80
!
! Writes solution as simple X-Y line output
!
!=============================================================================80
  subroutine write_soln_line(iteration, cells, x_cc, prim_cc, cons_cc)

    use set_precision, only : dp

    integer,                        intent(in) :: iteration, cells
    real(dp), dimension(cells+2),   intent(in) :: x_cc
    real(dp), dimension(3,cells+2), intent(in) :: prim_cc, cons_cc
!    real(dp), dimension(0:cells+1),   intent(in) :: x_cc
!    real(dp), dimension(3,0:cells+1), intent(in) :: prim_cc, cons_cc

    integer :: i, var, line_unit

    continue

    line_unit = find_available_unit()

    open(line_unit, file='q1d_lines.dat', access='append', status='old')

    write(line_unit,*) 'VARIABLES = "X_cc", "rho", "U", "P", "rho*U", "rho*E"'
    write(line_unit,*) 'ZONE DATAPACKING=BLOCK, I=', cells
    write(line_unit,*) 'T="', iteration,'" '
    do i = 2, cells+1!1, cells
      write(line_unit,*) x_cc(i)
    end do
    do var = 1,3
      do i = 2, cells+1!1, cells
        write(line_unit,*) prim_cc(var,i)
      end do
    end do
    do var = 2,3
      do i = 2, cells+1!1, cells
        write(line_unit,*) cons_cc(var,i)
!        write(line_unit,*) conserved_to_primitive_1D(prim_cc(var,i))
      end do
    end do

    close(line_unit)

  end subroutine write_soln_line

!=============================== write_entropy ===============================80
!
! Writes entropy vars as simple X-Y line output
!
!=============================================================================80
  subroutine write_entropy( cells, x_cc, prim_cc )

    use set_precision,   only : dp
    use set_constants,   only : half
    use fluid_constants, only : gamma, xgm1

    integer,                        intent(in) :: cells
    real(dp), dimension(cells+2),   intent(in) :: x_cc
    real(dp), dimension(3,cells+2), intent(in) :: prim_cc
!    real(dp), dimension(0:cells+1),   intent(in) :: x_cc
!    real(dp), dimension(3,0:cells+1), intent(in) :: prim_cc

    integer :: i, var, entropy_unit
    real(dp), dimension(cells+2)    :: entropy
    real(dp), dimension(3, cells+2) :: ent_var
!    real(dp), dimension(0:cells+1)    :: entropy
!    real(dp), dimension(3, 0:cells+1) :: ent_var

    continue

    entropy_unit = find_available_unit()

    open(entropy_unit, file='q1d_entropy.dat', status='replace')

    write(entropy_unit,*) 'VARIABLES = "X_cc", "entropy", "S1", "S2", "S3"'
    write(entropy_unit,*) 'ZONE DATAPACKING=BLOCK, I=', cells

    do i = 2, cells+1!1, cells
      write(entropy_unit,*) x_cc(i)
    end do

    do i = 2, cells+1!1, cells
      entropy(i) = log(prim_cc(3,i)/prim_cc(1,i)**gamma)

      ent_var(1,i) = (gamma - entropy(i))*xgm1                                &
                   - half*prim_cc(1,i)*prim_cc(2,i)**2/prim_cc(3,i)
      ent_var(2,i) = prim_cc(1,i)*prim_cc(2,i)/prim_cc(3,i)
      ent_var(3,i) = -prim_cc(1,i)/prim_cc(3,i)

      write(entropy_unit, *) entropy(i)
    end do

    do var = 1,3
      do i = 2, cells+1!1, cells
        write(entropy_unit,*) ent_var(var,i)
      end do
    end do

    close(entropy_unit)

  end subroutine write_entropy

!=============================================================================80
!
! Writes solution as cell output ( fancy color plot )
! FIXME: add this... need TECPLOT cell centered format
!
!=============================================================================80
  subroutine write_soln_cells()

  end subroutine write_soln_cells

  include 'find_available_unit.f90'
  include 'conserved_to_primitive_1D.f90'

end module write_soln
