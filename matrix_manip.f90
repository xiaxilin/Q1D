module matrix_manip

  implicit none

  private

  public :: triblocksolve
  public :: ludcmp
  public :: lubkslv
  public :: matrix_inv

!  public :: det_3x3
!  public :: det_4x4
!  public :: det_5x5

contains

!=============================== triblocksolve ===============================80
!
! Uses the Thomas algorithm to solve a generic block tridiagonal system
!
!=============================================================================80
  subroutine triblocksolve(neq, dof, LD, DD, UD, RHS, soln)

    use set_precision, only : dp

    implicit none

    integer,                          intent(in)    :: neq, dof
    real(dp), dimension(neq,neq,dof), intent(inout) :: LD, DD, UD
    real(dp), dimension(neq,dof),     intent(inout) :: RHS
    real(dp), dimension(neq,dof),     intent(out)   :: soln

    integer :: i

    real(dp), dimension(neq,dof) :: y
    real(dp), dimension(neq,neq) :: A, AT, AT_inv, temp

    continue

    A(:,:) = DD(:,:,1)
    y(:,1) = RHS(:,1)

    do i = 2, dof
      AT = transpose(A)
      call matrix_inv(neq, AT, AT_inv)

      DD(:,:,i-1) = transpose(AT_inv)
      temp = transpose( matmul( AT_inv, transpose(LD(:,:,i)) ) ) 

      A      = DD(:,:,i) - matmul(temp, UD(:,:,i-1))
      y(:,i) = RHS(:,i)  - matmul(temp, y(:,i-1))
    end do

    call matrix_inv(neq, A, AT_inv)

    soln(:,dof) = matmul(AT_inv, y(:,dof))

    do i = dof-1,1,-1
      soln(:,i) = matmul( DD(:,:,i), (y(:,i) - matmul(UD(:,:,i), soln(:,i+1))) )
    end do

  end subroutine triblocksolve

!================================= ludcmp ====================================80
!
! This subroutine performs LU decomposition without pivoting,
! most ( CFD related ) matrices passed to it will be diagonally dominant
!
!=============================================================================80
  subroutine ludcmp(neq, matrix, lower, upper)

    use set_precision, only : dp
    use set_constants, only : zero, one

    implicit none

    integer,                      intent(in)  :: neq
    real(dp), dimension(neq,neq), intent(in)  :: matrix
    real(dp), dimension(neq,neq), intent(out) :: upper, lower

    integer  :: i, j, k
    real(dp) :: factor

    continue

    lower = zero
    upper = zero

! set up diagonal of lower matrix
    do i = 1, neq
      lower(i,i) = one
    end do

! set first row of upper and first column of lower matrices
    upper(1,:) = matrix(1,:)
    lower(:,1) = matrix(:,1)/upper(1,1)

    do i = 2, neq-1
! set diagonal of upper matrix
      factor = zero
      do k = 1, i-1
        factor = factor + lower(i,k)*upper(k,i)
      end do
      upper(i,i) = matrix(i,i) - factor

! off diagonals of upper matrix
      do j = i+1, neq
        factor = zero
        do k = 1, i-1
          factor = factor + lower(i,k)*upper(k,j)
        end do
        upper(i,j) = matrix(i,j) - factor

! off diagonals of lower matrix
        factor = zero
        do k = 1, i-1
          factor = factor + lower(j,k)*upper(k,i)
        end do
        lower(j,i) = (matrix(j,i) - factor)/upper(i,i)
      end do
    end do

    factor = zero
    do k = 1, neq-1
      factor = factor + lower(neq,k)*upper(k,neq)
    end do

    upper(neq,neq) = matrix(neq,neq) - factor

  end subroutine ludcmp

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine lubkslv(neq, lower, upper, b, x)

    use set_precision, only : dp
    use set_constants, only : zero

    implicit none

    integer,                      intent(in)  :: neq
    real(dp), dimension(neq,neq), intent(in)  :: lower, upper
    real(dp), dimension(neq),     intent(in)  :: b
    real(dp), dimension(neq),     intent(out) :: x

    integer  :: i, j
    real(dp) :: factor
    real(dp), dimension(neq) :: y

    continue

! forward substitution
    y(1) = b(1)
    do i = 2, neq
      factor = zero
      do j = 1, i-1
        factor = factor + lower(i,j)*y(j)
      end do
      y(i) = b(i) - factor
    end do

! back substitution
    x(neq) = y(neq)/upper(neq,neq)
    do i = neq-1,1,-1
      factor = zero
      do j = i+1, neq
        factor = factor + upper(i,j)*x(j)
      end do
      x(i) = (y(i) - factor)/upper(i,i)
    end do

  end subroutine lubkslv

!=============================================================================80
!
!
!
!=============================================================================80

  subroutine matrix_inv(neq, matrix, inv)

    use set_precision, only : dp
    use set_constants, only : zero, one

    implicit none

    integer,                      intent(in)  :: neq
    real(dp), dimension(neq,neq), intent(in)  :: matrix
    real(dp), dimension(neq,neq), intent(out) :: inv

    integer :: i
    real(dp), dimension(neq) :: b, soln
    real(dp), dimension(neq,neq) :: l,u

    continue

! perform LU decomposition
    call ludcmp(neq, matrix, l, u)

! perfrom backsolve to get inverse
    do i = 1, neq
      b    = zero
      b(i) = one
      call lubkslv(neq, l, u, b, soln)
      inv(:,i) = soln
    end do

  end subroutine matrix_inv

end module matrix_manip
