module matrix_manip

  implicit none

  private

  public :: triblocksolve
  public :: pentablocksolve
  public :: ludcmp
  public :: lubkslv
  public :: matrix_inv
  public :: mat_inv_3x3

!  public :: det_3x3
!  public :: det_4x4
!  public :: det_5x5

contains

!============================== triblocksolve ================================80
!
! Uses the Thomas algorithm to solve a generic block tridiagonal system
! Destroys the DD, UD, RHS matrices
!
!=============================================================================80
  subroutine triblocksolve(neq, dof, LD, DD, UD, RHS, soln)

    use set_precision, only : dp

    implicit none

    integer,                          intent(in)    :: neq, dof
    real(dp), dimension(neq,neq,dof), intent(inout) :: LD, DD, UD
    real(dp), dimension(neq,dof),     intent(inout) :: RHS
    real(dp), dimension(neq,dof),     intent(out)   :: soln

    integer                      :: i
    real(dp), dimension(neq,neq) :: temp

    continue

! Normalize the first row...
!    call matrix_inv(neq, DD(:,:,1), temp)
    call mat_inv_3x3(DD(:,:,1), temp)

    DD(:,:,1) = matmul(temp,DD(:,:,1))
    UD(:,:,1) = matmul(temp,UD(:,:,1))
    RHS(:,1)  = matmul(temp,RHS(:,1))

! Loop to eliminate lower diagonal and then normalize the row
    do i = 2, dof
      DD(:,:,i) = DD(:,:,i) - matmul(LD(:,:,i), UD(:,:,i-1))
      RHS(:,i)  = RHS(:,i)  - matmul(LD(:,:,i), RHS(:,i-1))

!      call matrix_inv(neq, DD(:,:,i), temp)
      call mat_inv_3x3(DD(:,:,i), temp)

      DD(:,:,i) = matmul(temp,DD(:,:,i))
      UD(:,:,i) = matmul(temp,UD(:,:,i))
      RHS(:,i)  = matmul(temp,RHS(:,i))
    end do

! Back solve... since the diagonal is an identity matrix this is easy
    soln(:,dof) = RHS(:,dof)

    do i = dof-1,1,-1
      soln(:,i) = RHS(:,i) - matmul(UD(:,:,i), soln(:,i+1))
    end do

  end subroutine triblocksolve

!============================== pentablocksolve ==============================80
!
! Uses a modified Thomas algorithm to solve a generic block pentadiagonal system
!
!=============================================================================80
  subroutine pentablocksolve(neq, dof, LD2, LD, DD, UD, UD2, RHS, soln)

    use set_precision, only : dp

    implicit none

    integer,                          intent(in)    :: neq, dof
    real(dp), dimension(neq,neq,dof), intent(inout) :: LD2, LD, DD, UD, UD2
    real(dp), dimension(neq,dof),     intent(inout) :: RHS
    real(dp), dimension(neq,dof),     intent(out)   :: soln

    integer                      :: i
    real(dp), dimension(neq,neq) :: temp, temp2

    continue

! Eliminate the LD2 diagonal
    do i = 3, dof-1
      call mat_inv_3x3(LD(:,:,i-1), temp)
      temp2 = matmul(LD2(:,:,i), temp)

      LD(:,:,i) = LD(:,:,i) - matmul(temp2, DD(:,:,i-1))
      DD(:,:,i) = DD(:,:,i) - matmul(temp2, UD(:,:,i-1))
      UD(:,:,i) = UD(:,:,i) - matmul(temp2, UD2(:,:,i-1))
      RHS(:,i)  = RHS(:,i)  - matmul(temp2, RHS(:,i-1))
    end do

    i = dof
    call mat_inv_3x3(LD(:,:,i-1), temp)
    temp2 = matmul(LD2(:,:,i), temp)

    LD(:,:,i) = LD(:,:,i) - matmul(temp2, DD(:,:,i-1))
    DD(:,:,i) = DD(:,:,i) - matmul(temp2, UD(:,:,i-1))
    RHS(:,i)  = RHS(:,i)  - matmul(temp2, RHS(:,i-1))

! Eliminate the LD diagonal
    do i = 2, dof-1
      call mat_inv_3x3(DD(:,:,i-1), temp)
      temp2 = matmul(LD(:,:,i), temp )

      DD(:,:,i) = DD(:,:,i) - matmul(temp2, UD(:,:,i-1))
      UD(:,:,i) = UD(:,:,i) - matmul(temp2, UD2(:,:,i-1))
      RHS(:,i)  = RHS(:,i)  - matmul(temp2, RHS(:,i-1))

    end do

    i = dof
    call mat_inv_3x3(DD(:,:,i-1), temp)
    temp2 = matmul(LD(:,:,i), temp )

    DD(:,:,i) = DD(:,:,i) - matmul(temp2, UD(:,:,i-1))
    RHS(:,i)  = RHS(:,i)  - matmul(temp2, RHS(:,i-1))

! Back solve
    call mat_inv_3x3(DD(:,:,dof), temp)
    soln(:,dof)   = matmul(temp, RHS(:,dof))

    call mat_inv_3x3(DD(:,:,dof-1), temp)

    soln(:,dof-1) = matmul(temp, &
      RHS(:,dof-1) - matmul(UD(:,:,dof-1), soln(:,dof)))

    do i = dof-2,1,-1
      call mat_inv_3x3(DD(:,:,i), temp)     
      soln(:,i) = matmul(temp, RHS(:,i)                                        &
        - matmul(UD(:,:,i), soln(:,i+1)) - matmul(UD2(:,:,i), soln(:,i+2)))
    end do

  end subroutine pentablocksolve

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
!=============================== mat_inv_3x3 =================================80
!
! Calculates inverse of 3x3 matrix for speed!!!!
!
!=============================================================================80

  subroutine mat_inv_3x3(mat, inv)

    use set_precision, only : dp
    use set_constants, only : zero, one

    implicit none

    real(dp), dimension(3,3), intent(in)  :: mat
    real(dp), dimension(3,3), intent(out) :: inv

    real(dp) :: det

    continue

    inv(1,1) = mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
    inv(2,1) = mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
    inv(3,1) = mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)

    inv(1,2) = mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
    inv(2,2) = mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
    inv(3,2) = mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2)

    inv(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
    inv(2,3) = mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)
    inv(3,3) = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

    det = det_3x3(mat)

    if (det == zero) det = one

    inv = inv/det

  end subroutine mat_inv_3x3

  include 'det_3x3.f90'

end module matrix_manip
