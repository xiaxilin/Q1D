subroutine triblocksolve(neq, dof, L, D, U, RHS, x)

  use set_precision, only : dp

  implicit none

  integer,                          intent(in)    :: neq, dof
  real(dp), dimension(neq,neq,dof), intent(inout) :: L, D, U
  real(dp), dimension(neq,dof),     intent(inout) :: RHS
  real(dp), dimension(neq,dof),     intent(out)   :: x

  integer :: i

  real(dp), dimension(neq,dof) :: y
  real(dp), dimension(neq,neq) :: A, AT, AT_inv, temp

  continue

  A(:,:) = D(:,:,1)
  y(:,1) = RHS(:,1)

  do i = 2, dof
    AT = transpose(A)
    call matrix_inv(neq, AT, AT_inv)

    D(:,:,i-1) = transpose(AT_inv)
    temp = transpose( matmul( AT_inv, transpose(L(:,:,i)) ) ) 

    A      = D(:,:,i) - matmul(temp, U(:,:,i-1))
    y(:,i) = RHS(:,i) - matmul(temp, y(:,i-1))
  end do

  call matrix_inv(neq, A, AT_inv)

  x(:,dof) = matmul(AT_inv, y(:,dof))

  do i = dof-1,1,-1
    x(:,i) = matmul( D(:,:,i), (y(:,i) - matmul( U(:,:,i), x(:,i+1) )) )
  end do

end subroutine triblocksolve
