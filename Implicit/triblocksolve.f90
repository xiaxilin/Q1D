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
    soln(:,i) = matmul( DD(:,:,i), (y(:,i) - matmul( UD(:,:,i), soln(:,i+1) )) )
  end do

end subroutine triblocksolve
