!================================= det_3x3 ===================================80
!
! computes determinant of a 3x3 matrix
!
!=============================================================================80

pure function det_3x3( mat )

  use set_precision,  only : dp

  real(dp), dimension(3,3), intent(in) :: mat
  real(dp)                             :: det_3x3

  continue

  det_3x3 = mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)) &
          - mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1)) &
          + mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

end function det_3x3
