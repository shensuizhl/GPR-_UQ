module AlgebraDeterm
  implicit none
contains
! 求矩阵的Determinant值
 function Determinant(matrix)
  real*8:: Determinant
  real*8    :: matrix(:,:)
  real*8, allocatable :: ma(:,:)
  integer :: i,N
  N = size(matrix,1)
  allocate(ma(N,N))
  ma = matrix
  call Upper(ma)  
  Determinant = 1.0
  do i=1,N
    Determinant = Determinant*ma(i,i)
  end do
end function
! 求上三角矩阵的子程序
subroutine Upper(matrix)
  real*8    :: matrix(:,:)
  integer :: M,N
  integer :: I,J
  real*8 :: E
  M=size(matrix,1)
  N=size(matrix,2)
  do I=1,N-1
	do J=I+1,M		
	  E=matrix(J,I)/matrix(I,I)
      ! 用90的功能可以少一层循环
	  matrix(J,I:M)=matrix(J,I:M)-matrix(I,I:M)*E
	end do
  end do
  return
end subroutine Upper
end module

