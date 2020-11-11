module AlgebraDeterm
  implicit none
contains
! ������Determinantֵ
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
! �������Ǿ�����ӳ���
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
      ! ��90�Ĺ��ܿ�����һ��ѭ��
	  matrix(J,I:M)=matrix(J,I:M)-matrix(I,I:M)*E
	end do
  end do
  return
end subroutine Upper
end module

