!#####################################################################
!                                                                    #
!  The FORTRAN programs are written by Dan Lu at the Department of   #
!  Scientific Computing at the Florida State University under        #
!  supervision of Professor Ming Ye with support of NSF-EAR grant    #
!  0911074 and DOE-SBR grant DESC0002687.                            #
!  The codes implement the DREAM algorithm developed by Jasper A.    #
!  Vrugt for Markov Chain Monte Carlo simulation.                    #
!  There is no guarantee that this program will suit the user's      #
!  needs or goals, execute efficiently and without mishap on the     #
!  user's computer, exhibit no errors or bugs, or yield a            #
!  scientifically defensible result.                                 #
!  The user is welcome to make any modifications needed to suit      #
!  his/her interpretive and modeling needs.                          #
!                                                                    #
!  Dan Lu (demi.ludan@gmail.com)                                     # 
!  Florida state university, August, 2013                            #
!                                                                    #
!#####################################################################

module dream

   use rand_generator
   use pdf_density 
   use Linearinverse
   use AlgebraDeterm
   use multinorm
   !use run_mod
   implicit none

! -- Private variables

   private par_dim,ob_dim,max_gen,Z, fit,dim_pre,in_dim,seed,api!���Ŷ�


! -- Private subroutines           
 
   !private likelihood
!*****************************************************
!           PARAMETER DEFINITION
!*****************************************************

!  Number of parameters
   integer :: par_dim,Status
   integer::dim_pre,in_dim,seed!Ԥ��ĸ���,���εĸ���(�˴�����Ϊ3),����ֵά��,�������
   real::api!���Ŷ�
!  observed data
   character(50) :: ob_file
   integer :: ob_dim
   real*8, allocatable :: ob_data(:)   
   real*8,allocatable:: yi(:,:)!yi(in_dim,ob_dim)������������ob_data(ob_dim)��Ӧ
   real*8,allocatable:: yo(:,:)!yo(in_dim,dim_pre)��Ԥ�⴦��������

!  Maximum generations for a single chain
   integer :: max_gen
 


!  Markov chain 
   real*8, allocatable :: Z(:)
   real*8, allocatable :: fit(:,:)
   !real*8, allocatable :: Zp(:)
!ģ�ͼ���ֵÿ�β���sampleʱ�õ���
   !real*8,allocatable::m_comp(:,:,:)!ģ��ֵ
   real*8,allocatable::cov_0(:,:)!Э�������
   real*8,allocatable::invercov(:,:)!Э�������������
   real*8,allocatable::cov_1(:,:)!Ԥ�����1 
   real*8,allocatable::cov_2(:,:)!Ԥ�����2
   real*8,allocatable::b_pre(:)!ģ��Ԥ��ֵ���
!  Acceptance rate

!  Define the prior density type


   contains 

 
!**************************************************************************************
!
!                       FUNCTIONS FOR INITIALIZATING CHAINS
!
!**************************************************************************************
!--------------------------------------------------------------------------------------
!***  Read input parameters from the input file


   subroutine read_vars( )

      implicit none

      integer :: iunit
      integer :: i

      call get_unit(iunit)
      
      open(iunit, file = 'DreamTest.txt' )
      read(iunit,*) 
      read(iunit,*) par_dim 
      read(iunit,*) 
      allocate(z(par_dim))
      read(iunit,*) z(1:par_dim)
      read(iunit,*) 
      read(iunit,*) max_gen 
      read(iunit,*) 
      read(iunit,'(a50)') ob_file
      read(iunit,*) 
      read(iunit,*) ob_dim
      close(iunit)

   end subroutine

!**********************************************************
! import observation data

     subroutine import_observ()

      implicit none

      integer :: iunit, i

      allocate( ob_data( ob_dim ) )

      call get_unit( iunit )
      open(iunit, file = trim(ob_file) )

      !read(iunit,*)

      do i = 1, ob_dim
          

        read(iunit,*) ob_data(i)
        !read(iunit,*)
      end do

      close(iunit)

   end subroutine
!����۲�㼰Ԥ���
subroutine import_point( )
  integer::j
  integer::status=0
  open(12,file='input_point.txt',form="formatted",status="old")!����۲��
  read(12,"(I)") in_dim
  allocate(yi(in_dim,ob_dim))!yi(in_dim,ob_dim)������������ob_data(ob_dim)��Ӧ
  do j=1,ob_dim
          read(12,*) yi(1:in_dim,j)
          !write(*,*) yi(1:in_dim,j)
  end do
!200 format(<in_dim>(f5.1,1x))  !�˴���ʽ���ļ�����
  close(12)
  
  open(20,file='predict_point.txt',form="formatted",status="old")!����Ԥ���
  read(20,*) dim_pre,seed,api!��һ����dim_pre ,seed,api
!190 format(I2,1X,I5,1X,F5.2)  !�˴���ʽ���ļ�����
  allocate(yo(in_dim,dim_pre))!yo(in_dim,dim_pre)��Ԥ�⴦��������
  do j=1,dim_pre
      read(20,*) yo(1:in_dim,j)
  end do
  close(20)

end subroutine
!--------------------------------------------------------------------------------------
!***  Start markov chains

 


!**************************************************************************************
!
!                       FUNCTIONS TO CALCULATE LOG LIKELIHOOD  
!
!**************************************************************************************

!--------------------------------------------------------------------------------------
!***  Compute log likelihood function 


  !subroutine likelihood(r,m,z)
  ! !parameter:intgegr::r num_chain;
  ! !integer::m num_gen
  !    implicit none
  !    integer::r,m,i
  !    real*8  :: z(par_dim)! notice :there are some hyperparemeters
  !    real*8  :: pi
  !    real*8::s(ob_dim)
  !    !real::invm(ob_dim,ob_dim)
  !    real*8::Determ_1,Determ_2(1,1),c1(ob_dim,1),c2(1,ob_dim),c3(1,ob_dim)
  !    call cov0(r,m,z)
  !    call inverse(cov_0(:,:),invercov(:,:))!��������\
  !    !data z/0.5,0.3,0.1,38.4,13.8,0.1869,0.0890/
  !    !call Run_Model(z, par_dim, m_comp(:,r,m),ob_dim,m_pre(:),dim_pre,yo(:,:))
  !    !call comp_mod(r,m,z)!���Է���runmodel.f90
  !    Determ_1=Determinant(cov_0(:,:))!������ʽ��ֵ
  !    !Determ_2=0
  !    pi=3.1415926
  !    do i=1,ob_dim
  !         c1(i,1)=ob_data(i)-m_comp(i)
  !    end do
  !    c2=transpose(c1)
  !    c3=matmul(c2,invercov(:,:))
  !    Determ_2=matmul(c3, c1)
  !    fit(r,m)=-1.0/2.0*Determ_2(1,1)-1.0/2.0*log(Determ_1)-ob_dim/2.0*log(2*pi)
  ! end subroutine

!--------------------------------------------------------------------------------------
!***  Use GP algorithm 
subroutine Guass()
   implicit none
   !real*8::sortpre(max_gen,dim_pre)
   !real*8::low(dim_pre),upper(dim_pre)!�ѱ����������ڣ�
   integer::i,j,info
   real*8::s_pre(dim_pre),b_pre(dim_pre),av_b(dim_pre),var_b(dim_pre,dim_pre),ave(dim_pre)
   write(*,*) z
   allocate(cov_0(ob_dim,ob_dim))
   allocate(invercov(ob_dim,ob_dim))
   allocate(cov_1(dim_pre,ob_dim))!Ԥ�����1 
   allocate(cov_2(dim_pre,dim_pre))!Ԥ�����2
   call cov0(z)
   call inverse(cov_0(:,:),invercov(:,:))!��������
   call cov1(z)
   call cov2(z)
   av_b=avbfun()
   var_b=varbfun()
   write(*,*) Determinant(var_b)
   s_pre(:)=0
   open(70,file='resultoutput.txt',status="replace")
   do j = 1, max_gen!�����ĸ���
       call multinormal_sample(dim_pre,var_b,av_b,seed,b_pre,info)!�õ��ṹ���b��dim_preά
       write(70,'(<dim_pre>(1x,F12.5))') b_pre(:)
       s_pre(:)=b_pre(:)+s_pre(:)
   end do
   close(70)
   do i=1,dim_pre
     ave(i)=s_pre(i)/(max_gen)
     !write(*,*) 'ave=',ave(i)
   end do
end subroutine


subroutine cov0(z)
implicit none
integer::r,k,i,j,m
real*8::s
real*8::z(par_dim)!�˴�����Ϊ�������Σ����ڸ�������Ҫ��߳���������
!data lens /1000.0,500.0/
do i=1,ob_dim
        do j=1,ob_dim
            s=0.0
            do k=1,in_dim
                s=((yi(k,i)-yi(k,j)))**2+s
            end do
            !WRITE(*,*) s
            
            s=s/(2*z(1)*z(1))!yiΪ�����������ͨ�������ļ���ȡ
            if (i==j) then
                cov_0(i,j)=z(3)*z(3)*exp(-s)+z(2)*z(2)
            else
                cov_0(i,j)=z(3)*z(3)*exp(-s)
            end if
        end do
        !STOP
end do
end subroutine
subroutine cov1(z)
implicit none
integer::r,k,i,j,m
real*8::s
real*8::z(par_dim)
do i=1,dim_pre
        do j=1,ob_dim
            s=0
            do k=1,in_dim
                s=((yo(k,i)-yi(k,j)))**2+s
            end do
            s=s/(2*z(1)*z(1))!yiΪ�����������ͨ�������ļ���ȡ
            cov_1(i,j)=z(3)*z(3)*exp(-s)
        end do
    end do
end subroutine
subroutine cov2(z)
implicit none
integer::r,k,i,j,m
real*8::s
real*8::z(par_dim)!�˴�����Ϊ��ά���Σ����ڸ�������Ҫ��߳���������
do i=1,dim_pre
        do j=1,dim_pre
            s=0
            do k=1,in_dim
                s=((yo(k,i)-yo(k,j)))**2+s
            end do
            s=s/(2*z(1)*z(1))!yiΪ�����������ͨ�������ļ���ȡ
            cov_2(i,j)=z(3)*z(3)*exp(-s)
        end do
    end do
end subroutine

 !����ṹ���ƽ��ֵ
 function avbfun( )
 real*8::avbfun(dim_pre)!dim_preά
 real*8::s(ob_dim),c1(ob_dim,1),c2(dim_pre,ob_dim),c3(dim_pre,1)!ob_dim
 integer::r,m
 integer::i,j,k
 do i=1,ob_dim
    c1(i,1)=ob_data(i)
 end do
 c2=matmul(cov_1(:,:),invercov(:,:))
 c3=matmul(c2,c1)
 do i=1,dim_pre
     avbfun(i)=c3(i,1)
 end do
 end function
 
 !����ṹ���Э����
 
  function varbfun()
 real*8::varbfun(dim_pre,dim_pre),c1(dim_pre,ob_dim),c2(dim_pre,dim_pre)!dim_pre*dim_preά
 real*8::s(ob_dim)
 integer::i,j,k,n,r,m
 !varbfun
 c1=matmul(cov_1(:,:),invercov(:,:))!c1(dim_pre,ob_dim)
 c2=matmul(c1,transpose(cov_1(:,:)))!c2(dim_pre,dim_pre)
 varbfun(:,:)=cov_2(:,:)-c2(:,:)
 do i=1,dim_pre
     write(*,*) varbfun(i,i)
 end do
 end function
          
!function normal(mean,sigma) !������̫�����
!    implicit none 
!    integer :: flag 
!    double precision, parameter :: pi = 3.141592653589793239  
!    real :: u1, u2, y1, y2, normal, mean, sigma 
!    save flag 
!    data flag /0/ 
!    u1 = ran(); u2 = ran() 
!    if (flag.eq.0) then 
!      y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2) 
!      normal = mean + sigma*y1 
!      flag = 1 
!    else 
!      y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2) 
!      normal = mean + sigma*y2 
!      flag = 0 
!    endif  
! end function normal 

subroutine bubble_sort(a,n)
  implicit none
  integer :: n
  real*8::a(n),temp
  integer i,j
  do i=n-1,1,-1   ! ��ʼ��n-1�ε�ɨ��
    do j=1,i      ! һ��һ�Ե����Ƚϣ�i֮������ֲ��ñȽ�
    ! ���a(j) > a(j+1) �Ͱ���������ֵ����
      if ( a(j) > a(j+1) ) then
        temp=a(j)
        a(j)=a(j+1)
        a(j+1)=temp
      end if
    end do
  end do
  return
end subroutine
  subroutine get_unit ( iunit )

       implicit none

       integer ( kind = 4 ) i
       integer ( kind = 4 ) ios
       integer ( kind = 4 ) iunit
       logical lopen

       iunit = 0

       do i = 1, 99

          if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

             inquire ( unit = i, opened = lopen, iostat = ios )

             if ( ios == 0 ) then
                if ( .not. lopen ) then
                   iunit = i
                   return
                end if
             end if

          end if

       end do

       return

   end subroutine
end module