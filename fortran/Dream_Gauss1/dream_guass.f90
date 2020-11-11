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
   implicit none

! -- Private variables

   private par_dim,  range, ob_dim,                                                &
           nchain, max_gen, npair, jumpstep,  restart_flag,                        &
           nCR, CR, CR_ind, pCR, L_CR, Dis_CR,                                     &
           Gel_Ru_R, printstep, GR_count, conv_flag, GR_threshold,                 &
           jumprate_table, jumprate, jump_num, jump_dim,                           &
           Z, fit, Zp, AR,                                                         &   
           mnor_dim, uni_dim, mnor_ind, uni_ind, mnor_mean, cov_mat, parm, logdet, &
           prior_den,dim_pre,in_dim,seed,api!置信度


! -- Private subroutines           
 
   private likelihood,       & 
           get_prior_sample, &
           prior_density,    &
           gen_candidate,    &
           outofbound,       &
           init_CR,          &
           choose_CR,        &
           update_CR_dis,    &
           update_pCR,       &
           comp_std,         &
           choose_jumprate,  &
           comp_diff,        &
           Gel_Ru,           &
           Output,           &
           get_unit
           
           

!*****************************************************
!           PARAMETER DEFINITION
!*****************************************************

!  Number of parameters
   integer :: par_dim
   integer::dim_pre,in_dim,seed!预测的个数,超参的个数(此处设置为3),输入值维数,随机种子
   real::api!置信度
!  Range of each parameter
   real  :: range( 2, 1000 )

!  observed data
   character(50) :: ob_file
   integer :: ob_dim
   real*8, allocatable :: ob_data(:)   
   real,allocatable:: yi(:,:)!yi(in_dim,ob_dim)即输入坐标与ob_data(ob_dim)对应
   real,allocatable:: yo(:,:)!yo(in_dim,dim_pre)即预测处输入坐标
!  Number of Markov chains
   integer :: nchain   

!  Maximum generations for a single chain
   integer :: max_gen

!  Number of pairs of chains to generate candidate samples
   integer :: npair

!  Number of steps to reset jump rate Gamma
   integer :: jumpstep

!  Restart the chain from last generation of last time simulation or not
   character(10) :: restart_flag
   

!*** Other variables used in this module
!  Crossover parameters
   integer :: nCR
   real  :: CR(20)
   integer :: CR_ind
   real    :: pCR(20) 
   integer :: L_CR(20)
   real  :: Dis_CR(20)

!  Gelman Rubin R statistic used to check MC convergence
   real, allocatable :: Gel_Ru_R(:,:)
   integer :: printstep
   integer :: GR_count
   logical :: conv_flag
   real  :: GR_threshold 

!  Jump rate table
   real, allocatable :: jumprate_table(:)
   real  :: jumprate
   integer :: jump_num 
   integer, allocatable :: jump_dim(:)

!  Markov chain 
   real, allocatable :: Z(:,:,:)
   real, allocatable :: fit(:,:)
   real, allocatable :: Zp(:)
!模型计算值每次参数sample时得到的
   real,allocatable::m_comp(:,:,:)!模拟值
   real,allocatable::cov_0(:,:,:,:)!协方差矩阵
   real,allocatable::invercov(:,:,:,:)!协方差矩阵的逆矩阵
   real,allocatable::cov_1(:,:,:,:)!预测矩阵1 
   real,allocatable::cov_2(:,:,:,:)!预测矩阵2
   real,allocatable::m_pre(:,:,:)!模型预测值不加结构误差
   real,allocatable::pred(:,:,:)!预测值加结构误差
   real,allocatable::ave(:)!预测平均值
   real,allocatable::m_ave(:)
!  Acceptance rate
   real :: AR

!  Define the prior density type
   type density
      character*40 :: type                    ! prior density type
      integer      :: numrealpar              ! number of density parameters
      real         :: realpar(2)              ! array to save density parameters 
   end type density

!  Variables for multinormal density
   integer :: mnor_dim, uni_dim
   integer :: mnor_ind(1000), uni_ind(1000)
   real, allocatable :: mnor_mean(:)          ! mean vector
   real, allocatable :: cov_mat(:,:)          ! covariance matrix
   real, allocatable :: parm(:)               ! used for calculate inverse covariance matrix of multi-normal distribution 
   real :: logdet                             ! log determinant of covariance matrix

   type( density ), allocatable :: prior_den(:)


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
      read(iunit,*) range(1,1:par_dim) 
      read(iunit,*) 
      read(iunit,*) range(2,1:par_dim)  
      read(iunit,*) 
      read(iunit,*) max_gen 
      !write(*,*)max_gen
      read(iunit,*) 
      read(iunit,*) nchain
      read(iunit,*) 
      read(iunit,'(a50)') ob_file
      read(iunit,*) 
      read(iunit,*) ob_dim
      read(iunit,*) 
      read(iunit,*) npair 
      read(iunit,*) 
      read(iunit,*) nCR 
      read(iunit,*) 
      read(iunit,*) jumpstep 
      read(iunit,*) 
      read(iunit,*) printstep 
      read(iunit,*) 
      read(iunit,*) GR_threshold 
      read(iunit,*) 
      read(iunit,*) restart_flag
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
!导入观测点及预测点
subroutine import_point( )
  integer::j
  integer::status=0
  open(12,file='input_point.txt',form="formatted",status="old")!导入观测点
  read(12,"(I)") in_dim
  allocate(yi(in_dim,ob_dim))!yi(in_dim,ob_dim)即输入坐标与ob_data(ob_dim)对应
  do j=1,ob_dim
          read(12,200) yi(1:in_dim,j)
          !write(*,200) yi(1:10,j)
  end do
200 format(<in_dim>(f5.1,1x))  !此处格式视文件而定
  close(12)
  
  open(20,file='predict_point.txt',form="formatted",status="old")!导入预测点
  read(20,*) dim_pre,seed,api!第一行是dim_pre ,seed,api
190 format(I2,1X,I5,1X,F5.2)  !此处格式视文件而定
  allocate(yo(in_dim,dim_pre))!yo(in_dim,dim_pre)即预测处输入坐标
  do j=1,dim_pre
      read(20,200) yo(:,j)
      !write(*,*) yo(:,j)
  end do
  close(20)

end subroutine
!--------------------------------------------------------------------------------------
!***  Start markov chains

   subroutine init_chain( )

      implicit none

      integer :: i, j, iunit, dummyi

      allocate( Z(par_dim, nchain, max_gen))
      allocate(m_comp(ob_dim,nchain, max_gen))!模拟输出值
      allocate(cov_0(ob_dim,ob_dim,nchain, max_gen))!协方差矩阵
      allocate(invercov(ob_dim,ob_dim,nchain, max_gen))!协方差矩阵的逆矩阵
      allocate(cov_1(dim_pre,ob_dim,nchain, max_gen))!预测矩阵1 
      allocate(cov_2(dim_pre,dim_pre,nchain,max_gen))!预测矩阵2
      allocate( fit( nchain, max_gen ))!似然值
      allocate( m_pre(dim_pre,nchain, max_gen ))!模型预测值
      allocate( pred(dim_pre,nchain, max_gen ))!预测值
      allocate(ave(dim_pre))!平均值
      allocate(m_ave(dim_pre))!模拟平均值
      allocate( Zp( par_dim)) 
    ! if restart the chain, continue the chain from the last sample of last time simulation
       !write(*,*) 'please enter a seed number(0~12345678):'
       !read(*,*) seed
      if( trim(restart_flag) == 'yes' )then

         call get_unit( iunit )
         open(iunit, file = 'Restart.out' )
         read(iunit,*) 
         
         do i = 1, nchain
            read(iunit, *)dummyi, fit(i,1), Z(:,i,1)
         end do

         close(iunit)

    ! initialize the chain from a centrain prior distribution
      
      else

         do i =1, nchain
            Z(:,i,1) = get_prior_sample( )
            !m_comp(:,i,1)=cm(Z(:,i,1))!计算初始模拟值模型计算值M(:,x,t)(1…C)
            !call cov0(i,1,Z)
            !call cov1(i,1,Z)
            !call cov2(i,1,Z)
            !fit(i,1) = likelihood(Z(:,i,1))
            call likelihood(i,1,Z(:,i,1))
         end do

      end if

   end subroutine


!**************************************************************************************
!
!                       FUNCTIONS TO CALCULATE LOG LIKELIHOOD  
!
!**************************************************************************************

!--------------------------------------------------------------------------------------
!***  Compute log likelihood function 


   subroutine likelihood(r,m,z)
   !parameter:intgegr::r num_chain;
   !integer::m num_gen
      implicit none
      integer::r,m,i,j
      real  :: z(par_dim)! notice :there are some hyperparemeters
      real  :: pi
      real::s(ob_dim)
      !real::invm(ob_dim,ob_dim)
      real::Determ_1,Determ_2(1,1),c1(ob_dim,1),c2(1,ob_dim),c3(1,ob_dim)
      call cov0(r,m,z)
      call cov1(r,m,z)
      call cov2(r,m,z)
      call inverse(cov_0(:,:,r,m),invercov(:,:,r,m))!求矩阵的逆
      call comp_mod(r,m,z)
      Determ_1=Determinant(cov_0(:,:,r,m))!求行列式的值
      Determ_2=0
      pi=3.1415926
      do i=1,ob_dim
           c1(i,1)=ob_data(i)-m_comp(i,r,m)
      end do
      c2=transpose(c1)
      c3=matmul(c2,invercov(:,:,r,m))
      Determ_2=matmul(c3, c1)
      fit(r,m)=-1.0/2.0*Determ_2(1,1)-1.0/2.0*log(Determ_1)-ob_dim/2.0*log(2*pi)
     
      !-------- TEST --------------------
      !case 1: 10-D  bimodal distribution对于贝叶斯校准，只需对下面的likelihood计算公式进行计算

      !likelihood = log(1.0/3.0/(2.0*pi)*exp(-0.5d0*sum((Zx+5.0d0)**2)) + 2.0/3.0/(2.0*pi)*exp(-0.5d0*sum((Zx-5.0d0)**2)) )

   end subroutine


!**************************************************************************************
!
!                       FUNCTIONS ABOUT PRIOR INFORMATION 
!
!**************************************************************************************

!--------------------------------------------------------------------------------------
!***  Read prior information 


  subroutine init_prior( )

      implicit none
      
      integer :: iunit, i, cov_dim      
      allocate( prior_den( par_dim ) )
   
    ! number of parameters have multinormal distribution
      mnor_dim = 0
      
    ! number of parameters have univariate distribution  
      uni_dim = 0
    
      call get_unit(iunit) 

    ! In the file specify the prior distribution    
      open(iunit, file = 'prior.in' )

      read(iunit,*) 
      do i = 1, par_dim
         read(iunit, '(a20)') prior_den(i)%type
         read(iunit,*) prior_den(i)%numrealpar, prior_den(i)%realpar(1:prior_den(i)%numrealpar) 

         if( trim(prior_den(i)%type) == 'multinormal') then
            mnor_dim = mnor_dim + 1
            mnor_ind(mnor_dim) = i
         else
            uni_dim = uni_dim + 1
            uni_ind(uni_dim) = i
         end if
      end do 
      
      close(iunit)


      if( mnor_dim == 1 )then
         write(*,*) '*** ERROR: Only one parameter follows multinormal distribution. please use normal instead!'
         stop

      elseif( mnor_dim > 1 )then

         allocate( mnor_mean( mnor_dim ) )
         allocate( cov_mat(mnor_dim, mnor_dim) )
         allocate( parm( mnor_dim*(mnor_dim+3)/2+1 ) )
       
       ! read the covariance matrix of the multinormal density from matrix_file input block
         call get_unit(iunit) 
         
         open(iunit, file = 'covmatrix.dat' )
         
         read(iunit,*) 
         read(iunit,*) cov_dim

         if( cov_dim /= mnor_dim )then
	    write(*,*) '*** ERROR: Covariance Matrix dimension in matrix file must be the same with  &
	             the number of parameters having multinormal density type in MCMC_Prior_PDF block!'
            stop
         end if

         do i = 1, mnor_dim 
            read(iunit,*) cov_mat(i,:)
         end do
         
         close(iunit)

       ! collect the mean values of the multinormal desity
         do i = 1, mnor_dim
            mnor_mean(i) = prior_den(mnor_ind(i))%realpar(1)
         end do

         call setgmn( mnor_dim, mnor_mean, cov_mat, parm, logdet)
         
      end if


   end subroutine


!--------------------------------------------------------------------------------------
!***  Get a Sample from the specified prior distribution to start MCMC

   
   function get_prior_sample( )

      implicit none

      real :: get_prior_sample( par_dim )
      integer :: i
      real :: work(mnor_dim)
      real :: tmp(mnor_dim)


    ! generate samples follow univariate density    
      if( uni_dim > 0 )then  
         do i = 1, uni_dim
   
            select case( trim(prior_den(uni_ind(i))%type) )

               case('uniform')
                  get_prior_sample(uni_ind(i)) = genunf( prior_den(uni_ind(i))%realpar(1), &
                                                         prior_den(uni_ind(i))%realpar(2) ) 
         
               case('normal')
                  get_prior_sample(uni_ind(i)) = gennor( prior_den(uni_ind(i))%realpar(1), &
                                                         prior_den(uni_ind(i))%realpar(2) )

               case('beta')
                  get_prior_sample(uni_ind(i)) = genbet( prior_den(uni_ind(i))%realpar(1), &
                                                         prior_den(uni_ind(i))%realpar(2) )

               case('chi-square')
                  get_prior_sample(uni_ind(i)) = genchi( prior_den(uni_ind(i))%realpar(1) )
      
               case('inv-chi-square')
                  get_prior_sample(uni_ind(i)) = gengam( 0.5, 0.5*prior_den(uni_ind(i))%realpar(1) )
                  if( get_prior_sample(uni_ind(i)) .ne. 0.0 ) then
                     get_prior_sample(uni_ind(i)) = 1.0 / get_prior_sample(uni_ind(i))
                  end if
   
               case('scaled-inv-chi-square')
                  get_prior_sample(uni_ind(i)) = gengam( 0.5*prior_den(uni_ind(i))%realpar(1)  & 
                                                            *prior_den(uni_ind(i))%realpar(2), &
                                                         0.5*prior_den(uni_ind(i))%realpar(1) )

                  if( get_prior_sample(uni_ind(i)) .ne. 0.0 ) then
                     get_prior_sample(uni_ind(i)) = 1.0 / get_prior_sample(uni_ind(i))
                  end if

               case('gamma')
                  get_prior_sample(uni_ind(i)) = gengam( prior_den(uni_ind(i))%realpar(1), &
                                                         prior_den(uni_ind(i))%realpar(2) )
               
               case('inv-gamma')
                  get_prior_sample(uni_ind(i)) = gengam( prior_den(uni_ind(i))%realpar(1), &
                                                         prior_den(uni_ind(i))%realpar(2) )
                  if( get_prior_sample(uni_ind(i)) .ne. 0.0 ) then
                     get_prior_sample(uni_ind(i)) = 1.0 / get_prior_sample(uni_ind(i))
                  end if

               case('exponential')
                  get_prior_sample(uni_ind(i)) = genexp( prior_den(uni_ind(i))%realpar(1) )

               case default
                  write(*,*)'*** ERROR: unknown density type!' 
                  stop 

            end select
      
         end do

      end if

      
    ! generate samples follow multinormal density
      if( mnor_dim > 1 )then
      
         call genmn( parm, tmp, work )
         
         do i = 1, mnor_dim
            get_prior_sample(mnor_ind(i)) = tmp(i)
         end do
         
      end if

   end function


!--------------------------------------------------------------------------------------
!***  Compute density values of the specified prior density type to calculate Metropolis ratio 


   function prior_density( rval )

      implicit none

      real :: rval( par_dim )
      real :: prior_density
      integer :: i
      real :: tmp( mnor_dim )

      prior_density = 1.0
      

    ! calculate density values for univariate density     
      if( uni_dim > 0 )then
         do i = 1, uni_dim
   
            select case( trim(prior_den(uni_ind(i))%type) )
   
               case('uniform')
                  prior_density =   prior_density &
                                  * getpdfunf( prior_den(uni_ind(i))%realpar(1), &
                                               prior_den(uni_ind(i))%realpar(2), rval(uni_ind(i)) ) 
         
               case('normal')
                  prior_density =   prior_density &
                                  * getpdfnor( prior_den(uni_ind(i))%realpar(1), &
                                               prior_den(uni_ind(i))%realpar(2), rval(uni_ind(i)) ) 
   
               case('beta')
                  prior_density =   prior_density &
                                  * getpdfbet( prior_den(uni_ind(i))%realpar(1), &
                                               prior_den(uni_ind(i))%realpar(2), rval(uni_ind(i)) ) 
   
               case('chi-square')
                  prior_density =   prior_density &
                                  * getpdfchi( prior_den(uni_ind(i))%realpar(1), rval(uni_ind(i)) ) 
            
               case('inv-chi-square')
                  prior_density =   prior_density &
                                  * getpdfinvchi( prior_den(uni_ind(i))%realpar(1), rval(uni_ind(i)) ) 
   
               case('scaled-inv-chi-square')
                  prior_density =   prior_density &
                                  * getpdfscinvchi( prior_den(uni_ind(i))%realpar(1), &
                                                    prior_den(uni_ind(i))%realpar(2), rval(uni_ind(i)) ) 
   
               case('gamma')
                  prior_density =   prior_density &
                                  * getpdfgam( prior_den(uni_ind(i))%realpar(1), &
                                               prior_den(uni_ind(i))%realpar(2), rval(uni_ind(i)) ) 
            
               case('inv-gamma')
                  prior_density =   prior_density &
                                  * getpdfinvgam( prior_den(uni_ind(i))%realpar(1), &
                                                  prior_den(uni_ind(i))%realpar(2), rval(uni_ind(i)) ) 

               case('exponential')
                  prior_density =   prior_density &
                                  * getpdfexp( prior_den(uni_ind(i))%realpar(1), rval(uni_ind(i)) ) 
  
            end select

         end do

      end if
      

    ! calculate density values for multinormal density     
      if( mnor_dim > 1 )then
      
         do i = 1, mnor_dim
            tmp(i) = rval(mnor_ind(i))
         end do
         
         prior_density = getpdfmnor( mnor_dim, mnor_mean(1:mnor_dim), cov_mat, logdet, tmp ) 
         
      end if

   end function


!**************************************************************************************
!
!                       FUNCTIONS ABOUT DREAM ALGORITHM  
!
!**************************************************************************************

!--------------------------------------------------------------------------------------
!***  Use DREAM algorithm to get a candicate parameter sample 
!and use guass process algorithm coupling with dream.


   subroutine DREAM_Guass( )
     
     implicit none
     
     integer :: i, j, accept, Zp_count,k,ai,bi!ai,bi上下限的项数
     real  :: ratio
     real::sortpre(max_gen*nchain,dim_pre)
     real::low(dim_pre),upper(dim_pre)!把变量声明放在？
     real::s_pre(dim_pre),sm_pre(dim_pre)
     Zp_count = 0
     accept = 0
     s_pre(:)=0
     sm_pre(:)=0

   ! Initialize the CR values
     call init_CR()
    open(40,file='sample.txt',status="replace")
     do i = 2, max_gen
        do j = 1, nchain

         ! Choose a CR value
           call choose_CR( )
	   !write(*,*)'hello'
         ! Generate a candidate
           call gen_candidate( i, j )
           write(*,*) Zp
           Zp_count = Zp_count + 1
         ! Compute log likelihood function 
           call likelihood(j,i,Zp)
         !fit(j,i) = likelihood( Zp )
	
         ! Compute the metroplis ratio 

            ratio = min( exp( ( fit(j,i) + log( prior_density(Zp) ) ) - ( fit(j,i-1) + log( prior_density(Z(:,j,i-1)) ) ) ), 1.0 )
       !跟新样本及矩阵
	   if( ratio >= rand_uni01() )then
	      Z(:,j,i) = Zp
	      accept = accept + 1
	   else
	      Z(:,j,i) = Z(:,j,i-1)
	      fit(j,i) = fit(j,i-1)
          cov_0(:,:,j,i)=cov_0(:,:,j,i-1)
          cov_1(:,:,j,i)=cov_1(:,:,j,i-1)
          cov_2(:,:,j,i)=cov_2(:,:,j,i-1)
          invercov(:,:,j,i)=invercov(:,:,j,i-1)
          m_comp(:,j,i)=m_comp(:,j,i-1)
          !更新矩阵否？
       end if
       write(*,'(a9,I4,I4,A3)') 'sample(',j,i,')'
      
       write(40,'(<par_dim>F8.5)') Z(:,j,i)
       
      !此处可对所得样本模型值进行计算
        call comp_pre(j,i)
        call predict(j,i)
        !write(*,*) 'pridiction is'
        !write(*,*) m_pre(:,j,i)
        !write(*,*) m_comp(:,j,i)
        !计算总和
        s_pre(:)=pred(:,j,i)+s_pre(:)
        sm_pre(:)=m_pre(:,j,i)+sm_pre(:)
        
        !write(*,*) 'pridiction is'
        !write(*,*) pred(:,j,i),s_pre(:)
        !得到预测值序列,调换行列减少后面排序运算量
        do k=1,dim_pre
            sortpre(j-nchain+i*nchain,k)=pred(k,j,i)
            !write(*,*) 'pred(k,j,i)'
            !write(*,*) sortpre(j-nchain+i*nchain,k)
        end do
         ! Update CR distance
           if( conv_flag .eqv. .false. .and. nCR > 1 )then  
              call update_CR_dis( i, j )
           end if
        end do

      ! Update the multinomial distribution of CR
        if( conv_flag .eqv. .false. .and. nCR > 1 .and. mod(i,10)==0 ) then
           call update_pCR( i )
        end if

      ! Compute Gelman Rubin R and export result
        if( mod(i, printstep) == 0 )then
           call Gel_Ru(i)
           call Output_1( i-printstep+1, i )
        end if
        
      ! Outlier test
        if( conv_flag .eqv. .false. .and. mod(i,10) == 0 )then
           call outliertest( i )
        end if        

     end do 
     close(40)
   ! Compute the acceptance rate

     AR = dfloat(accept) / dfloat(Zp_count)         
     !计算均值
     do i=1,dim_pre
     ave(i)=s_pre(i)/(max_gen*nchain-nchain)
     m_ave(i)=sm_pre(i)/(max_gen*nchain-nchain)
     !write(*,*) 'ave=',ave(i)
     end do
     !与95%置信区间
     open(30,file='result.txt',status="replace")
      write(30,'(4a10)') 'ave','mod_ave','low','upper'
     do i=1,dim_pre
         call bubble_sort(sortpre(:,i),max_gen*nchain)
         ai=int(max_gen*nchain*(0.5-api/2))!api为置信值如95%
         bi=int(max_gen*nchain*(0.5+api/2))
         !write(*,*) ai,bi
         low(i)=sortpre(ai,i)
         upper(i)=sortpre(bi,i)
         !write(*,*) yo(1,i),'=',low(i),upper(i)
         write(30,*) ave(i),m_ave(i),low(i),upper(i)
         !write(*,*) sortpre(:,i)
         !可进行输出
     end do
     close(30)
   end subroutine
   
!--------------------------------------------------------------------------------------
!***  Generate candidate parameter samples


   subroutine gen_candidate( gen_num, chain_num )

      implicit none

      integer :: gen_num, chain_num
      integer :: R(2, npair)
      real  :: noise_e(par_dim), b, eps(par_dim)
      integer :: i
     
     
      b=0.0  ! used to calculate e follow U(-b,b)


    ! Pick n pairs of other chains for crossover
      do i = 1, npair
         R(:,i) = rand_R( nchain )
         do while ( R(1,i) == R(2,i) .or. R(1,i) == chain_num .or. R(2,i) == chain_num )
            R(:,i) = rand_R( nchain )
         end do
      end do
	
    ! Determine a jump rate
      call choose_jumprate( gen_num )

    ! Calculate e in eq.4 of Vrugt et al. 2009
      noise_e = b * (2 * rand_uni01_dim(par_dim) - 1) 

    ! Get epsilon value from multinormal distribution                      
      do i =1, par_dim
         eps(i) = gennor(0.0, 1.0e-10)
      end do

    ! Generate the candidate sample based on eq.4 of Vrugt et al. 2009
      Zp = Z(:,chain_num,gen_num-1) + (1.0 + noise_e) * jumprate * comp_diff(gen_num, chain_num, R) + eps

    ! Check whether candidate sample is out of parameter bound
      call outofbound( )

   end subroutine
   
   
   
!--------------------------------------------------------------------------------------
!*** Out of parameter bound test

!  Test whether generated candidate sample is out of the bound set by users
!  if it is out of bound then fold the sample into the bound


   subroutine outofbound( )
   
     implicit none
   
     integer :: i
     
   
     do i = 1, par_dim
   
        if( Zp(i) < range(1,i) )then
            Zp(i) = range(2,i) - range(1,i) + Zp(i)
        
        else if( Zp(i) > range(2,i) )then
            Zp(i) = range(1,i) - range(2,i) + Zp(i)

        end if

! Just in case the sample is still outside bound

        if( Zp(i) < range(1,i) .or. Zp(i) > range(2,i) )then
            Zp(i) = range(1,i) + rand_uni01() * ( range(2,i) - range(1,i) )
        
        end if
     
     end do

   
   end subroutine
   

!********************* Compute crossover probability CR *******************************
!--------------------------------------------------------------------------------------
!*** initialize the CR values

   
   subroutine init_CR( )

      implicit none

      integer :: i
      

      do i = 1, nCR

         CR(i) = dfloat(i) / dfloat(nCR)  
         pCR(i) = 1.0d0 / dfloat(nCR)
         L_CR(i) = 1
         Dis_CR(i) = 1.0d0

      end do

      pCR(nCR) = 1.0d0 - sum(pCR(1:nCR-1))

   end subroutine
   
   
!--------------------------------------------------------------------------------------
!***  choose a CR value


   subroutine choose_CR()

      implicit none
      
      integer :: tmp_ind(nCR) 
      integer :: i
      

      if( nCR == 1 ) then
         CR_ind = 1
         
      else
 
         call genmul(1, pCR(1:nCR), nCR, tmp_ind)

         do i = 1, nCR
            if( tmp_ind(i) == 1 ) then
               CR_ind = i
               exit
            end if
         end do

      end if

   end subroutine


!--------------------------------------------------------------------------------------
!***  update CR distance


  subroutine update_CR_dis( gen_num, chain_num ) 

      implicit none
      integer :: gen_num, chain_num
      real  :: std( par_dim )
      integer :: i
      

    ! Compute standard deviation for all parameters
      std = comp_std( gen_num )

      L_CR(CR_ind) = L_CR(CR_ind) + 1

      do i = 1, par_dim 

         Dis_CR(CR_ind) = Dis_CR(CR_ind) + ( Z(i, chain_num, gen_num ) - Z(i, chain_num, gen_num-1) )**2 / std(i)**2

      end do 

   end subroutine


!--------------------------------------------------------------------------------------
!***  update CR probabilities 

 
    subroutine update_pCR( gen_num )

      implicit none
      
      integer :: i, gen_num


      do i = 1, nCR-1

         pCR(i) = (Dis_CR(i)/L_CR(i)) / sum(Dis_CR(1:nCR) / L_CR(1:nCR)) 

      end do

      pCR(nCR) = 1.0d0 - sum( pCR(1:nCR-1) )
        
   end subroutine
   

!--------------------------------------------------------------------------------------
!***  compute standard deviation 

 
  function comp_std( gen_num )
 
      implicit none

      integer :: gen_num
      real  :: comp_std( par_dim )
      real  :: mean( par_dim )
      integer :: i
      

      do i = 1, par_dim

         mean(i) = sum( sum( Z(i,:,1:gen_num),1 ), 1) / nchain / gen_num
         comp_std(i) = sqrt(sum( sum( (Z(i,:,1:gen_num) - mean(i))**2, 1), 1 ) / (nchain*gen_num-1) )
       
      end do

   end function


!************************* Compute jump rate Gamma ************************************
!--------------------------------------------------------------------------------------
!***  Initiailize the jump rate table


   subroutine init_jumprate( )

      implicit none

      integer :: i
      
      allocate( jump_dim( par_dim ) )
      allocate( jumprate_table( par_dim ) )
      
 
      do i = 1, par_dim
         jumprate_table(i) = 2.38d0 / sqrt( 2.0d0 * npair * i )
      end do

   end subroutine


!--------------------------------------------------------------------------------------
!***  Choose a jump rate from the jump rate table
 
 
    subroutine choose_jumprate( gen_num )

      implicit none
   
      integer :: i
      integer :: gen_num

    ! Determine the dimensions that will be updated
      jump_num = 0
      jump_dim = 0
      

      do i = 1, par_dim
         if( rand_CR( ) > 1-CR(CR_ind) )then
            jump_num = jump_num + 1
            jump_dim(jump_num) = i
         end if
      end do

    ! Calculate general jump rate        
      if( jump_num == 0 )then
         jumprate = 0.0
      else
         jumprate = jumprate_table(jump_num)
      endif

    ! If parameter dimension is 1, 2, or 3, fix the jump rate to 0.6
      if( par_dim == 1 .or. par_dim == 2 .or. par_dim == 3 )then
         jumprate = 0.6
      end if

    ! Determine if do a long jump 
      if( mod(gen_num-1,jumpstep) == 0 )then
         jumprate = 0.98
         return
      end if

   end subroutine



!--------------------------------------------------------------------------------------
!***  Calculate the differential evoluation


   function comp_diff( gen_num, chain_num, R ) 

      implicit none

      integer :: gen_num, chain_num
      real    :: comp_diff( par_dim )
      integer :: R(2, npair)
      integer :: i, j
    
    ! Do differential evolution
      Zp = Z(:,chain_num,gen_num-1)

    ! Produce the difference of the pairs used for population evolution
      comp_diff = 0.0
      
      do i = 1, npair
         do j = 1, jump_num
            comp_diff(jump_dim(j)) = comp_diff(jump_dim(j)) + &
                                   (Z(jump_dim(j),R(1,i),gen_num-1)-Z(jump_dim(j),R(2,i),gen_num-1))         
         end do
      end do

   end function


!**************************************************************************************
!
!               FUNCTIONS ABOUT CONVERGENCE CRETERIA GELMAN-RUBIN R
!
!**************************************************************************************

!--------------------------------------------------------------------------------------
!***  Initialize the Gelman Rubin statistic


   subroutine init_GR()

      implicit none

      allocate( Gel_Ru_R( par_dim, max_gen ) )

      GR_count = 0
      conv_flag = .false.

   end subroutine


!--------------------------------------------------------------------------------------
!***  Compute Gelman Rubin statistics R used to check convergence


   subroutine Gel_Ru( gen_num )

      implicit none

      integer :: gen_num
      real  :: mean_chain( nchain )
      real  :: mean_all
      real  :: S( nchain )
      real  :: W_var, B_var, Var 
      integer :: i, j, ind0


      GR_count = GR_count + 1
      ind0 = 0.5d0*gen_num

      do i = 1, par_dim
    
         do j = 1, nchain
            mean_chain(j) = sum( Z(i,j,ind0:gen_num) ) / dfloat(ind0)
         end do

         mean_all = sum( mean_chain ) / nchain
 
         B_var = dfloat(ind0) / dfloat( nchain - 1 ) * sum( (mean_chain - mean_all)**2 )   

         do j = 1, nchain
            S(j) = sum( (Z(i,j,ind0:gen_num) - mean_chain(j))**2 ) / (ind0 - 1.0d0)
         end do
 
         W_var = sum(S) / nchain

         Var = dfloat(ind0-1)/dfloat(ind0) * W_var + B_var / dfloat(ind0)   

         Gel_Ru_R( i, GR_count ) = sqrt( Var / W_var )

      end do

      conv_flag = .true.
      
      do i = 1, par_dim
         if( Gel_Ru_R( i, GR_count ) > GR_threshold )then
            conv_flag = .false.
            exit
         end if
      end do

      if (conv_flag .eqv. .true.) write(*,*)'Converged at iteration: ',gen_num

   end subroutine 


!**************************************************************************************
!
!                       FUNCTIONS FOR EXPORTING RESULTS 
!
!**************************************************************************************

!--------------------------------------------------------------------------------------
!***  Write parameter samples into file

   subroutine Output_1( ind1, ind2 )

      implicit none

      integer :: ind1, ind2
      integer :: iunit, i, j 
      character*10 :: ic
   

!     write parameter samples of all chains
      do i = 1, nchain

       ! open file
         if( ind1 == 1 )then
            write(ic,'(i2.2)') i
            call get_unit( iunit )
            open(iunit, file = 'ParSamples._chain'//trim(ic),STATUS='REPLACE')
            write(iunit,11)trim(ic)
         else
            write(ic,'(i2.2)') i
            call get_unit( iunit )
            open(iunit, file = 'ParSamples._chain'//trim(ic),STATUS='OLD', position = 'append') 
         end if

       ! Write loglikelihood function and parameter samples  
         do j = ind1, ind2
            write(iunit,10)j, fit(i,j), Z(:,i,j) 
         end do

       ! close file
         close(iunit)   

      end do


!     write Gelman-Rubin R
    ! open file
      if( ind1 == 1 )then
         call get_unit( iunit )
         open( iunit, file = 'Gel_Ru.out',STATUS='REPLACE' )
         write(iunit,21)
      else
         open( iunit, file = 'Gel_Ru.out', STATUS='OLD', position = 'append'  ) 
      end if

    ! Write Gelman-Rubin statistic R
      write(iunit,20) printstep*GR_count, Gel_Ru_R(:,GR_count)

    ! close file
      close(iunit)


11    format(1x,'"MONITORED PARAMETER SAMPLE VALUES AND ASSOCIATED LOG LIKELIHOOD FUNCTION VALUES FOR CHAIN # ',A,'"')
21    format(1x,'"MONITORED PARAMETER INTER-CHAINS GELMAN RUBIN R STATISTIC"')
10    format(1x,i7,6x,es14.7,6x,1000(es14.7,2x)) 
20    format(1x,i9,6x,1000(f14.4,2x)) 

   
   end subroutine


!--------------------------------------------------------------------------------------
!***  Write the last parameter samples into restart file used for the next time simulation


   subroutine OutputRestart( )
   
      implicit none
   
      integer :: iunit, i, j 
   
  
      call get_unit( iunit )
      
      open( iunit, file = 'Restart.out' )
      
      write(iunit,11)

      do i = 1, nchain
         write(iunit,10)i, fit(i, max_gen), Z(:,i,max_gen)
      end do
      
      close(iunit)

         
10    format(1x,i7,7x,es14.7,6x,1000(es14.7,2x)) 
11    format(1x,'"PARAMETER VALUES OF THE LAST GENERATION FOR ALL CHAINS FOR USE IN RESTARTING THE CHAINS"')     


    ! Print out the acceptance rate on creen.   
      write(*,*)'The acceptance rate is: ', AR

   end subroutine


!**************************************************************************************
!
!                       FUNCTIONS TO DETECT OUTLIER CHAIN 
!
!**************************************************************************************   
   

!--------------------------------------------------------------------------------------
! ******* Test outlier chain in burn-in period 

  subroutine outliertest( gen_num )

      implicit none

      integer :: gen_num

      real  :: avg(nchain)
      real  :: avg_tmp(nchain)
      integer :: i, j
      integer :: ind1, ind2
      real  :: IQR, Q1, Q3
      integer :: best(1)


      do i = 1,nchain
         avg(i) = sum( fit(i, gen_num/2:gen_num) ) / size( fit(i, gen_num/2:gen_num) )
      end do 

      best = maxloc( avg )
      avg_tmp = avg

      call sort( nchain, avg_tmp )

      ind1 = nchain * 0.25 + 1
      ind2 = nchain * 0.75 + 1
     
      Q1 = avg_tmp(ind1)
      Q3 = avg_tmp(ind2)

      IQR = Q3 - Q1
  
      do i = 1, nchain
         if( avg(i) < Q1 - 2.0 * IQR )then
            Z(:,i,gen_num) = Z(:,best(1),gen_num)
            fit(i, gen_num/2:gen_num) = fit(best(1), gen_num/2:gen_num)
            write(*,201)i,gen_num,best(1)
         end if
      end do
      
201 FORMAT(2X,'Chain ',i2,' is an outlier chain at iteration ',i10, / &
          '  its samples at this iteration are replaced by those from the chain ',i2,' with the largest log likelihood function.'/)

   end subroutine outliertest

!--------------------------------------------------------------------------------------
!****** Sorting an array in ascending order ************
!****** called by subroutine outliertest ************** 

   SUBROUTINE SORT(n,a)
   
   implicit none

   integer :: i, j, k
   integer :: n
   real  :: a(n), tmp
   
   do i = 1, n-1
      do j = i+1, n
         if( a(i) > a(j) )then
            tmp = a(i)
            a(i) = a(j)
            a(j) = tmp
         end if
      end do
   end do
   
   end subroutine sort


!**********************************************************
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
   !--------------------------------------------------------------------------------------
!****** 模拟值计算************
 function cm(z,x)!compute modeling value,z donates parameters matrix
    real::cm
    real::z(par_dim)
    real::x(in_dim)!自变量的取值
    cm=x(1)*z(4)!在地下水中此为地下水模型
 end function
 subroutine comp_mod(r,m,z)!计算模拟值子程序，对应r,chain;m,gen
 implicit none
 integer::r,m,i
 real::z(par_dim)
 do i=1,ob_dim
 m_comp(i,r,m)=cm(z,yi(:,i))!ob_dim维
 end do
 end subroutine
 
 
!--------------------------------------------------------------------------------------
!****** 矩阵计算************
subroutine cov0(r,m,z)
implicit none
integer::r,k,i,j,m
real::s
real::z(3)!此处定义为三个超参，对于各向异性要提高超参数个数
do i=1,ob_dim
        do j=1,ob_dim
            s=0
            do k=1,in_dim!输入变量维数,简单例子为1
                s=(yi(k,i)-yi(k,j))**2+s!yi为输入变量，可通过输入文件读取
            end do
            if (i==j) then
                cov_0(i,j,r, m)=z(1)*exp(-1/(z(2)*z(2))*s)+z(3)
            else
                cov_0(i,j,r, m)=z(1)*exp(-1/(z(2)*z(2))*s)
            end if
        end do
    end do
end subroutine
subroutine cov1(r,m,z)
implicit none
integer::r,k,i,j,m
real::s
real::z(3)
do i=1,dim_pre
        do j=1,ob_dim
            s=0
            do k=1,in_dim
                s=(yo(k,i)-yi(k,j))**2+s!yi为输入变量，y2为预测输入变量，可通过输入文件读取
            end do
            if (i==j) then
                cov_1(i,j,r, m)=z(1)*exp(-1/(z(2)*z(2))*s)+z(3)
            else
                cov_1(i,j,r, m)=z(1)*exp(-1/(z(2)*z(2))*s)
            end if
        end do
    end do
end subroutine
subroutine cov2(r,m,z)
implicit none
integer::r,k,i,j,m
real::s
real::z(3)!此处定义为三维超参，对于各向异性要提高超参数个数
do i=1,dim_pre
        do j=1,dim_pre
            s=0
            do k=1,in_dim!输入变量维数,简单例子为1
                s=(yo(k,i)-yo(k,j))**2+s!yi为输入变量，可通过输入文件读取
            end do
            if (i==j) then
                cov_2(i,j,r,m)=z(1)*exp(-1/(z(2)*z(2))*s)+z(3)
            else
                cov_2(i,j,r,m)=z(1)*exp(-1/(z(2)*z(2))*s)
            end if
        end do
    end do
end subroutine
!*********预测计算************
subroutine comp_pre(r,m)!计算预测值子程序
 implicit none
 integer::r,m,i
 do i=1,dim_pre
    m_pre(i,r,m)=cm(Z(:,r,m),yo(:,i))!dim_pre维
 end do
 end subroutine
 
 subroutine predict(r,m)!计算最终预测值
 implicit none
 integer::r,m
 real::ita(dim_pre)!dim_pre维测量误差
 real::av_b(dim_pre)!dim_pre维
 real(kind=8)::var_b(dim_pre,dim_pre)!dim_pre*dim_pre维，协方差矩阵
 real(kind=8)::b(dim_pre)!dim_pre维
 integer::seed,i
 av_b=avbfun(r,m)
 var_b=real(varbfun(r,m),kind=8)
 !do i=1,dim_pre
 !write(*,*) var_b(1:dim_pre,i)
 !end do
 call multinormal_sample(dim_pre,var_b,real(av_b,kind=8),seed,b)!得到结构误差b，dim_pre维
 do i=1,dim_pre
 ita(i)=normal(0.0,Z(3,r,m))!得到测量误差ita,dim_pre维
 end do
 write(*,*) 'b='
 write(*,*) b(:)
 pred(:,r,m)=real(b(:),kind=4)+m_pre(:,r,m)+ita(:)!不同点的预测值
 end subroutine
 
 !计算结构误差平均值
 function avbfun(r,m)
 real::avbfun(dim_pre)!dim_pre维
 real::s(ob_dim),c1(ob_dim,1),c2(dim_pre,ob_dim),c3(dim_pre,1)!ob_dim
 integer::r,m
 integer::i,j,k
 do i=1,ob_dim
 c1(i,1)=ob_data(i)-m_comp(i,r,m)
 end do
 c2=matmul(cov_1(:,:,r,m),invercov(:,:,r,m))
 c3=matmul(c2,c1)
 do i=1,dim_pre
     avbfun(i)=c3(i,1)
 end do
 !do i=1,dim_pre
 !    s=0
 !    avbfun(i)=0
 !    do j=1,ob_dim
 !        do k=1,ob_dim
 !         s(j)=cov_1(i,k,r,m)*cov_0(k,j,r,m)+s(j)
 !        end do
 !        avbfun(i)=s(j)*(ob_data(j)-m_comp(j,r,m))+avbfun(i)
 !    end do
 !end do
 end function
 
 !计算结构误差协方差
 function varbfun(r,m)
 real::varbfun(dim_pre,dim_pre),c1(dim_pre,ob_dim),c2(dim_pre,dim_pre)!dim_pre*dim_pre维
 real::s(ob_dim)
 integer::i,j,k,n,r,m
 !varbfun
 c1=matmul(cov_1(:,:,r,m),invercov(:,:,r,m))!c1(dim_pre,ob_dim)
 c2=matmul(c1,transpose(cov_1(:,:,r,m)))!c2(dim_pre,dim_pre)
 varbfun(:,:)=cov_2(:,:,r,m)-c2(:,:)
 !do j=1,dim_pre
 !    do i=1,dim_pre
 !         s=0
 !         varbfun(i,j)=0
 !         do k=1,ob_dim
 !              do n=1,ob_dim
 !              s(k)=cov_1(j,n,r,m)*invercov(n,k,r,m)+s(k)
 !              end do
 !              do n=1,ob_dim
 !              varbfun(i,j)=s(k)*cov_1(i,n,r,m)+varbfun(i,j)
 !              end do
 !         end do
 !         varbfun(i,j)=cov_2(i,j,r,m)-varbfun(i,j)
 !    end do
 !end do
 end function
          
  function ran()   !returns random number between 0 - 1  均匀分布随机数
    implicit none 
    integer , save :: flag = 0
    real :: ran 
    if(flag==0) then 
      call random_seed()
      flag = 1 
    endif 
    call random_number(ran)     ! built in fortran 90 random number function  
   end function ran
   
    function normal(mean,sigma) !产生正太随机数
    implicit none 
    integer :: flag 
    double precision, parameter :: pi = 3.141592653589793239  
    real :: u1, u2, y1, y2, normal, mean, sigma 
    save flag 
    data flag /0/ 
    u1 = ran(); u2 = ran() 
    if (flag.eq.0) then 
      y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2) 
      normal = mean + sigma*y1 
      flag = 1 
    else 
      y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2) 
      normal = mean + sigma*y2 
      flag = 0 
    endif  
    end function normal 
!function r8_uniform_01 ( seed )!产生8字节随机数，用以multinor_module多元正态随机数的产生
!!  Parameters:
!!
!!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!!    NOT be 0. On output, SEED has been updated.
!!
!!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!!    strictly between 0 and 1.
!  integer ( kind = 4 ) k
!  real ( kind = 8 ) r8_uniform_01
!  integer ( kind = 4 ) seed
!
!  if ( seed == 0 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
!    write ( *, '(a)' ) '  Input value of SEED = 0.'
!    stop
!  end if
!
!  k = seed / 127773
!  
!  seed = 16807 * ( seed - k * 127773 ) - k * 2836
!  
!  if ( seed < 0 ) then
!    seed = seed + 2147483647
!  end if
!!
!!  Although SEED can be represented exactly as a 32 bit integer,
!!  it generally cannot be represented exactly as a 32 bit real number!
!!
!  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10
!
!  return
!end function
!****** 计算confidance band************
!对预测样本进行排序
subroutine bubble_sort(a,n)
  implicit none
  integer :: n
  real::a(n),temp
  integer i,j
  do i=n-1,1,-1   ! 开始做n-1次的扫瞄
    do j=1,i      ! 一对一对的来比较，i之后的数字不用比较
    ! 如果a(j) > a(j+1) 就把这两个数值交换
      if ( a(j) > a(j+1) ) then
        temp=a(j)
        a(j)=a(j+1)
        a(j+1)=temp
      end if
    end do
  end do
  return
end subroutine
end module