module run_mod
    implicit none
    contains
    subroutine Run_Model(para, m, head,d,p_head,p,yi,yo)
    integer::i,m,d,p
    real*8::para(m),head(d),p_head(p)
    real*8::yi(1,d),yo(1,p)
    do i=1,d
        head(i)=yi(1,i)*para(4)
    end do
    do i=1,p
        p_head(i)=yo(1,i)*para(4)
    end do
    end subroutine
    !subroutine Run_Model(para, m, head,d,p_head,p)
    !    !para 参数
    !    !d参数维数
    !    !head 观测点输出水头
    !    !d观测点个数
    !    !p_head预测点水头
    !    !p预测点个数
    !    Use dflib
    !    !Use Model_data2
    !    Use Model_data3
    !    !--
    !    implicit none
    !    !定义变量
    !    integer i, j, m, d,p, chd(360, 3), n, re
    !    real*8 para(m), K1(r1, r2), K3_tran(r1, r2), vk1(r1, r2), vk2(r1, r2),p_head(p), head(d), log_k1(r1, r2)
    !    character*25 headm, Fname
    !    real*8 sigma, l
    !    logical alive
    !    Fname='M4.nam'
    !    !----
    !    do i=1,120
    !      chd(i, 1)=1
    !      chd(i, 2)=i
    !      chd(i, 3)=1
    !    end do
    !    do i=121, 240
    !      chd(i, 1)=2
    !      chd(i, 2)=i-120
    !      chd(i, 3)=1
    !    end do
    !    do i=241, 360
    !      chd(i, 1)=3
    !      chd(i, 2)=i-240
    !      chd(i, 3)=1
    !    end do
    !    !----read the alternative parameters
    !    !--par1-recharge ratio
    !    open(unit=19, file='M4.rch', recl=50, pad='yes')
    !    write(unit=19, fmt='(i1, i2)') 1, 0
    !    write(unit=19, fmt='(i1)') 1
    !    write(unit=19, fmt='(a8, f14.10)') 'constant', para(1)*9.33333333E-4
    !    close (19)
    !    !--
    !    !--par2-constant head
    !    open(unit=16, file='M4.chd', recl=50, pad='yes')
    !    write(unit=16, fmt='(i3)') 360
    !    write(unit=16, fmt='(i3, i2)') 360, 0
    !    do j=1, 360
    !      write(unit=16, fmt='(i1, i4, i2, f10.6, f10.6)') chd(j, 1), chd(j, 2), chd(j, 3), para(2), para(2)
    !    end do
    !    close(16)
    !    !--par3-riverbed conductance
    !    open(unit=15, file='M4.riv', recl=50, pad='yes')
    !    write(unit=15, fmt='(i3, i2)') 120, 0
    !    write(unit=15, fmt='(i3, i2)') 120, 0
    !    do j=1, 120
    !      write(unit=15, fmt='(i1, i4, i4, f5.1, f12.6, f5.1)') 1, j, 200, 35.0, para(3), 30.0
    !    end do
    !    close(15)
    !    !--
    !    !--for K field
    !    sigma=para(4)
    !    l=para(5)
    !    !--
    !    call K_Kriging(sigma, l, log_k1, r1, r2)
    !    !--
    !    do i=1, r1
    !      do j=1, r2
    !        K1(i, j) = exp(log_k1(i, j))
    !      end do
    !    end do
    !    !--
    !    !-K3_tran means transmisivity of layer 3
    !    !--
    !    open (unit=131, file='k3_tran.txt', recl=5000, blank='null', position='rewind', pad='yes')
    !    do i=1, r1
    !      read(unit=131, fmt='(<r2>f14.6)') k3_tran(i, 1: r2)
    !    end do
    !    close (131)
    !    !--
    !    do i=1, r1
    !      do j=1, r2
    !        vk1(i, j)=1.0/(175.0/k1(i, j)+250.0)
    !        vk2(i, j)=1.0/(250.0+ 2000.0/k3_tran(i, j))
    !      end do
    !    end do
    !    !--
    !    open(unit=14, file='M4.bcf6', recl=5000, pad='yes')
    !    write(unit=14, fmt='(i1, f8.3, i2, 2f4.1, i2)') 0, 999.999, 0, 1.0, 1.0, 1
    !    write(unit=14, fmt='(a2, a3, a3)') '01', '02', '02'
    !    write(unit=14, fmt='(a8, i2)') 'constant', 1
    !    write(unit=14, fmt='(a8, f4.1, a7, i2)') 'internal', 1.0, '(free)', 0
    !    !--
    !    do i=1, r1
    !      write(unit=14, fmt='(<r2>f16.6)') K1(i, 1: r2)
    !    end do
    !    !--
    !    write(unit=14, fmt='(a8, f4.1, a7, i2)') 'internal', 1.0, '(free)', 0
    !    do i=1, r1
    !      write(unit=14, fmt='(<r2>f12.8)') vk1(i, 1: r2)
    !    end do
    !    !--
    !    write(unit=14, fmt='(a8, f4.1)') 'constant', 0.5
    !    write(unit=14, fmt='(a8, f4.1, a7, i2)') 'internal', 1.0, '(free)', 0
    !    do i=1, r1
    !      write(unit=14, fmt='(<r2>f12.8)') vk2(i, 1: r2)
    !    end do
    !    !--
    !    write(unit=14, fmt='(a8, f4.1, a7, i2)') 'internal', 1.0, '(free)', 0
    !    do i=1, r1
    !      write(unit=14, fmt='(<r2>f14.6)') K3_tran(i, 1: r2)
    !    end do
    !    !--
    !    close(14)
    !    !----end for parameters reading
    !    !--
    !      n=RunQQ ('mf2005.exe', Fname)
    !    !--
    !    !----read data from output files
    !    !--
    !    !inquire (file = 'head.txt', exist = alive)
    !      open(unit=30, file='head.txt', recl=50, pad='yes')
    !      read(unit=30, fmt='(a10)') headm
    !      do i=1, d
    !        read(unit=30, fmt='(f20.6)') head(i)
    !      end do
    !      do i=1,p
    !        read(unit=30, fmt='(f20.6)') p_head(i)
    !      end do 
    !      close(30)
    !      !--delete the output file
    !      !re=systemqq('del head.txt')
    !end subroutine
    !!----------------------------------------------------------------------------------------------------
    !Subroutine K_Kriging(sigma, l, field, c1, c2)
    !Use dflib
    !Implicit none
    !!--
    !integer i, j, k, c1, c2, n
    !real*8 field(c1, c2), sigma, l, samples(c1*c2)
    !character*20 head
    !logical alive
    !!--
    !open(unit=110, file='sgsim.par', recl=50, pad='yes')
    !write(unit=110, fmt='(a19)') 'START OF PARAMETERS'
    !write(unit=110, fmt='(a16)') 'lnK_measures.txt'
    !write(unit=110, fmt='(i1, i2, i2, i2, i2, i2)') 1, 2, 0, 3, 0, 0
    !write(unit=110, fmt='(D10.2, D10.2)') -1.0E16, 1.0E16
    !write(unit=110, fmt='(i1)') 0
    !write(unit=110, fmt='(a9)') 'sgsim.trn'
    !write(unit=110, fmt='(i1)') 0
    !write(unit=110, fmt='(a12)') 'histsmth.out'
    !write(unit=110, fmt='(i1, i2)') 0, 0
    !write(unit=110, fmt='(f3.1, f5.1)') 0.0, 15.0
    !write(unit=110, fmt='(i1, f4.1)') 1, 0.0
    !write(unit=110, fmt='(i1, f5.1)') 1, 15.0
    !write(unit=110, fmt='(i1)') 0
    !write(unit=110, fmt='(a9)') 'sgsim.dbg'
    !write(unit=110, fmt='(a9)') 'sgsim.txt'
    !write(unit=110, fmt='(i1)') 1
    !write(unit=110, fmt='(i3, f5.1, f5.1)') 200, 12.5, 25.0
    !write(unit=110, fmt='(i3, f5.1, f5.1)') 120, 12.5, 25.0
    !write(unit=110, fmt='(i1, f5.1, f5.1)') 1, 12.5, 25.0
    !write(unit=110, fmt='(i5)') 69069
    !write(unit=110, fmt='(i1, i3)') 0, 48
    !write(unit=110, fmt='(i1)') 3
    !write(unit=110, fmt='(i1)') 1
    !write(unit=110, fmt='(i1, i2)') 1, 3
    !write(unit=110, fmt='(i1)') 0
    !write(unit=110, fmt='(f6.1, 2f7.1)') 2000.0, 2000.0, 2000.0
    !write(unit=110, fmt='(f3.1, 2f4.1)') 0.0, 0.0, 0.0
    !write(unit=110, fmt='(i3, i4, i2)') 120, 200, 1
    !write(unit=110, fmt='(i1, f4.1, f4.1)') 0, 0.6, 1.0
    !write(unit=110, fmt='(a10)') 'nodata.txt'
    !write(unit=110, fmt='(i1)') 0
    !write(unit=110, fmt='(i1, f4.1)') 1, 0.0
    !write(unit=110, fmt='(i1, f10.6, 3f4.1)') 2, sigma, 0.0, 0.0, 0.0
    !write(unit=110, fmt='(f12.4, f12.4, f5.1)') 3.0*l, 3.0*l, 25.0
    !close(110)
    !!--
    !n=RunQQ ('sgsim.exe', 'sgsim.par')
    !!--
    !inquire (file = 'sgsim.txt', exist = alive)
    !if (alive) then
    !!--
    !  open (unit=111, file='sgsim.txt', recl=25, blank='null', position='rewind', pad='yes')
    !  do i=1, 3
    !    read(unit=111, fmt='(a20)') head
    !  end do
    !!--
    !  do i=1, c1*c2
    !    read(unit=111, fmt='(f14.8)') samples(i)
    !  end do
    !  close(111)
    !!--
    !  n=systemqq('del sgsim.txt')
    !!--
    !  do i=1, c1
    !    do j=1, c2
    !      k=(120-i)*c2 + j
    !      field(i, j)=samples(k)
    !    end do
    !  end do
    !!--
    !else
    !  write(*, *) 'The sgsim.txt does not exist!!!'
    !  field=exp(1.0)*sqrt(sigma)
    !end if
    !!--
    !End subroutine
end module
    
    
    
