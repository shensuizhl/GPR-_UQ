Module Model_data1
implicit none
!--
integer, parameter :: r1=50
integer, parameter :: r2=100
!--
End Module
module run_mod
    implicit none
    real*8::lens(2)
    data lens /500.0,250.0/
    contains
Subroutine init_model()
use Model_data1
Use dflib
implicit none
integer,parameter::inputfileid=10
integer,parameter::outputfileid=20
real::k1(r2,r1),lnk1(r2,r1)
integer::i,j,k,n
integer::status=0
real::buffer
character*100  Fname,line
!interpolation the value of lnK1
Fname='kb2d.par'!插值点需要改变对f1进行更新。
n=RunQQ ('kb2d.exe', Fname)
!layer1 lnK
!lnK1.txt->translnK1.txt
open(unit=inputfileid,file='kb2d.out',form="formatted",status="old")
open(unit=outputfileid,file='translnK1.txt',status="replace")
do i=1,4
    read(inputfileid,*)
end do
do i=1,50
    do j=1,100
       read(inputfileid,*) buffer
       lnk1(j,i)=buffer
       k1(j,i)=exp(lnk1(j,i))
       write(outputfileid,"(1Xi4,1Xi4,2Xf12.8)") (10*j-5),(10*i-5),lnk1(j,i)
    end do
end do
close(inputfileid)
close(outputfileid)
!rewrite the file of *.lpf
open(30,file='virtual1.lpf',form="formatted",status="replace")
!layer2 k
!k2.txt
open(63,file='k2.txt',status="old")
write(30,"(I2,1X,F6.1,1X,I1)") 14,-888.0,0
write(30,*) '1 1 1 '
write(30,*) '0 0 0 '
write(30,*) '-1.0 -1.0 -1.0 '
write(30,*) '1 1 1 '
write(30,*) '0 0 0 '
!layer1
write(30,*) 'INTERNAL 1 (free) -1'
do j=1,50
    do i=1,10
    write(30,"(10(f5.1,1X))") k1((10*i-9):(10*i),j)
    end do
end do
write(30,*) 'CONSTANT 1.0'
write(30,*) 'CONSTANT 3.0'
!layer2
write(30,*) 'INTERNAL 1 (free) -1'
do j=1,500
    read(63,"(a100)") line
    !write(*,*) line(:)
    write(30,"(a100)") line
end do
write(30,*) 'CONSTANT 1.0'
write(30,*) 'CONSTANT 3.0'
!layer3
write(30,*) 'CONSTANT 10.0'
write(30,*) 'CONSTANT 1.0'
write(30,*) 'CONSTANT 3.0'
close(63)
close(30)
end Subroutine
!run model
Subroutine Run_Model(para,d, m_c,mi,p_c,p,locat_p)
   !para 参数
    !d参数维数
    !m_c 观测点输出浓度
    !mi观测点个数
    !p_c预测点浓度
    !p预测点个数
    !locat_m观测点坐标
    !locat_p预测点坐标
use Model_data1
Use dflib
implicit none
integer,parameter::inputfileid=10
integer,parameter::outputfileid=20
integer::i,j,k,n,bufferint,fli,ki(12),ii(12),ji(12),ti,d,mi,p
integer::status=0
real::buffer,bufferl(4)
real*8::para(8),locat_p(2,p)
!real::line(10)
character*100  Fname,line
real*8::concentration(27),m_c(mi),p_c(p)
!write the parameter file
!porosity ->*.btn file
open(34,file='virtual10.btn',status="old")
open(35,file='virtual1.btn',status="replace")
do i=1,12
read(34,'(a100)') line
write(35,'(T1,a)') trim(line)
end do 
read(34,*) line
write(35,"(i10,f10.7)") 0,para(4)
do i=1,11
read(34,'(a100)') line
!write(*,*) trim(line)
write(35,'(T1,a)') trim(line)
end do 
!T->*.obs file
ti=nint(para(8))
read(34,*) line
write(35,"(f10.1)") (ti+100.0)

do i=1,35
read(34,'(a100)') line
!write(*,*) trim(line)
write(35,'(T1,a)') trim(line)
end do 
close(35)
close(34)
!Q->*.ssm file
!read the source point from *.txt
open(37,file='source_point.txt',status="old")
do i=1,12
    read(37,*)  ki(i),ii(i),ji(i)
end do
close(37)
open(38,file='virtual1.ssm',status="replace")
write(38,'(T1,a)') ' F F F F F F'
write(38,'(T1,i10)') 312
do j=1,3
    write(38,'(T1,i10)') 12
    do i=1,12
        write(38,"(T1,3i10,f10.3,i10,f10.3)") ki(i),ii(i),ji(i),para(5)*1000,15,para(5)*1000
    end do
end do
close(38)
!αL、 αT->*.dsp file
open(36,file='virtual1.dsp',status="replace")
write(36,*) 0,para(6)
write(36,*) 0,10.0
write(36,*) 0,40.0
write(36,*) 0,para(7)/para(6)
write(36,*) 0,para(7)/para(6)
write(36,*) 0,0
close(36)

Fname='virtual1.nam'
n=RunQQ ('mf2005.exe', Fname)
Fname='virtual1mt.nam'
n=RunQQ ('mt3dms5b.exe', Fname)
!get the data of observation
open(unit=50, file="mt3d001.obs", recl=500, pad='yes')
do i=1,2*ti+1
   read(50,*) line
end do
read(50,*) bufferint,buffer,concentration(1:16)
read(50,*) concentration(17:27)
close(50)
!lnk of observation
open(61,file='xypoint.txt',status="old")
!open(62,file='xylnk.txt',status="replace")
open(65,file='xyconcentration.txt',status="replace")
do k=1,27
   read(61,*) i,j
   !write(62,"(1Xi4,1Xi4,2Xf12.8)") (10*j-5),(10*i-5),lnk1(j,i)
   write(65,*) (10*j-5),(10*i-5),concentration(k)
end do
close(61)
!close(62)    
close(65)
!p_c=concentration(1:9)
m_c=concentration(1:27)
!get the poluted field
!need to revise pm_test.nam file
open(66,file='pm_test.nam',status="replace")
   write(66,'(T1,a)') 'mt3d001.ucn'
   write(66,'(T1,a)') '2'
   write(66,'(T1,a)') 'mt3d.cnf'
   write(66,'(T1,a)') 'y'
   write(66,'(T1,a)') '0 0 0'
   write(66,'(T1,a)') '0'
   write(66,*) ti+100-int((ti+100.0)/360)*360,1,int((ti+100.0)/360)+1
   !write(*,*) ti+100-int((ti+100.0)/360)*360,1,int((ti+100.0)/360)
   write(66,'(T1,a)') '1 1 1'
   write(66,'(T1,a)') '100 50 1'
   write(66,'(T1,a)') 'pm_test.dat'
   write(66,'(T1,a)') '2'
   write(66,'(T1,a)') 'n'
close(66) 
!get pm_test.dat file
fli=systemqq('pmnew.exe<pm_test.nam')
open(67,file='pm_test.dat',status="old")
do i=1,5000
   read(67,*)  bufferl(1),bufferl(2),bufferl(3),bufferl(4)
   do j=1,p
       if ((bufferl(1)==locat_p(1,j)).and.(bufferl(2)==locat_p(2,j))) then
           p_c(j)=bufferl(4)
           !write(*,*) p_c(j)
       end if
   end do
end do
close(67) 
!stop
end subroutine
end module