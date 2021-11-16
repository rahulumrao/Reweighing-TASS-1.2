!-----------------------------------------------------------------------------------
! This module calculates the statistical error in Molecular Dynamics Simulation run
! https://doi.org/10.1063/1.449733
! Written by Vaishali (vthakkur@iitk.ac.in)
! Modularized by Rahul Verma (vrahul@iitk.ac.in)
!-----------------------------------------------------------------------------------
MODULE Error_estimation
USE GetSteps
USE Input_file
CONTAINS
!===================================================================================
SUBROUTINE statistical_error(cv,nr,u,ncv,mtd,code_name,nblocks)
IMPLICIT NONE
INTEGER                 :: i,j,ir,nr,i_md,u,ncv,dummy,xi,yi
INTEGER                 :: nblocks,bs,skip_steps,nsteps,md_steps
REAL*8,ALLOCATABLE      :: mean(:),var(:),delta_G(:),xmean(:),cv1(:),sum2(:)
REAL*8                  :: cv(nr,ncv,*),pcons(nr),kcons(nr)
REAL*8                  :: b,dum,xcor,dr,dum2,dummy1,sum_mean,k
CHARACTER*5             :: mtd
CHARACTER(LEN=10)       :: code_name
CHARACTER(LEN=50)       :: filename(nr),filename_mtd(2,nr)
REAL*8, PARAMETER       :: kb=1.98E-3 !kcal-1 mol-1 K-1 
!===================================================================================
!bs=30000 !block_size
!skip_steps=15000
!----------------------------------------------------------
CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)
!----------------------------------------------------------
PRINT*,'Calculating Statistical Error ! Please Wait ...!'
!nsteps=5000
OPEN(21,file="variance.dat")
ALLOCATE(var(nr),xmean(999999),cv1(999999),sum2(nr))
!----------------------------------------------------------
102 FORMAT (I2,2X,F16.2,F16.2,2X,I8,4X,A)

replica_loop: DO ir=1,nr

OPEN(11,file=filename(ir),status="old")
CALL get_steps(11,md_steps)
skip_steps = md_steps*(14/100) ! 14% initial steps are skipped
!nblocks=INT((md_steps-skip_steps)/bs)
bs=INT((md_steps-skip_steps)/nblocks)

b=0;dum=0.0d0;xmean=0.0d0
      do j=1,nblocks
sum_mean = 0.0
xi = j*bs
yi = xi+1-bs
         do i_md=yi,xi                          !loop over steps in block-j
!          READ(11,*)dummy,dummy,dummy1,dummy1,dummy1,dummy1,dummy1,cv1(i_md)
!          WRITE(6,*)ir,pcons(ir),i_md,cv(ir,u,i_md),cv1(i_md)
         sum_mean = sum_mean + cv(ir,u,i_md) 
!         sum_mean = sum_mean + cv1(i_md) 
         end do
         xmean(j)=sum_mean/dfloat(bs)
         b=b+1
         dum=dum+xmean(j)
      end do
      dum=dum/dfloat(nblocks)
!----------------------------------------------------------
dummy=0.0d0
      do j=1,nblocks
      var(ir)=dummy+(dum-xmean(j))**2
      end do
      var(ir)=var(ir)/dfloat(nblocks-1)
!      write(21,*)ir,xcor(ir),var(ir)
      write(21,*)ir,pcons(ir),var(ir)
!      write(*,*)ir,pcons(ir),var(ir)
    CLOSE(11)  
END DO replica_loop

!----------------------------------------------------------
CLOSE(21)

OPEN(21,file="delta_G.dat")
103 FORMAT (I2,F8.2,2X,F8.2,3X,F8.6)
WRITE(21,'(A)')'#Statistical Error in MD (https://doi.org/10.1063/1.449733)'
WRITE(21,'(A,A,A,A)')'#nr  Umb_mean delta_G sqrt(delta_G)'
ALLOCATE(delta_G(nr))

!----------------------------------------------------------
 sum2(1)=0.0d0     !nr is the referenece point
 DO ir=nr,2,-1
   sum2(ir)=(var(ir)*kcons(ir)**2+var(ir-1)*kcons(ir-1)**2)*0.25d0
   dr=(pcons(ir)-pcons(ir-1))**2
   sum2(ir)=sum2(ir)*dr
 END DO 
 dum2=0.0d0
 DO ir=1,nr
    dum2=dum2+sum2(ir)
    delta_G(ir)=dum2
    WRITE(21,103)ir,pcons(ir),delta_G(ir),dsqrt(delta_G(ir))
 END DO
!-----------------------------------------------------------
END SUBROUTINE Statistical_error
END MODULE error_estimation
