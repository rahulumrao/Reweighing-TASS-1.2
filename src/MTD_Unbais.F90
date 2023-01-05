MODULE MTD_Unbais
!=========================================================================================================!
! This Subroutine calculates the vbias and ct of MetaD bias which then used to calculate bias probability !
!=========================================================================================================!
USE Input_file        
USE GetSteps
USE GetFileName
USE Error_msg
CONTAINS
!---------------------------------------------------------------------------------------------------------------------------!
!Compute WT-MTD unbiased potential
!
SUBROUTINE mtd_unbiased(au_to_kcal,bias_fact,kt,kb,md_steps,mtd_steps,w_cv,w_hill &
           & ,ncv,t_min,t_max,max_step,gridmin,gridmax,griddif,vbias,ct,nbin,m,mtd,code_name,nr)
IMPLICIT NONE
INTEGER :: j,i_mtd,mtd_steps,i_md,md_steps,mtd_max,m,ir,nr,vt_max
INTEGER :: i_s1,w_cv,w_hill,ncv,t_max,t_min,nbin(ncv)
REAL*8  :: bias_fact,kt,kb,ktb,au_to_kcal
REAL*8  :: dummy1,diff_s2,ds2,ss,hh,num,den,alpha,dum
REAL*8,ALLOCATABLE             :: width(:,:),ht(:,:),hill(:,:),cv(:,:,:)
REAL*8,DIMENSION(ncv)          :: gridmin,gridmax,griddif
REAL*8,DIMENSION(nbin(m))      :: grid
REAL*8,ALLOCATABLE             :: ct(:,:),fes_2D(:)
REAL*8,ALLOCATABLE             :: vbias(:,:),tmp(:,:)
REAL*8                         :: kcons(nr),pcons(nr)
REAL*8,PARAMETER               :: kj_to_kcal = 0.239006
CHARACTER(LEN=5)               :: mtd
CHARACTER(LEN=10)              :: code_name
CHARACTER(LEN=50)              :: filename_loc
CHARACTER(LEN=50)              :: filename(nr),filename_mtd(2,nr)
LOGICAL                        :: max_step
!! #### m = METADYNAMICS CV INDEX
!=============================================================================================================================!
CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)

DO ir = 1,nr
OPEN(22,FILE=filename(ir))
OPEN(23,FILE=filename_mtd(1,ir))
CALL get_steps(22,md_steps)
CALL get_steps(23,mtd_steps)
!--------------------------------------------------------------------------------------------------!
IF (MOD(md_steps,mtd_steps) .ne. 0.0) THEN
PRINT*, "**********************************************************************"
PRINT*, "         ERROR !! MISSING DATA IN colvar OR parvar/hill FILES.        "
PRINT*, "**********************************************************************"
STOP ; ENDIF
CLOSE(22) ; CLOSE(23)
ENDDO
!--------------------------------------------------------------------------------------------------!
   kt = kb*kt
  ktb = kt*bias_fact
alpha = bias_fact/(bias_fact-1.D0)
!--------------------------------------------------------------------------------------------------!
!ALLOCATE(width(nr,mtd_steps),ht(nr,mtd_steps),hill(nr,mtd_steps))
!ALLOCATE(cv(nr,ncv,md_steps))
ALLOCATE(width(nr,999999),ht(nr,999999),hill(nr,999999))
ALLOCATE(cv(nr,ncv,999999))
!--------------------------------------------------------------------------------------------------!
replica_loop: DO ir = 1,nr
IF (code_name .eq. 'CPMD') THEN
OPEN(22,FILE=filename_mtd(1,ir))
OPEN(23,FILE=filename_mtd(2,ir))
CALL get_steps(22,mtd_steps)
!--------------------------------------------------------------------------------------------------!

   DO i_mtd=1,mtd_steps
     READ(22,*) dummy1,dummy1,width(ir,i_mtd),ht(ir,i_mtd)
     READ(23,*) dummy1,hill(ir,i_mtd)
     ht(ir,i_mtd) = ht(ir,i_mtd)*au_to_kcal       ! au_to_kcal 
   END DO
!--------------------------------------------------------------------------------------------------!
ELSEIF (code_name .eq. 'PLUMED') THEN
OPEN(22,FILE=filename_mtd(1,ir))
OPEN(23,FILE=filename(ir))
CALL get_steps(22,mtd_steps)
CALL get_steps(23,md_steps)

   DO i_mtd=1,mtd_steps
     READ(22,*) dummy1,hill(ir,i_mtd),width(ir,i_mtd),ht(ir,i_mtd)
!     WRITE(6,*) dummy1,hill(ir,i_mtd),width(ir,i_mtd),ht(ir,i_mtd)
!     IF(hill(ir,i_mtd) .lt. 0.0) hill(ir,i_mtd) = 0.0d0 
     ht(ir,i_mtd) = ht(ir,i_mtd)*kj_to_kcal        ! kj_to_kcal
   END DO

   DO i_md = 1,md_steps
!    READ(23,*)dummy1,dummy1,cv(ir,u,i_md),dummy1,cv(ir,m,i_md)
     READ(23,*)dummy1,(cv(ir,j,i_md) ,j=1,ncv)
!     IF (cv(ir,m,i_md) .lt. 0.0) cv(ir,m,i_md) = 0.0d0
!     IF (cv(ir,m,i_md) .gt. 2.0) cv(ir,m,i_md) = 2.0d0
!-------------------------------------------------------------------------
DO j = 1,ncv    ! error if cv range is not correct
IF ( cv(ir,j,i_md) .lt. gridmin(j)) CALL cv_error (j,'min') 
IF ( cv(ir,j,i_md) .gt. gridmax(j)) CALL cv_error (j,'max')
ENDDO
   ENDDO
ENDIF
!--------------------------------------------------------------------------------------------------!
ENDDO replica_loop
CLOSE(22) ; CLOSE(23)
!-------------------------------------------------------------------------

DO i_s1 = 1,nbin(m)                                            ! Metadynamics CV bins
    grid(i_s1) = gridmin(m) + DFLOAT(i_s1 - 1)*griddif(m)
ENDDO

!--------------------------------------------------------------------------------------------------!
ALLOCATE(fes_2D(nbin(m)))
!----------------------------------computing ct factor-------------------------------------------!

DO ir = 1,nr
OPEN(22,FILE=filename_mtd(1,ir))
CALL get_steps(22,mtd_steps)            ! Get the MD steps from each window
CLOSE(22)
!PRINT*,'mtd_steps:',mtd_steps

fes_2D = 0.0D0
CALL get_filename('ct_test.dat_',filename_loc,ir)
OPEN(21,FILE=filename_loc,STATUS='unknown') ! open idividual files for every umberlla window
 DO i_mtd=1,mtd_steps
      ds2 = width(ir,i_mtd)*width(ir,i_mtd)
       ss = hill(ir,i_mtd)
       hh = ht(ir,i_mtd)/alpha
num = 0.D0; den = 0.D0 
   DO i_s1=1,nbin(m)  ! Metadynamics CV bins
          diff_s2 = grid(i_s1) - ss
          diff_s2 = diff_s2*diff_s2*0.5D0
     fes_2D(i_s1) = fes_2D(i_s1)+hh*DEXP(-diff_s2/ds2)
             num  = num + DEXP(alpha*fes_2D(i_s1)/kt)
             den  = den + DEXP((alpha -1.0)*fes_2D(i_s1)/kt)
  END DO
       ct(ir,i_mtd) = kt*DLOG(num/den)
       WRITE(21,'(I10,F16.8)')i_mtd,ct(ir,i_mtd)
!PRINT*,ir,i_mtd,ct(ir,i_mtd)
 END DO

CLOSE(21)
ENDDO
DEALLOCATE(fes_2D)
!  call allocate(tmp(size_old+n_enlarge))
!  tmp(1:size_old) = array(1:size_old)
!  call move_alloc(tmp, array)
!--------------------------------------------------------------------------------------------------!
WRITE(*,'(A,/)')'ct values written in ct_test.dat_$'
!!@calculate v(s,t)
!------------------------------------------------------------------------------------------------------!
101 FORMAT (A8,1X,A14,4X,A20,1X,A20,20X,A10,1X,A8,1X,A8)
102 FORMAT (I4,5X,F8.2,11X,F10.2,10X,A30,10X,I8,2X,I8,2X,I8)
write(6,101)'#window ','umbrella_mean','umbrella_k(kcal/mol)','CV_VAL_file',"MD Steps","t_min","t_max"
!------------------------------------------------------------------------------------------------------!

vt_max = t_max !md_steps
DO ir = 1,nr
OPEN(23,FILE=filename(ir))
CALL get_steps(23,md_steps)             ! Get the MD steps from each window
CLOSE(23)

CALL max_t(max_step,t_min,t_max,vt_max,md_steps)   ! Computing t_max [If t_max is not defined then t_max=md_steps]
WRITE(6,102)ir, pcons(ir), kcons(ir), filename(ir),md_steps,t_min,t_max

CALL get_filename('vbias_test.dat_',filename_loc,ir)
OPEN(21,FILE=filename_loc,STATUS='unknown')

 DO i_md = t_min,t_max
    mtd_max = ((i_md - 1)*w_cv/w_hill) + 1
         ss = cv(ir,m,i_md)
    IF(MOD((i_md - 1)*w_cv,w_hill) .eq. 0) mtd_max = mtd_max - 1

dum = 0.d0
DO i_mtd = 1,mtd_max
     ds2 = width(ir,i_mtd)*width(ir,i_mtd)
      hh = ht(ir,i_mtd)/alpha
 diff_s2 = cv(ir,m,i_md) - hill(ir,i_mtd)
 diff_s2 = diff_s2*diff_s2*0.5D0
     dum = dum + hh*DEXP(-diff_s2/ds2)
END DO
   vbias(ir,i_md) = dum
    WRITE(21,'(I10,F16.8)')i_md,vbias(ir,i_md)
!    WRITE(6,*)ir,vbias(ir,i_md)
 END DO
CLOSE(21)
ENDDO
WRITE(*,'(A,/)')'vbias values written in vbias_test.dat_$'
!--------------------------------------------------------------------------------------------------!
END SUBROUTINE mtd_unbiased
END MODULE MTD_Unbais
!---------------------------------------------------------------------------------------------------------------------------!
