!--------------------------------------------------------------------------------
! THIS SUBROUTINE COMPUTE'S MEAN FORCE AT EACH UMBRELLA WINDOW AND THEN 
! INTEGRATE THAT MEAN FORCE TO CALCULATE FREE ENERGY (1D) ALONG UMBRELLA CV
! u --> UMBRELLA CV INDEX ; nr --> NUMBER OF UMBRELLA WINDOWS
! pcons --> POSITION OF UMBRELL ; kcons --> Kappa VALUE IN EACH UMBRELLA WINDOW
! dfds = delF/delS ; fes --> FREE ENERGY ALONG UMBRELLA CV
! THE 2D FREE ENERGY IS COMPUTED USING NEW ALGORITM OF HISTOGRAMMING
! Re-WRITTEN BY : Rahul Verma (vrahul@iitk.ac.in)
!--------------------------------------------------------------------------------
MODULE MeanForce
USE Input_file        
USE GetSteps
USE GetFileName
USE MTD_Unbais
CONTAINS
SUBROUTINE mean_force(max_step,u,cv_num,prob_nD,ncv,cv,nr,kt,nbin,t_min,t_max,pcons,kcons, &
                 & gridmin,gridmax,griddif,fes,mtd,w_cv,w_hill,ct,vbias,code_name)
IMPLICIT NONE
INTEGER                 :: i,j,i_md,i_mtd,t_min,t_max,nr
INTEGER                 :: ncv,w_cv,w_hill,md_steps,prob_nD
INTEGER                 :: nbin(*),u,cv_num(*),vt_max,cv1,cv2
INTEGER                 :: i_s1,i_s2,index1,index2,ir
REAL                    :: s1,s2
REAL*8                  :: diff_s,den,num,kt,dum
REAL*8                  :: cut_min, cut_max, a, b
REAL*8                  :: gridmin(*),gridmax(*),griddif(*)
REAL*8                  :: vbias(nr,*),ct(nr,*),cv(nr,ncv,*)
REAL*8                  :: pcons(nr),kcons(nr)
REAL*8,ALLOCATABLE      :: dummy(:,:,:),prob(:),prob_2D(:,:),fes_2d(:,:)
REAL*8,ALLOCATABLE      :: dfds(:,:),av_dfds(:),norm(:),fes(:),diff(:)
REAL*8, PARAMETER       :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER       :: au_to_kcal = 627.51
REAL*8, PARAMETER       :: kj_to_kcal = 0.239006

LOGICAL                 :: max_step
!LOGICAL                 :: read_ct,read_vbias
CHARACTER*5             :: mtd
!CHARACTER(LEN=50)       :: filename_loc
CHARACTER(LEN=10)       :: code_name
CHARACTER(LEN=50)       :: filename(nr),filename_mtd(2,nr)

IF(nbin(1) .eq. 0 ) STOP "ERROR : NUMBER OF BINS CAN NOT BE ZERO"

kt = kt*kb
PRINT*, "NOTE: kB T0 (kcal/mol)        = ", kt

CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)
!------------------------------------------------------------------------------------------------------!
101 FORMAT (A8,1X,A14,4X,A20,1X,A20,20X,A10,1X,A8,1X,A8)
102 FORMAT (I4,5X,F8.2,11X,F10.2,10X,A30,10X,I8,2X,I8,2X,I8)
write(*,101)'#replica','umbrella_mean','umbrella_k(kcal/mol)','CV_VAL_file',"MD Steps","t_min","t_max"
!------------------------------------------------------------------------------------------------------!
cv2 = cv_num(prob_nD)           ! cv 2 index
ALLOCATE (dfds(nr,9999999))
ALLOCATE (av_dfds(nr))
ALLOCATE (prob_2D(nr,nbin(cv2)))
ALLOCATE (fes_2D(nr,nbin(cv2)))
ALLOCATE (norm(nr))
ALLOCATE (diff(nr))
!------------------------------------------------------------------------------------------------------!

vt_max = t_max
replica_loop: DO i = 1,nr
OPEN(11,FILE=filename(i),status='old')
CALL get_steps(11,md_steps)
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)   ! get t_max

WRITE(*,102)i, pcons(i), kcons(i), filename(i),md_steps,t_min,t_max
!------------------------------------------------------------------------------------------------------!
! 1D Mean Force !
!---------------!
DO i_md=1,md_steps
  diff_s = cv(i,u,i_md) - pcons(i)
  dfds(i,i_md) = -diff_s*kcons(i)                 ! calculating df/ds
ENDDO
!
IF (mtd .eq. 'y') THEN
den = 0.d0 ; num = 0.d0 ; dum = 0.d0
!!REWEIGHTING THE TIME DEPENDEND BIAS
DO i_md = t_min,t_max
  i_mtd = ((i_md-1)*w_cv/w_hill) + 1
    dum = vbias(i,i_md) - ct(i,i_mtd)
    num = num + dfds(i,i_md)*DEXP(dum/kt)
    den = den + DEXP(dum/kt)
ENDDO
ELSEIF (mtd .eq. 'n') THEN
den = 0.d0 ; num = 0.d0 ; dum = 0.d0
   DO i_md = t_min,t_max
     num = num + dfds(i,i_md)
     den = den + DEXP(dum/kt)
   ENDDO
!PRINT*,'num =  ',num
ENDIF
   av_dfds(i) = num/den                            ! average of df/ds along every umbrella
ENDDO replica_loop  ! end loop for each replica
!----------------------------------------------------------------------------
103 FORMAT (A,I1,/,A,I1,A,I1,A,/,'Please Wait ...')
WRITE(6,103)"1D Free Energy along CV:",u, &
            "2D Free Energy along CV:",u,' and CV:',cv2,' will be computed'
!----------------------------------------------------------------------------
OPEN(12,FILE='av_dfds.dat')
DO i = 1,nr
WRITE(12,*)pcons(i),av_dfds(i)
ENDDO
!------------------------------------------------------------------------------------------------------!
ALLOCATE(prob(nbin(1)))

OPEN(13,FILE="free_energy.dat")
!
fes = 0.d0 ; num = 0.d0
!Integration using TrapeZoidal Rule 
DO i = 1,nr-1
   dum = pcons(i+1) - pcons(i)
   num = num + dum*(av_dfds(i+1) + av_dfds(i))
   fes(i+1) = num*0.5d0
   WRITE(13,*)pcons(i+1),fes(i+1)
ENDDO
!
DEALLOCATE(av_dfds,dfds)
CLOSE(11)
!------------------------------------------------------------------------------------------------------!
! 2D Mean Force   !
!-----------------!
prob_2D = 0.0d0
! Calculating difference between two adjecent umbrella for normalization of probability
DO i = 1,nr
IF (i .eq. 1) diff(i) = pcons(2) - pcons(1) !griddif(u) 
IF (i .gt. 1) diff(i) = pcons(i) - pcons(i-1)
ENDDO
! Looping over every Umbrella Window
r_loop: DO i = 1,nr 
norm(i) = 0.d0 ; den = 0.d0
open(11,file=filename(i),status='old')
md_steps = 0
CALL get_steps(11,md_steps)
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)   ! get t_max
! Finding the index along CV1 by checking the value of cv at every step
  DO i_md = t_min,t_max                 ! Over MD steps (if t_max not defined then t_max=md_steps)
      DO j = 1,nr                       ! Over the index of every umbrella mean position
        IF (j .eq. 1) THEN
           a = (pcons(j+1) - pcons(j))/2.0
           cut_min = pcons(1)
           cut_max = pcons(1) + (a - 0.01)
         ELSEIF (j .gt. 1 .and. j .lt. nr)THEN
           a = (pcons(j+1) - pcons(j))/2.0
           b = (pcons(j) - pcons(j-1))/2.0
           cut_min = pcons(j) - b  
           cut_max = pcons(j) + a
         ELSEIF (j .eq. nr) THEN
           b = (pcons(j) - pcons(j-1))/2.0
           cut_min = pcons(nr) - b
           cut_max = pcons(nr)
         ENDIF

        IF (cv(i,u,i_md) .ge. cut_min .and. cv(i,u,i_md) .le. cut_max)THEN
        index1 = j
        ENDIF
     ENDDO
! Calculating index along CV2 using nomral histogramming
        index2 = NINT((cv(i,cv2,i_md) - gridmin(cv2))/griddif(cv2)) + 1
! Histogramming        
     IF((index1 .gt. 0 .and. index2 .gt. 0) .and. (index1 .le. nr .and. index2 .le. nbin(cv2))) THEN
         prob_2D(index1,index2) = prob_2D(index1,index2) + 1.0
     ENDIF
  ENDDO

DO index1 = 1,nr
  DO index2 = 1,nbin(cv2)
    den = den + prob_2D(index1,index2)
  ENDDO
ENDDO
!Normalization
norm(i) = den*diff(i)*griddif(cv2)
norm(i) = 1.0d0/norm(i)
!PRINT*,i,norm(i),den,diff(i)
CLOSE(11)
END DO r_loop ! END replica_loop
!!-----------------------------------------------------------------------------------------------------
OPEN(14,FILE='free_energy_2D.dat',STATUS='unknown')
DO i_s1 = 1,nr
   s1 = pcons(i_s1)
   DO i_s2 = 1,nbin(cv2)
      s2 = FLOAT(i_s2-1)*griddif(cv2) + gridmin(cv2)
      num = prob_2D(i_s1,i_s2)*norm(i_s1) 
      fes_2D(i_s1,i_s2) = -kt*dlog(max(num,1E-32)) + fes(i_s1)
      IF (fes_2D(i_s1,i_s2) .gt. 50.0 ) THEN 
      WRITE(14,'(2E16.8,2X,A,1E16.8)')s1, s2,'+Infinity'
      ELSE
      WRITE(14,'(5E16.8)')s1, s2,fes_2D(i_s1,i_s2) 
      ENDIF
   END DO
   WRITE(14,*)
END DO
!------------------------------------------------------------------------------------------------------!
WRITE(6,'(A)')"Free Energies are Written in file -> 'free_energy.dat' and 'free_energy_2D.dat'"   
CLOSE(12) ; CLOSE(13) ; CLOSE(14)
END SUBROUTINE 
END MODULE MeanForce

