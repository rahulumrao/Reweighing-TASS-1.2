MODULE US_Prob
USE MTD_Potential        
USE Input_file        
USE GetSteps
USE GetFileName
!USE wham_code
CONTAINS
!================================================== US Probability =========================================================!
!Compute 1 Dimensional Probability along Umbrella Coordinate
!
SUBROUTINE oneD_prob(ii,nr,u,m,mtd,max_step,t_min,t_max,md_steps,den,prob,prob_mtd,ncv,cv,nbin &
           & ,gridmin,griddif,norm,code_name)
IMPLICIT NONE
INTEGER :: ii,ir,nr,u,m,i_md,t_min,vt_max,t_max,index1,md_steps,i_s1,i_s2,ncv,nbin(*)
REAL*8  :: den,dum,s1,s2,prob(nbin(u)),gridmin(*),griddif(*),cv(nr,ncv,*)
REAL*8  :: prob_mtd(nr,nbin(u),nbin(m))
REAL*8  :: prob_1D(nr,99999),pcons(nr),kcons(nr),norm(*)
LOGICAL :: max_step
CHARACTER*5  :: mtd
CHARACTER*10 :: code_name
CHARACTER(LEN=50) :: filename_loc
CHARACTER(LEN=50) :: filename(nr),filename_mtd(2,nr)
!-----------------------------------------------------------------------------------!
!ALLOCATE(filename(nr))
!ALLOCATE(filename_mtd(nr))
!-----------------------------------------------------------------------------------!
CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)
vt_max = t_max
!PRINT*,'t_min = ',t_min
!-----------------------------------------------------------------------------------!
IF (ii .eq. u .and. mtd .eq. 'y') THEN    ! on;y if MetaD bias is enabled
OPEN(22,FILE=filename(ir))
DO ir = 1,nr
CALL get_steps(22,md_steps)
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)  ! get t_max
!-----------------------------------------------------------------------------------!
CALL get_filename('PROB.dat_',filename_loc,ir) ! individual file to write probability
OPEN(2,FILE=filename_loc,STATUS='unknown')

  DO i_s1 = 1,nbin(u)
   s1 = DFLOAT(i_s1-1)*griddif(u) + gridmin(u)
    dum = 0.0d0 ; den = 0.0 ;  prob_1D = 0.d0
     DO i_s2 = 1,nbin(m)
      s2 = DFLOAT(i_s2-1)*griddif(m) + gridmin(m)
      dum = dum + prob_mtd(ir,i_s1,i_s2)          ! integrating MetaD biased probability to get prob along umbrella CV
! WRITE(2,*)s1,s2,prob_mtd(ir,i_s1,i_s2)*norm(ir)
     ENDDO
  prob_1D(ir,i_s1) = dum*griddif(m)
  den = den + dum*griddif(m)
  WRITE(2,*)s1,prob_1D(ir,i_s1)*norm(ir)
  ENDDO
ENDDO
ENDIF
CLOSE(22)
!-----------------------------------------------------------------------------------!
IF (ii .eq. u .and. mtd .eq. 'n') THEN
DO ir = 1,nr
OPEN(22,FILE=filename(ir))
CALL get_steps(22,md_steps)
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)
!PRINT*,'t_max = ',t_max

prob=0.0d0 ; den = 0.d0
 DO i_md=t_min,t_max
 index1 = nint((cv(ir,u,i_md) - gridmin(u))/griddif(u)) +1
!-----------------------------------------------------------------------------------!
IF (index1 .gt. nbin(u)) THEN
PRINT*, "            ***************************************************"
PRINT*, "                ERROR !! UMBRELLA CV RANGE IS NOT CORRECT"
PRINT*, "            ***************************************************"
STOP ; ENDIF
!-----------------------------------------------------------------------------------!
    prob(index1) = prob(index1) + 1.d0
 END DO
den=0.0
DO index1 = 1,nbin(u)
    den = den + prob(index1)
END DO

CALL get_filename('PROB.dat_',filename_loc,ir)
OPEN(2,FILE=filename_loc,STATUS='unknown')
DO i_s1 = 1,nbin(u)
     s1 = DFLOAT(i_s1-1)*griddif(u) + gridmin(u)
     IF (mtd .eq. 'y') THEN
     WRITE(2,*)s1,prob(i_s1)/(den*griddif(u)*griddif(m))
     ELSE
     prob_1D(ir,i_s1) = prob(i_s1)/(den*griddif(u))
     WRITE(2,*)s1,prob_1D(ir,i_s1)
     ENDIF  
ENDDO ; ENDDO ; ENDIF
!-----------------------------------------------------------------------------------!
CLOSE(2) ; CLOSE(22)
WRITE(*,'(A)')"Unbiased 1D distribution along US  written in 'PROB.dat!'"
END SUBROUTINE
END MODULE US_Prob
