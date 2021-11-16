MODULE US_TEMP
USE Input_file
USE GetSteps
USE GetFileName
USE Error_msg
CONTAINS
!=================================================== US vs TEMP Probability ================================================!
!Compute Unbiased Distribution along Umbrella Coordinate vs One of the Temperature CV
!
SUBROUTINE twoD_temp_prob(jj,u,nr,max_step,t_min,t_max,md_steps,prob_2D,ncv,cv,nbin &
                         & ,gridmin,griddif,mtd,code_name)
IMPLICIT NONE
INTEGER :: jj,u,nr,ir,t,i_md,t_min,t_max,vt_max,md_steps
INTEGER :: i_s1,i_s2,ncv,nbin(*),index1,index2
REAL*8  :: dum,den,s1,s2
REAL*8  :: prob_2D(nbin(u),nbin(jj)),gridmin(*),griddif(*),cv(nr,ncv,*)
REAL*8  :: pcons(nr),kcons(nr)
LOGICAL :: max_step
CHARACTER(LEN=5)   :: mtd
CHARACTER(LEN=50)  :: filename_loc
CHARACTER(LEN=10)  :: code_name
CHARACTER(LEN=50)  :: filename(nr),filename_mtd(2,nr)

CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)
101 FORMAT (A,I1,A,I1,A,/,'Please Wait ...')
WRITE(6,101)"2D Probabilities along CV:",u,' and CV:',jj,' are being computed'
den = 0.0d0 ; t = jj
vt_max = t_max
!PRINT*,'t_min = ',t_min
DO ir = 1,nr
OPEN(22,FILE=filename(ir))
CALL get_steps(22,md_steps)
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)

prob_2D = 0

  DO i_md=t_min,t_max

index1 = nint((cv(ir,u,i_md)-gridmin(u))/griddif(u)) + 1
index2 = nint((cv(ir,t,i_md)-gridmin(t))/griddif(t)) + 1
!----------------------------------------------------------------------------
IF (index1 .gt. nbin(u) .or. index2 .gt. nbin(t)) THEN
PRINT*, "            ***************************************************"
PRINT*, "                ERROR !! ONE OF THE CV RANGE IS NOT CORRECT"
PRINT*, "            ***************************************************"
STOP ; ENDIF
!----------------------------------------------------------------------------
       prob_2D(index1,index2) = prob_2D(index1,index2) + 1.0
  END DO
!----------------------------------------------------------------------------
DO index1=1,nbin(u)
  DO index2=1,nbin(t)
    den=den+prob_2D(index1,index2)
  END DO
END DO
!----------------------------------------------------------------------------
dum = den*griddif(u)*griddif(t) ;  den = 1.D0/dum
CALL get_filename ('PROB_2D.dat_',filename_loc,ir)
OPEN(2,FILE=filename_loc,STATUS='unknown')
DO i_s1 = 1,nbin(u)
        s1 = DFLOAT(i_s1-1)*griddif(u)+gridmin(u)
  DO i_s2 = 1,nbin(t)
        s2 = DFLOAT(i_s2-1)*griddif(t)+gridmin(t)
        prob_2D(i_s1,i_s2) = prob_2D(i_s1,i_s2)*den
    WRITE(2,'(3E16.8)')s1,s2,prob_2D(i_s1,i_s2)
  END DO
    WRITE(2,*)
END DO
ENDDO
WRITE(*,'(A)')'Unbiased 2D distribution along US vs TEMP  written in PROB_2D.dat'
CLOSE(2) ; CLOSE(22)
END SUBROUTINE
END MODULE US_TEMP
