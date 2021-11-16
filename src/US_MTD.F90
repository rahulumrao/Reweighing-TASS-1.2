MODULE US_MTD
USE Input_file        
USE GetSteps
USE GetFileName
USE Error_msg
CONTAINS
!================================================== US vs MTD Probability ================================================!
!Compute Unbiased Distribution Along Umbrella Coordinate vs MTD coordinate
!
SUBROUTINE twoD_prob(ii,jj,u,m,nr,kt,w_cv,w_hill,ct,vbias,max_step,t_min,t_max,md_steps,den,prob_mtd &
                    & ,ncv,cv,nbin,gridmin,gridmax,griddif,norm,mtd,code_name)

IMPLICIT NONE
INTEGER :: ii,jj,u,m,i_md,i_mtd,t_min,t_max,index1,md_steps,i_s1,i_s2,ncv,nbin(*),w_cv,w_hill,indx(2)
INTEGER :: ir,nr,vt_max,ios
REAL*8  :: dum,den,s1,s2,kt,norm(nr)
REAL*8  :: prob_mtd(nr,nbin(u),nbin(m)),prob_2D(nbin(u),nbin(m))
REAL*8  :: gridmin(*),gridmax(*),griddif(*),cv(nr,ncv,*),vbias(nr,*),ct(nr,*)
LOGICAL :: max_step
REAL*8  :: pcons(nr),kcons(nr)
CHARACTER(LEN=10)  :: code_name
CHARACTER(LEN=5)   :: mtd
CHARACTER(LEN=50)  :: filename_loc
CHARACTER(LEN=50)  :: filename(nr),filename_mtd(2,nr)
!--------------------------------------------------------------------------------------------------
CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)
!PRINT*,filename(1),code_name
!=====================================================================!
IF (ii .eq. u .and. jj .eq. m ) THEN
vt_max = t_max        
DO ir = 1,nr
OPEN(22,FILE=filename(ir))
CALL get_steps(22,md_steps)
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)
!WRITE(6,'(A,I5,2X,A,I5)')'Minumum MD Steps =>',t_min,'; Maximum MD Steps =>',t_max
!--------------------------------------------------------------------------------------------------
DO i_md=t_min,t_max

 indx(u) = nint((cv(ir,u,i_md)-gridmin(u))/griddif(u)) + 1
 indx(m) = nint((cv(ir,m,i_md)-gridmin(m))/griddif(m)) + 1
!----------------------------------------------------------------------------
IF (indx(u) .gt. nbin(u) .or. indx(m) .gt. nbin(m)) THEN
PRINT*, "            ***************************************************"
PRINT*, "                ERROR !! ONE OF THE CV RANGE IS NOT CORRECT"
PRINT*, "            ***************************************************"
STOP ; ENDIF
END DO
!----------------------------------------------------------------------------
CALL get_filename ('PROB_2D.dat_',filename_loc,ir)   ! get individual file to write probability for every umbrella
OPEN(2,FILE=filename_loc,STATUS='unknown')
DO i_s1 = 1,nbin(u)
          s1 = DFLOAT(i_s1-1)*griddif(u)+gridmin(u)
  DO i_s2 = 1,nbin(m)
          s2 = DFLOAT(i_s2-1)*griddif(m)+gridmin(m)
          prob_2D(i_s1,i_s2) = prob_2D(i_s1,i_s2)*den
    WRITE(2,'(3E16.8)')s1,s2,prob_mtd(ir,i_s1,i_s2)*norm(ir)
  END DO
    WRITE(2,*)
END DO
ENDDO
WRITE(*,'(A)')'Unbiased 2D distribution along US vs MTD  written in Pu_2D.dat'
ENDIF
CLOSE(2) ; CLOSE(22)
END SUBROUTINE
END MODULE US_MTD
