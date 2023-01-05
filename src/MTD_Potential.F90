MODULE MTD_Potential
!======================================================================================================!
! This Subroutine calculates the probability and stores in 'prob_mtd(ir,indx(u),indx(m))'
! which is then used in other routines to get the biased prob along other dimensions.
!======================================================================================================!
USE Input_file
USE GetSteps
USE Error_msg
CONTAINS
SUBROUTINE mtd_pot(md_steps,mtd_steps,w_cv,w_hill,t_min,t_max,max_step,gridmin,gridmax,griddif,vbias, &
                 ct,m,u,ncv,kt,nbin,cv,den,prob_mtd,norm,ir,nr,mtd,code_name)

IMPLICIT NONE
INTEGER           :: i_md,md_steps,i_mtd,mtd_steps,t_min,t_max,m,u,w_cv,w_hill
INTEGER           :: ncv,indx(ncv),nbin(*),j,ir,nr,vt_max
REAL*8            :: kt,den,dum,ct(nr,*),gridmin(*),gridmax(*),griddif(*),vbias(nr,*)
REAL*8            :: cv(nr,ncv,*),prob_mtd(nr,nbin(u),nbin(m))
REAL*8            :: pcons(nr),kcons(nr),norm(nr)
CHARACTER(LEN=5)  :: mtd
CHARACTER(LEN=10) :: code_name
CHARACTER(LEN=50) :: filename(nr),filename_mtd(2,nr)
LOGICAL           :: max_step
!======================================================================================================!

CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)

!calculate prob (unbiased from MTD potential)
prob_mtd = 0.0
DO ir = 1,nr
OPEN(22,FILE=filename(ir))
OPEN(23,FILE=filename_mtd(1,ir))
CALL get_steps(22,md_steps)
CALL get_steps(23,mtd_steps)
vt_max = t_max 
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)   ! get t_max
den = 0.0 ; dum = 0.0 
DO i_md=t_min,t_max
!-------------------------------------------------------------------------
DO j = 1,ncv ! error mesaage if the cv range is not correct
IF ( cv(ir,j,i_md) .lt. gridmin(j)) CALL cv_error (j,'min')
IF ( cv(ir,j,i_md) .gt. gridmax(j)) CALL cv_error (j,'max')
ENDDO
!-------------------------------------------------------------------------
     indx(u) = nint((cv(ir,u,i_md)-gridmin(u))/griddif(u)) + 1
     indx(m) = nint((cv(ir,m,i_md)-gridmin(m))/griddif(m)) + 1
       i_mtd = ((i_md - 1)*w_cv/w_hill) + 1
         dum = vbias(ir,i_md) - ct(ir,i_mtd)
         dum = dexp(dum/kt)
 prob_mtd(ir,indx(u),indx(m)) = prob_mtd(ir,indx(u),indx(m)) + dum
         den = den + dum
!PRINT*,ir,prob_mtd(ir,indx(u),indx(m))
END DO
norm(ir) = 1.0d0/(den*griddif(u)*griddif(m))  ! normalization 
END DO
CLOSE(22) ; CLOSE(23)
END SUBROUTINE mtd_pot
END MODULE MTD_Potential
!======================================================================================================!
