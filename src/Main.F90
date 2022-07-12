PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------------------------------------------------------!
!FORTRAN PROGRAM WRITTEN TO COMPUTE UNBIASED DISTRIBUTION (1D AND 2D)FROM TASS SIMULATION ALONG USER DEFINED          !
!COLLECTIVE COORDINATES..                                                                                             !
!ORIGINAL CODE WRITTEN BY Shalini Awasthi (ashalini@iitk.ac.in)                                                       !
!MODIFIED BY Rahul Verma (vrahul@iitk.ac.in)                                                                          !
!                                                                                                                     !
!kt0 = System Temeprature ; kt = Extended CV Temperature ; bias_fact = Biased Factor for MTD ; ct = ct Factor         !
!t_min = Minimum MD steps ; t_max = Maximum MD steps ; narg = Argumets ; ncv = Numer of CV ; v_baias = Total MTD bias !
!cv_mtd = m = MTD CV index ; u_cv = u = Umbrella CV index ; t_cv = t = Temeprature CV index                           !
!UCV = U cv argument ; MTD = MTD cv argument ; Prob_nD = Dimension of Unbiased Probability                            !
!CV_num = Probability is the Dimension of ; pfrqMD = Print Frequency argument ; w_cv Print Frequency in cvmdck File   !
!dtMTD  = Print Frequency argument for MTD bias added ; w_hill = Print Frequency in colvar File                       !
!gridmin = Minimum Grid Size = gridmax = Maximum Grid Size ; griddif = Grid Difference                                !
!width = Hill Width of Gaussian Bias in MTD ; ht = Hill Height of Gaussian Bias in MTD ; ht = MTD CV Displacement     !
!kb = Boltzman Constant in A.U. ; prob = 1D Probability ; prob_2D = 2D Probability ; prob_mtd = mtd_biased Probability!
!pmf = MeanForce based Reweighting ; ProbT = Generate Unbiased Probability ; spline = B-Spline Interpolation          !
!filename = cv Filename (COLVAR for PLUMED and cvmdck_mtd for CPMD)                                                   !                       
!filename_mtd = mtd data Filename (HILLS for PLUMED and parvar_mtd & colvar_mtd for CPMD)                             !
!---------------------------------------------------------------------------------------------------------------------!
USE ansi_colors
USE GetSteps
USE GetFileName
USE Input_file
USE Error_msg
USE MTD_Unbais
USE MTD_Potential
USE US_Prob
USE US_MTD
USE US_TEMP
USE MeanForce
USE Error_estimation
!USE wham_code
USE B_Spline
IMPLICIT NONE
REAL*8              :: dummy11,den,kt0,kt,bias_fact
REAL*8, ALLOCATABLE :: prob(:),fes(:),pcons(:),kcons(:)
REAL*8, ALLOCATABLE :: dummy(:,:,:),cv(:,:,:),prob_2D(:,:),prob_mtd(:,:,:)
REAL*8, ALLOCATABLE :: gridmin(:),gridmax(:),griddif(:),vbias(:,:),ct(:,:),norm(:)
INTEGER,ALLOCATABLE :: nbin(:),indx(:),t(:),t_cv(:),cv_num(:)
INTEGER :: md_steps,mtd_steps,dummy1,i,j,k,t_min,t_max,narg
INTEGER :: i_md,ncv,w_hill,w_cv,prob_nD,cv_mtd,cv_us!,cv_num(3)
INTEGER :: ii,jj,kk,u,m,ir,nr,ios,nblock
LOGICAL :: pmf,probT,spline,inpgrid,read_ct,read_vbias,max_step,stat_error
CHARACTER(LEN=5)  :: mtd,tool
CHARACTER(LEN=10) :: code_name
CHARACTER(LEN=120):: arg 
CHARACTER(LEN=50) :: filename_loc
CHARACTER(LEN=80) :: filenm
CHARACTER(LEN=50),ALLOCATABLE :: filename(:),filename_mtd(:,:)
REAL*8, PARAMETER :: kb = 1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51 
REAL*8, PARAMETER :: kj_to_kcal = 0.239006 
!-----------------------------------------------!
md_steps    = 9999999  ; mtd_steps  = 999999    !
kt0         = 300.D0   ; kt         = 300.D0    !
t_min       = 1        ; t_max      = 7000      !
w_hill      = 0        ; w_cv       = 0         !
pmf         = .FALSE.  ; inpgrid    = .FALSE.   !
read_ct     = .FALSE.  ; read_vbias = .FALSE.   !
mtd         = 'n'      ; tool       = 'pmf'     !
probT       = .FALSE.  ; spline     = .FALSE.   !
narg        = IARGC()  ; bias_fact  = 1500.D0   !
max_step    = .FALSE.  ; stat_error = .FALSE.   !
!!----------------------------------------------!
!PRINT*,"!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!"
!PRINT*,"ARGUMENTS IN THE RUN FILE ARE CASE SENSITIVE"
!!----------------------------------------------------!
k = 0
DO i = 1,999
READ(5,'(A)',ERR=999,END=999)filenm(1:80)

IF(INDEX(filenm,'SYSTEM TEMP') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) kt0 !sys_t
ELSEIF(INDEX(filenm,'CV TEMP') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) kt !cv_temp
ELSEIF(INDEX(filenm,'CODE NAME') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) code_name
        IF(code_name .ne. 'PLUMED') THEN
        IF(code_name .ne. 'CPMD') THEN
        CALL prog_error
        ENDIF ; ENDIF
ELSEIF(INDEX(filenm,'TMIN') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) t_min
ELSEIF(INDEX(filenm,'TMAX') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) t_max
        max_step=.TRUE.
ELSEIF(INDEX(filenm,'NUMBER OF CV') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) ncv
!------------------------------------------------!
ALLOCATE (cv_num(ncv))
ALLOCATE (gridmin(ncv),gridmax(ncv),griddif(ncv))
!------------------------------------------------!
ELSEIF(INDEX(filenm,'UCV COLUMN') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) cv_us
ELSEIF(INDEX(filenm,'MTD ENABLED') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) mtd
ELSEIF(INDEX(filenm,'MTD CV COLUMN') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) cv_mtd
ELSEIF(INDEX(filenm,'MTD BIAS FACTOR') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) bias_fact
ELSEIF(INDEX(filenm,'REWEIGHTING TOOL') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) tool
        IF (tool .eq. 'pmf') pmf=.TRUE.
        IF (tool .eq. 'prob') probT=.TRUE.
ELSEIF(INDEX(filenm,'PROBABILITY DIMENSION') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) prob_nd
        IF (prob_nD .gt. 3) STOP  "MORE THAN 3 DIMENSION IS NOT IMPLIMENTED"
ELSEIF(INDEX(filenm,'PROBABILITY CV INDEX') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) cv_num(1:prob_nd)
ELSEIF(INDEX(filenm,'NUMBER OF UMBRELLA') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) nr
ELSEIF(INDEX(filenm,'GRIDS') .ne. 0) THEN
        DO j = 1,ncv
        READ(5,*,ERR=999,END=999) gridmin(j),gridmax(j),griddif(j)
        k = k + 1
        ENDDO
ELSEIF(INDEX(filenm,'CV PRINT FREQUENCY') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) w_cv
ELSEIF(INDEX(filenm,'MTD PRINT FREQUENCY') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) w_hill
ELSEIF(INDEX(filenm,'STATISTICAL ERROR BLOCK SIZE') .ne. 0) THEN
        READ(5,*,ERR=999,END=999) nblock
        stat_error=.TRUE.
ELSEIF(INDEX(filenm,'B-SPLINE INTERPOLATION') .ne. 0) THEN
        spline=.TRUE.
ELSEIF(INDEX(filenm,'READ CT') .ne. 0) THEN
        read_ct=.TRUE.
ELSEIF(INDEX(filenm,'READ VBIAS') .ne. 0) THEN
        read_vbias=.TRUE.
        
ENDIF
END DO
999 CONTINUE
IF (k .lt. ncv) call grid_error ! IF ncv is not eq to GRID size in input
!=======================================================================!
IF (code_name .eq. 'CPMD')   OPEN(11,FILE='cvmdck_mtd',STATUS='unknown')
IF (code_name .eq. 'PLUMED') OPEN(11,FILE='COLVAR',STATUS='unknown')
!=======================================================================!
IF (mtd .eq. 'y') THEN
IF (code_name .eq. 'CPMD')   OPEN(12,FILE='parvar_mtd',STATUS='unknown')
IF (code_name .eq. 'CPMD')  OPEN(13,FILE='colvar_mtd',STATUS='unknown')
IF (code_name .eq. 'PLUMED')OPEN(13,FILE='HILLS',STATUS='unknown')
ENDIF
!=======================================================================!
IF (mtd .eq. 'y' .and. cv_mtd .eq. 0) STOP "***ERROR !! PLEASE SPECIFY METADYNAMICS CV COORDINATE INDEX"
WRITE(*,'(A)')'!------------------------------------------------------------------------------------'
     WRITE(*,'(A,2X,A)')color('! MOLECULAR DYNAMICS PACKAGE FOR SAMPLING      =',c_green),color(code_name, c_red)
     WRITE(*,'(A,2X,I2)')color('! No. of Umbrella Windows                      = ',c_green),nr
  !IF(probT)   WRITE(*,'(A,I10)')'! No: of MD  steps                            =',md_steps
     WRITE(*,'(A,I8)')color('! No. of min MD  steps                         =',c_green),t_min
  IF (max_step) WRITE(*,'(A,I9)')color('! No. of max MD  steps                         =',c_green),t_max
     WRITE(*,'(A,2X,I3)')color('! No. of CV                                    =',c_green),ncv
     WRITE(*,'(A,2X,I3)')color('! Umbrella CV Column                           =',c_green),cv_us
  IF (mtd .eq. 'n') THEN
     WRITE(*,'(A,4X,A)')color('! Metadynamics is Enabled                      =',c_green),'NO'
  ELSEIF(mtd .eq. 'y') THEN
     WRITE(*,'(A,5X,A)')color('! Metadynamics is Enabled                      =',c_green),'YES'
     WRITE(*,'(A,2X,I5)')color('! Metadynamics CV Column                       =',c_green),cv_mtd
  ENDIF
     WRITE(*,'(A,F10.2)')color('! Physical System Temperature T0 (K)           =',c_green),kt0
     WRITE(*,'(A,F11.2)')color('! CV Temperature T (K)                         =',c_green),kt
  IF (mtd .eq. 'y')  WRITE(*,'(A,F11.2)')color('! Bias Factor (K)                              =',c_green),bias_fact
     WRITE(*,'(A,I6)')color('! Print Freq. in cv file                       =',c_green),w_cv
     IF(mtd.eq.'y') WRITE(*,'(A,I7)')color('! Freq. of Hill Update                         =',c_green),w_hill
  IF (pmf) THEN
     WRITE(*,'(A,2X,A)')color('! ####### FREE ENERGY WILL BE CALCULATED USING MEAN FORCE METHOD #######',c_red)
  ELSEIF (probT)THEN
     WRITE(*,'(A,I6)')color('! Dimension of Probability                     =',c_green),prob_nD
     WRITE(*,'(A,3I7)')color('! Probability Index                            =',c_green),cv_num(1:prob_nD)
     WRITE(*,'(A,2X,A)')color('!####### PROBABILITIES WILL BE GENERATED THROUGH UNBIASING #######',c_red)
  ENDIF
WRITE(*,'(A)')'!------------------------------------------------------------------------------------'
ii = cv_num(1); jj = cv_num(2) ; kk = cv_num(3)
u  = cv_us    ; m  = cv_mtd  !! #### m = METADYNAMICS CV INDEX !! #### u = UMBRELLA CV INDEX
!===========================================================================================================!
IF(probT) THEN
WRITE(*,'(A85)')'=========================================================================================='! 
IF (ii .ne. u .and. jj .ne. u .and. kk .ne. u) THEN
WRITE(*,'(10X,A)') "!!SORRY!! PLEASE SPECIFY UMBRELLA CV INDEX PROPERLY" 
WRITE(*,'(A85)')'=========================================================================================='!
STOP ; ENDIF
IF (mtd .eq. 'y' .and. cv_mtd .eq. cv_us) THEN
!IF (ii .ne. m .and. jj .ne. m .and. kk .ne. m) THEN 
WRITE(*,'(10X,A)')"METADYNAMICS ENABLED, CAN'T UNBIAS PROBABILITY WITHOUT UNBIASING 'MTD'"
WRITE(*,'(A85)')'=========================================================================================='!
STOP
ENDIF ; ENDIF 
!-----------------------------------------------------------------------------------------------------!
ALLOCATE(filename(nr))          ; ALLOCATE(filename_mtd(2,nr))
ALLOCATE(cv(nr,ncv,md_steps))   ; ALLOCATE(dummy(nr,ncv,md_steps))
ALLOCATE(nbin(ncv))             ; ALLOCATE(vbias(nr,md_steps))
ALLOCATE(ct(nr,mtd_steps))      ; ALLOCATE(t(ncv))
ALLOCATE(kcons(nr))             ; ALLOCATE(pcons(nr))
ALLOCATE(t_cv(ncv))             ; ALLOCATE(norm(nr))
!-----------------------------------------------------------------------------------------------------!
OPEN(10,FILE='input.inp',STATUS='old',IOSTAT=ios,ACTION='read')
IF (ios .ne. 0) STOP "ERROR : FILE input.inp doesn't exist..!"

!OPEN(30,FILE='data_md.rst',FORM="UNFORMATTED")
101 FORMAT (I10,10F16.6)
102 FORMAT (F10.6,10F16.6)

DO ir = 1,nr
!-----------------------------------------------------------------------------------------------------!
   IF(code_name .eq. 'CPMD') THEN
     READ(10,*)pcons(ir),kcons(ir)
      kcons(ir)=kcons(ir)*au_to_kcal
      READ(10,'(A)')filename(ir)
     IF(mtd .eq. 'y')THEN
       READ(10,'(A)')filename_mtd(1,ir)
       READ(10,'(A)')filename_mtd(2,ir)
     ENDIF
CALL get_filename('cv.dat_',filename_loc,ir)
OPEN(14,FILE=filename_loc,STATUS='unknown')
OPEN(11,FILE=filename(ir),STATUS='old',IOSTAT=ios)
CALL get_steps(11,md_steps)
DO i_md=1,md_steps
   READ(11,*)dummy1,dummy1,(dummy(ir,j,i_md), j=1,ncv),(cv(ir,j,i_md) ,j=1,ncv)
!  WRITE(30)(cv(ir,j,i_md) ,j=1,ncv) ! Writing restart file
!   WRITE(*,*)dummy1,(dummy(ir,j,i_md), j=1,ncv),(cv(ir,j,i_md) ,j=1,ncv)
!-------------------------------------------------------------------------
DO j = 1,ncv
IF ( cv(ir,j,i_md) .lt. gridmin(j)) CALL cv_error (j)
IF ( cv(ir,j,i_md) .gt. gridmax(j)) CALL cv_error (j)
ENDDO
!-------------------------------------------------------------------------
   WRITE(14,101)dummy1,(cv(ir,j,i_md) ,j=1,ncv)
END DO
WRITE(6,'(A,I2,3X,A,I10)')'! No. of MD STEPS in umb:',ir,'=',md_steps
CLOSE(11);CLOSE(14)
!-------------------------------------------------------------------------
   ELSEIF (code_name .eq. 'PLUMED') THEN
     READ(10,*)pcons(ir),kcons(ir)
     kcons(ir)=kcons(ir)*kj_to_kcal
     READ(10,'(A)')filename(ir)
   IF(mtd .eq. 'y') THEN
     READ(10,'(A)')filename_mtd(1,ir)
   ENDIF

CALL get_filename('cv.dat_',filename_loc,ir)
OPEN(14,FILE=filename_loc,STATUS='unknown')
OPEN(11,FILE=filename(ir),STATUS='old',IOSTAT=ios)
CALL get_steps(11,md_steps)
DO i_md=1,md_steps
!   READ(11,*)dummy11,dummy11,cv(ir,u,i_md),dummy11,cv(ir,m,i_md)
   READ(11,*)dummy11,(cv(ir,j,i_md) ,j=1,ncv)
!-------------------------------------------------------------------------
DO j = 1,ncv
IF ( cv(ir,j,i_md) .lt. gridmin(j)) CALL cv_error (j)
IF ( cv(ir,j,i_md) .gt. gridmax(j)) CALL cv_error (j)
ENDDO
!-------------------------------------------------------------------------
!   WRITE(*,*)dummy11,(cv(j,i_md) ,j=1,ncv)
   WRITE(14,102)dummy11,(cv(ir,j,i_md) ,j=1,ncv)
END DO
WRITE(6,'(A,I2,3X,A,I10)')'! No. of MD STEPS in umb:',ir,'=',md_steps
CLOSE(11);CLOSE(14)
ENDIF
ENDDO
CLOSE(10)
!-----------------------------------------------------------------------------------------------------!
WRITE(*,'(A85)')'=========================================================================================='!
  DO i = 1,ncv
   nbin(i) = NINT((gridmax(i)-gridmin(i))/griddif(i)) + 1
  ENDDO

j = 0
WRITE(*,'(9X,4A9)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
WRITE(*,'(A10,3F8.4,I10)')'US   COORD:', gridmin(u),gridmax(u),griddif(u),nbin(u)
IF (mtd .eq. 'y') THEN
WRITE(*,'(A10,3F8.4,I10)')'MTD  COORD:', gridmin(m),gridmax(m),griddif(m),nbin(m)
ENDIF
DO i = 1,ncv ; IF (i .ne. u .and. i .ne. m) THEN
WRITE(*,'(A10,3F8.4,2I10)')'TASS COORD:', gridmin(i),gridmax(i),griddif(i),nbin(i)
j = j + 1 ; t_cv(j) = i
ENDIF ; ENDDO
WRITE(*,'(A85)')'=========================================================================================='!
t(1:j) = t_cv(1:j) !; PRINT*,t(1:j)    !# t_cv TASS CV INDEX
!n1 = nbin(1) ;n2 = nbin(2) !;n3 = nbin(3) !;n4 = nbin(4) 
!---------------------------------------------------------------------------------------------------------------------------!
IF (stat_error) CALL statistical_error(cv,nr,u,ncv,mtd,code_name,nblock)
!---------------------------------------------------------------------------------------------------------------------------!
ALLOCATE(prob((nbin(u))))
IF (prob_nd .eq. 2 .and. mtd .eq. 'y') THEN
ALLOCATE(prob_mtd(nr,nbin(u),nbin(m)))
ELSEIF (prob_nd .eq. 2 .and. mtd .ne. 'y') THEN
m = t(1)
ALLOCATE(prob_2D(nbin(u),nbin(m)))
ENDIF
!prob = 0.d0 ; prob_2D = 0.d0 ; prob_mtd = 0.d0
IF (jj .eq. 0 .and. kk .eq. 0) jj = 1 ; kk = 1
!ALLOCATE(prob_3D(nbin(ii),nbin(jj),nbin(kk)))
<<<<<<< HEAD
CLOSE (11) ; CLOSE(12) ; CLOSE(13)
=======
CLOSE(11) ; CLOSE(12) ; CLOSE(13)
>>>>>>> 022596b624f83121171634353808cf574fd6d83e
!---------------------------------------------------------------------------------------------------------------------------!
IF(mtd .eq. 'y') THEN

CALL mtd_unbiased(au_to_kcal,bias_fact,kt,kb,md_steps,mtd_steps,w_cv,w_hill &
           & ,ncv,t_min,t_max,gridmin,gridmax,griddif,vbias,ct,nbin,m,mtd,code_name,nr)

CALL mtd_pot(md_steps,mtd_steps,w_cv,w_hill,t_min,t_max,gridmin,gridmax,griddif,vbias, &
           & ct,m,u,ncv,kt,nbin,cv,den,prob_mtd,norm,ir,nr,mtd,code_name)
ENDIF
IF (pmf) THEN
ALLOCATE (fes(nr))
  CALL mean_force(max_step,u,cv_num,prob_nD,ncv,cv,nr,kt,nbin,t_min,t_max,pcons,kcons, &
                 & gridmin,gridmax,griddif,fes,mtd,w_cv,w_hill,ct,vbias,code_name)
  IF(spline) CALL bspline(nr,pcons,fes)
DEALLOCATE(pcons,fes)
ELSE
!---------------------------------------------------------------------------------------------------------------------------!
!Computing 1-Dimensional Probability along Umbrella Coordinate
IF (Prob_nD .eq. 1) THEN
CALL oneD_prob(ii,nr,u,m,mtd,max_step,t_min,t_max,md_steps,den,prob,prob_mtd,ncv,cv,nbin &
               & ,gridmin,griddif,norm,code_name)
DEALLOCATE(prob)
!---------------------------------------------------------------------------------------------------------------------------!
!Computing 2-Dimensional Probability along Umbrella/MTD (if MTD is enabled) or US/TEMP
ELSEIF (Prob_nd .eq. 2) THEN
IF (mtd .eq. 'y') THEN
CALL twoD_prob(ii,jj,u,m,nr,kt,w_cv,w_hill,ct,vbias,max_step,t_min,t_max,md_steps,den,prob_mtd, &
              & ncv,cv,nbin,gridmin,gridmax,griddif,norm,mtd,code_name)

ELSEIF (mtd .eq. 'n') THEN
CALL twoD_temp_prob(jj,u,nr,max_step,t_min,t_max,md_steps,prob_2D,ncv,cv,nbin, &
              & gridmin,griddif,mtd,code_name)
!---------------------------------------------------------------------------------------------------------------------------!
DEALLOCATE(prob_2D) ; DEALLOCATE(prob_mtd)
DEALLOCATE(ct,cv,vbias,nbin)
ENDIF ; ENDIF ; ENDIF
!===========================================================================================================================!
END PROGRAM WSMTD_rw_2D
