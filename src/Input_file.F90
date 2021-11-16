MODULE Input_file
CONTAINS
!
SUBROUTINE file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)
IMPLICIT NONE
INTEGER                 :: ir,nr,ios
REAL*8                  :: pcons(*),kcons(*)
REAL*8, PARAMETER       :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER       :: au_to_kcal = 627.51
REAL*8, PARAMETER       :: kj_to_kcal = 0.239006
CHARACTER(LEN=10)       :: code_name
CHARACTER(LEN=5)        :: mtd
CHARACTER(LEN=50)       :: filename(*),filename_mtd(2,*)

OPEN(10,FILE='input.inp',STATUS='old',IOSTAT=ios)
IF (ios .lt. 0) STOP "ERROR : FILE input.inp doesn't exist..!"
!=====================================================================!
DO ir = 1,nr
   IF(code_name .eq. 'CPMD') THEN
     READ(10,*)pcons(ir),kcons(ir)
     kcons(ir) = kcons(ir)*au_to_kcal
     READ(10,'(A)')filename(ir)
     IF(mtd .eq. 'y')THEN
       READ(10,'(A)')filename_mtd(1,ir)
       READ(10,'(A)')filename_mtd(2,ir)
     ENDIF
!---------------------------------------------------------------------!     
   ELSEIF (code_name .eq. 'PLUMED') THEN
     READ(10,*)pcons(ir),kcons(ir)
     kcons(ir) = kcons(ir)*kj_to_kcal
     READ(10,'(A)')filename(ir)
   IF(mtd .eq. 'y') THEN
     READ(10,'(A)')filename_mtd(1,ir)
   ENDIF
   ENDIF
!PRINT*,filename(ir)
ENDDO
!=====================================================================!
CLOSE(10)
END SUBROUTINE        
END MODULE       

