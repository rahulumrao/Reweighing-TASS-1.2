!*************************************************************************************************!
! THIS IS WHAM CODE TO COMPUTE FREE ENERGIES FROM PRECONPUTED PROBABILITIES.
! IT IS CAPABLE OF PATCHING 1D/2D FREE ENERGY DEPENDING WHAT ARE THE DIMENSION OF PROBABILITIES
! WRITTEN BY :: Shalini Awasthi
! MODULAR CODE WRITTEN by :: Rahul Verma
!*************************************************************************************************!
PROGRAM fes_calc
IMPLICIT NONE
INTEGER :: i, iter, ncv, umbr_n, nbin1, nbin2
INTEGER :: i_umbr, i_s1, i_s2
INTEGER :: rank, gleng1_min,gleng1_max,gleng2,ngrid
INTEGER, ALLOCATABLE :: nbin(:)
INTEGER, ALLOCATABLE :: nmax(:)
REAL*8 :: kt, toler, dummy, cnvg
REAL*8 :: dummy_1,dummy_2,dummy_3
REAL*8, ALLOCATABLE :: umbr_mean(:), umbr_k(:),a(:)
REAL*8, ALLOCATABLE :: grid0(:,:), v(:), prob(:,:),biased_prob(:,:,:),grid(:,:)
REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51
LOGICAL :: parent
CHARACTER (LEN=50) :: cvfile,outputfile,f1,f3

! INITIALIZE MPI
CALL MPI_Start
CALL Set_Parent(parent)
!
IF(PARENT)THEN
  OPEN (UNIT=1, FILE='whaminput',STATUS='old') 
  READ(1,*)toler, umbr_n        ! read : tolerance  , no. of umbrella 
  READ(1,*)ncv, kt              ! read : no. of cvs , temperature at which prob are generated
END IF
  kt=kb*kt
!
ALLOCATE(nbin(2))
!broadcast ncv
CALL IBcast(ncv,1)
CALL IBcast(umbr_n,1)
CALL RBcast(kt,1)

!allocate grid0
allocate(grid(3,ncv))
allocate(grid0(3,ncv))
allocate(a(umbr_n))
allocate(umbr_mean(umbr_n))
allocate(umbr_k(umbr_n))
allocate(nmax(umbr_n))
a=1.d0

IF (parent) THEN
 DO i=1,ncv
    READ(1,*)grid0(1:3,i)       ! read : grid_min, grid_max, grid_width
    WRITE(*,'(i10,3f16.6)')i,grid0(1:3,i)
 END DO
END IF

IF (parent) THEN
  IF (ncv.le.2) THEN
   DO i = 1, ncv        
     nbin(i)=nint((grid0(2,i)-grid0(1,i))/grid0(3,i))+1  ! compute : total bins
     IF (i .eq. 1) nbin1=nbin(1) ; nbin2=1               ! if probs are 1D then setting nbin2 to 1 for further
     IF (i .eq. 2) nbin1=nbin(1) ; nbin2=nbin(2)
   ENDDO
  ELSE IF (ncv.ge.3)THEN
STOP '3 or more CVs not implemented'
  END IF
END IF

!broadcast grids and bin info
CALL RBcast(grid0,3*ncv)
CALL IBcast(nbin1,1)
CALL IBcast(nbin2,1)
allocate(biased_prob(nbin1,nbin2,umbr_n))
allocate(prob(nbin1,nbin2))
! READING PROBABILITY FROM EACH UMBRELLA WINDOW
IF (parent) THEN
   DO i_umbr=1, umbr_n 
      WRITE(*,*) 'umbrella simulation #', i_umbr
      !reads force constant, r0, max points in that umbrella
      READ(1,*) umbr_mean(i_umbr),umbr_k(i_umbr),nmax(i_umbr)  ! reads : umbrella mean ,force constant(kcal/mol), MD steps in umbrella
      !reads probability file name
      IF (ncv .eq. 1 ) CALL get_filename('PROB.dat_',cvfile,i_umbr)
      IF (ncv .eq. 2 ) CALL get_filename('PROB_2D.dat_',cvfile,i_umbr)
      WRITE(*,'(A,1x,A)')'PROB FILE ::',cvfile
      f1 = cvfile 
      OPEN(UNIT=3, FILE=f1, STATUS='old' )
      DO i_s1=1,nbin1 !US
        IF(ncv .eq. 1) THEN
         READ(3,*)dummy,biased_prob(i_s1,1,i_umbr)
        ELSEIF(ncv .eq. 2) THEN
         DO i_s2=1,nbin2 !MTD
         READ(3,*)dummy,dummy,biased_prob(i_s1,i_s2,i_umbr)
         END DO
        ENDIF
      END DO
   ENDDO
END IF

call RBcast(toler,1)
call RBcast(umbr_mean,umbr_n)
call RBcast(umbr_k,umbr_n)
call IBcast(nmax,umbr_n)
call RBcast(biased_prob,nbin1*nbin2*umbr_n)

! PERFORMS WHAM
 IF (parent)  WRITE(*,*) 'wham begins'
! DISTRIBUTE DATA TO EVERY PROCESSOR
 CALL DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)

 IF (ncv .eq. 2) gleng2=nint((grid(2,2)-grid(1,2))/grid(3,2))+1

 WRITE(*,*)'new_grid', gleng1_min, gleng1_max,gleng2, rank, grid(1,1)
  iter=0
  scf_loop : do
       iter=iter+1
       IF (ncv .eq. 1) THEN
       call wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,1)
       ELSEIF (ncv .eq. 2) THEN
       call wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
       ENDIF
       IF (parent) WRITE(*,*) 'iteration #', iter, 'convergence =',cnvg
       IF (mod(iter,100) .eq. 0 ) then
       call print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
       ENDIF
       IF((cnvg.lt.toler).or.(iter.ge.200000))then
       IF (parent) write(*,*)'** convergence achieved **'
       EXIT scf_loop
       END IF
   END DO scf_loop

!prints the free energy.
call print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
call MPI_Stop
end program 

!***********************************************!


 subroutine print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
!prints pmf 
 implicit none
 integer:: i_s1, i_s2
 integer:: nbin1, nbin2, ncv
 real*8 :: prob(nbin1,nbin2), grid0(3,ncv)
 real*8 :: s1, s2, kt, dum
 logical:: parent
 real*8, allocatable :: prob1(:,:)
 character (len=50):: f2
! 
 ALLOCATE(prob1(nbin1,nbin2))
 prob1=0.0
 CALL GlobSumR(prob,prob1,nbin1*nbin2)
 IF (parent)THEN
 f2= 'free_energy_wham'
 OPEN( unit =7 , file = f2, status =  'unknown' )
 IF(ncv .eq. 2 ) THEN
 DO i_s1=1,nbin1 !US cv
  s1=DFLOAT(i_s1-1)*grid0(3,1)+grid0(1,1)
  DO i_s2=1,nbin2 !MTD cv
     s2=DFLOAT(i_s2-1)*grid0(3,2)+grid0(1,2)
     dum= -kt*DLOG(prob1(i_s1,i_s2))
     WRITE(7,'(3E16.8)') s1, s2, dum
  ENDDO
  WRITE(7,*)
 ENDDO
 ELSEIF(ncv .eq. 1 ) THEN
 DO i_s1=1,nbin1 !US cv
  s1=DFLOAT(i_s1-1)*grid0(3,1)+grid0(1,1) 
  i_s2=1
     dum= -kt*DLOG(prob1(i_s1,i_s2))
     WRITE(7,'(3E16.8)') s1, dum
  ENDDO
 ENDIF
 WRITE(*,*) 'free energy written in ',f2
 CLOSE(7)
 END IF
 endsubroutine
!**************************************************************************************************!


subroutine wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
!performs wham scf.
implicit none
integer:: i_s1, i_s2, i_umbr, nbin1, nbin2, umbr_n,ncv
integer:: nmax(umbr_n) 
real*8 :: umbr_k(umbr_n), umbr_mean(umbr_n),dum
real*8 :: num, den, dummy_v, dummy_s1, avg, del_s1, kt, dummy, cnvg
real*8 :: prob(nbin1,nbin2),a(umbr_n)
real*8 :: grid0(3,ncv), biased_prob(nbin1,nbin2,umbr_n),grid(3,ncv)
integer:: rank, gleng1_min,gleng1_max,gleng2,ngrid
real*8,allocatable :: dummy_a1(:),dummy_a(:) 

allocate(dummy_a(umbr_n))
allocate (dummy_a1(umbr_n) )

dummy_a = 0.0d0

!calculates probability at each grid_point.
do i_s1 =gleng1_min,gleng1_max !over US cv
do i_s2 =1,gleng2 !over MTD cv
   num = 0.0d0
   den = 0.0d0
   dummy_s1 = grid0(1,1)+dfloat(i_s1-1)*grid0(3,1)

   !calculates probability.
   do i_umbr=1,umbr_n
      del_s1=dummy_s1-umbr_mean(i_umbr)
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      num=num+dfloat(nmax(i_umbr))*biased_prob(i_s1,i_s2,i_umbr)
      den=den+dfloat(nmax(i_umbr))*a(i_umbr)*dummy_v
   enddo

   prob(i_s1,i_s2)=num/den
   if(prob(i_s1,i_s2).ne.prob(i_s1,i_s2)) prob(i_s1,i_s2)=1.0D-16 !remove NaN
   if(prob(i_s1,i_s2)+1.eq.prob(i_s1,i_s2)) prob(i_s1,i_s2)=1.0D-16 !remove infinity

!calculate a.
   IF (ncv .eq. 1) dum=grid0(3,1)
   IF (ncv .eq. 2) dum=grid0(3,1)*grid0(3,2)
   do i_umbr=1,umbr_n
      del_s1=dummy_s1 - umbr_mean(i_umbr)
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      dummy_a(i_umbr)=dummy_a(i_umbr) +dum*dummy_v*prob(i_s1,i_s2)
   enddo
enddo !end of MTD cv loop
enddo !end of US cv loop

!=======================================================================================

CALL GlobSumR(dummy_a,dummy_a1,umbr_n)
do i_umbr=1,umbr_n
dummy_a(i_umbr)=dummy_a1(i_umbr)
end do
!========================================================================================

!finds convergence and update a.
 avg = 0.0d0
 cnvg = 0.0d0
 do i_umbr=1,umbr_n
 dummy_a(i_umbr) = 1.0d0/dummy_a(i_umbr)
 cnvg = cnvg + dabs(dlog(dummy_a(i_umbr))-dlog(a(i_umbr)))
 a(i_umbr) = dummy_a(i_umbr)
 enddo
 cnvg = kt*cnvg
 end subroutine
!**************************************************************************************************!



SUBROUTINE DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)
!Distribute X grid over processors by mapping grid0 to grid 
IMPLICIT NONE
INTEGER :: ncv,rank
REAL*8 :: grid0(3,ncv), grid(3,ncv)
!
INTEGER :: i,ncpu,icpu,ngrids,ngrids_m,ngrids_y,ngrids_z,ngrids_o
INTEGER :: gleng1_min, gleng1_max

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

write(*,*)'NCPUS',ncpu
DO i=1,ncv
  grid(1:3,i)=grid0(1:3,i)
END DO
rank=0

IF(ncv.EQ.1)THEN
 ngrids_y=1
 ngrids_z=1
ELSE IF(ncv.EQ.2)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=1
ELSE IF(ncv.EQ.3)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=(nint((grid0(2,3)-grid0(1,3))/grid0(3,3))+1)
END IF

if(ncpu.eq.1) then 
gleng1_min=1
gleng1_max=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
end if

!Distribute X grids
if(icpu.eq.0)WRITE(*,'(3A12,3A16)') 'CPU','CV', 'GRID SIZE', 'GRID MIN', 'GRID MAX', 'GRID BIN'
CALL Sync_procs
IF(ncpu.GT.1)THEN
  ngrids=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  write(*,*)'NEW_GRID',ngrids,icpu,ncpu
  ngrids_o=ngrids
  ngrids=ngrids/ncpu
  IF(icpu.eq.ncpu-1)THEN
    ngrids_m=ngrids+mod(ngrids_o,ncpu)
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids_m-1)*grid0(3,1)
  ELSE
    ngrids_m=ngrids
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids-1)*grid0(3,1)
  END IF
  CALL Sync_procs
  WRITE(*,'(3I12,3F16.6)') icpu, 1, ngrids_m, grid(1,1), grid(2,1), grid(3,1)
  rank=ngrids_z*ngrids_y*ngrids*icpu
  gleng1_min=ngrids*icpu+1
  gleng1_max=ngrids*(icpu+1)
  if(icpu.eq.ncpu-1) gleng1_max=ngrids*(icpu+1)+mod(ngrids_o,ncpu)
END IF 
END 
!****************************************************************************************!



!****************************************************************************************!

SUBROUTINE MPI_Start()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
  call MPI_INIT(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE MPI_Stop()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
call MPI_FINALIZE(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE get_ncpu(ncpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: ncpu, i_err
!ncpu=1
!#if defined (_PARALLEL)
call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE get_cpuid(icpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: icpu, i_err
!icpu=0
!#if defined (_PARALLEL)
call MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE IBcast(myint,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: leng, myint(*), i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE RBcast(myreal,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: myreal(*)
INTEGER :: leng, i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE Sync_Procs
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER i_err
!#if defined (_PARALLEL)
call MPI_Barrier(MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE Set_Parent(parent)
IMPLICIT NONE
INCLUDE 'mpif.h'
LOGICAL :: parent
INTEGER :: icpu, i_err
parent=.false.
!icpu=0
!#if defined (_PARALLEL)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
IF(icpu.eq.0)parent=.true.
END
!****************************************************************************************!

SUBROUTINE GlobSumR(myreal_in,myreal_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
REAL*8 :: myreal_in(*), myreal_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE GlobSumI(myint_in,myint_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
INTEGER :: myint_in(*), myint_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
call MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!
SUBROUTINE get_filename(head,filename,ir)
IMPLICIT NONE
CHARACTER (LEN=*) :: head
CHARACTER (LEN=*) :: filename
INTEGER :: ir

 IF(ir.LT.10)THEN
   WRITE(FILENAME,'(A,I1)')TRIM(HEAD),ir
 ELSE IF(ir.LT.100)THEN
   WRITE(FILENAME,'(A,I2)')TRIM(HEAD),ir
 ELSE
   PRINT *,'error! GET_FILANAME ERROR: ir=',ir
   STOP
 END IF

END
!****************************************************************************************!
