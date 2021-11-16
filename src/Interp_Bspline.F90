!*****************************************************************************************
    PROGRAM interp_bspline
    USE GetSteps
    USE bspline_module
    USE bspline_kinds_module, only: wp, ip

    IMPLICIT NONE

    INTEGER(ip)              :: nx           !! number of points in x
    INTEGER(ip)              :: nxv          !! number of points to evaluate interpolant
    INTEGER(ip),PARAMETER    :: kx    = 4    !! order in x
    INTEGER(ip),PARAMETER    :: iknot = 0    !! automatically select the knots

    INTEGER                  :: n,ios,spoints
    INTEGER(ip)              :: i,j,idx,iflag,inbvx,iloy
    REAL(wp),ALLOCATABLE     :: x(:),f1(:)
    REAL(wp),ALLOCATABLE     :: xval(:)
    REAL(wp),ALLOCATABLE     :: tx(:)
    REAL(wp)                 :: val,tru,err,errmax,nb
    REAL(wp),DIMENSION(3*kx) :: w1_1d !! work array
    LOGICAL                  :: extrap
    CHARACTER*50             :: filename,outfile

    idx = 0
    WRITE(6,'(A)') "ENTER THE DATA FILENAME :"
    READ(5,*) filename
    OPEN(20,FILE=TRIM(filename),STATUS='old',IOSTAT=ios)
    IF (ios .ne. 0) STOP '"ERROR : Input File does not exist"'
    WRITE(6,'(A)') "ENTER THE OUTPUT FILENAME :"
    READ(5,*) outfile
    WRITE(6,'(A)') "ORDER OF INTERPOLATION (number of points to evaluate interpolant)"
    READ(5,*)spoints
    
    CALL get_steps(20,n)
    nx = n ; nxv = spoints*n

    ALLOCATE(x(nx))
    ALLOCATE(xval(nxv))
    ALLOCATE(tx(nx+kx))
    ALLOCATE(f1(nx))
    
    DO i = 1 ,nx
      READ(20,*)x(i),f1(i)
    !  WRITE(*,*)x(i),f1(i)
    ENDDO
    CLOSE(20)

    nb = (x(nx) - x(1))/nxv
    DO i = 1,nxv
      xval(i) = x(1)+nb*i
    ENDDO

    WRITE(6,'(I6,2X,A,A)')nx,"INITIAL DATA POINTS IN FILE :",filename
    101 FORMAT (A,I6,1X,A,1X,F4.2,2X,A,2X,F4.2,A)
    WRITE(*,101)"INTERPOLATION OF",nxv,"POINTS WLL BE DONE BETWEEN : [ ",x(1),'<-->',x(nx),' ]'
    !have to set these before the first evaluate call:
    inbvx = 1 ; iloy  = 1 ; iflag = 0
    ! initialize
    CALL db1ink(x,nx,f1,kx,iknot,tx,f1,iflag)

    IF (iflag/=0) THEN
        WRITE(*,*) 'Error initializing 1D spline: '//get_status_message(iflag)
    END IF
        
       OPEN(21,FILE=outfile)
        errmax = 0.0_wp
        DO i=1,nxv
           CALL db1val(xval(i),idx,tx,nx,kx,f1,val,iflag,inbvx,w1_1d,extrap=extrap)
           WRITE(21,*)xval(i), val
        END DO

        WRITE(*,*) ''
        WRITE(*,'(A,A)') '1D: INTERPOLATED DATA WRITTEN IN FILE : ', outfile
        WRITE(*,*) ''

DEALLOCATE(x)
DEALLOCATE(xval)
DEALLOCATE(f1)
CLOSE(21)
END
