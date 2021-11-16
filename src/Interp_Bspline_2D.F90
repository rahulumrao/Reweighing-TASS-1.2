!*****************************************************************************************
    PROGRAM interp_bspline_2D
    USE GetSteps
    USE bspline_module
    USE bspline_kinds_module, only: wp, ip

    IMPLICIT NONE

    INTEGER(ip)              :: nx,ny            !! number of points in x
    INTEGER(ip)              :: nxv,nyv          !! number of points to evaluate interpolant
    INTEGER(ip),PARAMETER    :: kx    = 4    !! order in x
    INTEGER(ip),PARAMETER    :: ky    = 4    !! order in y
    INTEGER(ip),PARAMETER    :: iknot = 0    !! automatically select the knots

    INTEGER                  :: n,m,ios,spoints
    INTEGER(ip)              :: i,j,k,idx,idy,iflag,inbvx,inbvy,iloy,iloz
    REAL(wp),ALLOCATABLE     :: x(:),y(:),f2(:,:)
    REAL(wp),ALLOCATABLE     :: xval(:),yval(:)
    REAL(wp),ALLOCATABLE     :: tx(:),ty(:)
    REAL(wp)                 :: val(2),tru,err,errmax,nb,nc
    REAL(wp),DIMENSION(ky)   :: w1_2d
    REAL(wp),DIMENSION(3*max(kx,ky)) :: w2_2d !! work array
    LOGICAL                  :: extrap
    CHARACTER*50             :: filename,outfile

    idx = 0
    WRITE(6,'(A)') "ENTER THE DATA FILENAME :"
    READ(5,*) filename
    OPEN(20,FILE=TRIM(filename),STATUS='old',IOSTAT=ios)
    IF (ios .ne. 0) STOP '"ERROR : Input File does not exist"'
    WRITE(6,'(A)') "ENTER THE OUTPUT FILENAME :"
    READ(5,*) outfile
    WRITE(6,'(A)') "BIN SIZE IN x-axis:"
    READ(5,*) nx
    WRITE(6,'(A)') "ORDER OF INTERPOLATION (number of points to evaluate interpolant)"
    READ(5,*)spoints
    
    CALL get_steps(20,n)
    ny = (n-nx)/nx ; nxv = spoints*nx ; nyv = spoints*ny
    PRINT*,nx 
    ALLOCATE(x(nx),y(ny))
    ALLOCATE(xval(nxv),yval(nyv))
    ALLOCATE(tx(nx+kx),ty(ny+ky))
    ALLOCATE(f2(nx,ny))
    
    DO i = 1 ,nx
      DO j = 1,ny
        READ(20,*)x(i),y(j),f2(i,j)
!        WRITE(6,*)x(i),y(j),f2(i,j)
      ENDDO
      READ(20,*)
    ENDDO
    CLOSE(20)

    nb = (x(nx) - x(1))/nxv
    DO i = 1,nxv
      xval(i) = x(1)+nb*i
    ENDDO

    nc = (y(ny) - y(1))/nyv
    DO i = 1,nyv
      yval(i) = y(1)+nc*i
    ENDDO

    WRITE(6,'(I6,2X,A,A)')nx,"INITIAL DATA POINTS IN FILE :",filename
    101 FORMAT (A,I6,1X,A,1X,F4.2,2X,A,2X,F4.2,A)
    WRITE(*,101)"INTERPOLATION OF",nxv,"POINTS WLL BE DONE BETWEEN : [ ",x(1),'<-->',x(nx),' ]'
    !have to set these before the first evaluate call:
    inbvx = 1 ; iloy  = 1 ; iflag = 0
    inbvy = 1 ; iloz  = 1
    ! initialize
    CALL db2ink(x,nx,y,ny,f2,kx,ky,iknot,tx,ty,f2,iflag)

    IF (iflag/=0) THEN
        WRITE(*,*) 'Error initializing 1D spline: '//get_status_message(iflag)
    END IF
        
       OPEN(21,FILE=outfile)
        errmax = 0.0_wp
        DO i=1,nxv
          DO j=1,nyv
           CALL db2val(xval(i),yval(j),idx,idy,tx,ty,nx,ny,kx,ky,f2,val(2),iflag,inbvx,inbvy,iloy,w1_2d,w2_2d,extrap=extrap)
           WRITE(21,*)xval(i),yval(j),val(2)
          ENDDO
          WRITE(21,*)
        END DO

        WRITE(*,*) ''
        WRITE(*,'(A,A)') '1D: INTERPOLATED DATA WRITTEN IN FILE : ', outfile
        WRITE(*,*) ''

DEALLOCATE(x,y)
DEALLOCATE(xval,yval)
DEALLOCATE(f2)
CLOSE(21)
END
