MODULE error_msg
CONTAINS
!-------------------------------------------------------------------------
SUBROUTINE cv_error (j,grid)
IMPLICIT NONE
INTEGER :: j
CHARACTER*3  :: grid

WRITE(6,'(A)')
WRITE(6,'(A,I2,A,A,A,2X,A)')"ERROR:: CV ",j," grid(",grid,") OUT of RANGE"
WRITE(6,'(A)')
STOP

END SUBROUTINE cv_error
!-------------------------------------------------------------------------
SUBROUTINE step_error (t_min,t_max)
IMPLICIT NONE
INTEGER :: t_min,t_max

IF (t_min .ge. t_max) THEN
WRITE(6,'(A)')
WRITE(6,'(A)')"ERROR:: tmax CAN NOT BE LESS THAN tmin"
WRITE(6,'(A)')
STOP
ENDIF

END SUBROUTINE step_error
!-------------------------------------------------------------------------
SUBROUTINE prog_error

WRITE(6,'(A)')
WRITE(6,'(A)')"ERROR :: THIS REWEIGHTING CODE ONLY WORKS WITH CPMD or PLUMED" 
WRITE(6,'(A)')"PLEASE REFER TO THE github PAGE : https://github.com/rahulumrao/Reweighing-TASS-1.1/" 
WRITE(6,'(A)')
STOP
END SUBROUTINE prog_error
!-------------------------------------------------------------------------
SUBROUTINE eof_error

WRITE(6,'(A)')
WRITE(6,'(A)')"ERROR :: EOF is RAISED, THE input() REACHES THE END OF YOUR CODE WITHOUT READING ANY DATA "
WRITE(6,'(A)')
STOP
END SUBROUTINE eof_error
!-------------------------------------------------------------------------
SUBROUTINE grid_error
WRITE(6,'(A)')
WRITE(6,'(A)')"ERROR :: NCV is NOT EQUAL TO GRIDS IN INPUT FILE"
WRITE(6,'(A)')
STOP
END SUBROUTINE grid_error
!-------------------------------------------------------------------------
SUBROUTINE max_t(max_step,t_min,t_max,vt_max,md_steps)
IMPLICIT NONE
INTEGER :: i,nr,t_max,t_min,vt_max,md_steps
LOGICAL :: max_step

IF (max_step) THEN
  t_max = vt_max
  IF(md_steps .le. vt_max) THEN
    t_max = md_steps
  ENDIF
ELSE
  t_max = md_steps
ENDIF
END SUBROUTINE max_t

END MODULE        
