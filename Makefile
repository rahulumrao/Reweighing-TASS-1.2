#Makefile
#---------------------------------------------
# Written by - Rahul Verma (vrahul@iitk.ac.in)
#---------------------------------------------

#.RECIPEPREFIX+=

default:
	@echo "Type..."
	@echo "make install   : create executabls"
	@echo "make bspline   : compile B-spline modules"
	@echo "make clean     : remove object and mod files"
	@echo "make distclean : clean the directory"

TOPDIR=/home/vrahul/My_Program/Rahul/Probability/Reweighing-TASS_1.2
LIBDIR=$(TOPDIR)/lib
SRCDIR=$(TOPDIR)/src
BINDIR=$(TOPDIR)/bin
BSPDIR=$(TOPDIR)/bspline-fortran
PYPDIR=$(BSPDIR)/src/tests/pyplot-fortran
BSPFLAGS=$(TOPDIR)/bspline-fortran/build/libbspline-fortran.a

F90=gfortran
#FC=mpif90
#FCFLAGS=-g3 -fcheck=all -fbacktrace

TARGET=tass.x
EXE=1d_bspline.x
EXE2=2d_bspline.x


SPOBJECT=bspline_kinds_module.o bspline_sub_module.o bspline_oo_module.o bspline_module.o pyplot_module.o
OBJECTS=Ansi_Colors.o GetSteps.o GetFileName.o Input_file.o Error_msg.o MTD_Unbais.o MTD_Potential.o US_Prob.o \
	 US_MTD.o US_TEMP.o Mean_Force.o Stat_Error_Estimate.o \
	 $(SPOBJECT) B_Spline.o

1d_bspline.x		: $(SPOBJECT) GetSteps.o		   ; $(F90) -o $(EXE) $(SRCDIR)/Interp_Bspline.F90 GetSteps.o $(SPOBJECT)
2d_bspline.x		: $(SPOBJECT) GetSteps.o		   ; $(F90) -o $(EXE2) $(SRCDIR)/Interp_Bspline_2D.F90 GetSteps.o $(SPOBJECT)
tass.x                    : $(OBJECTS) 				   ; $(F90) -o $(TARGET) $(FCFLAGS) $(SRCDIR)/Main.F90 $(OBJECTS)

Ansi_Colors.o		:   $(SRCDIR)/Ansi_Colors.F90              ; $(F90) -c $(SRCDIR)/Ansi_Colors.F90
GetSteps.o  		:   $(SRCDIR)/GetSteps.F90 	           ; $(F90) -c $(SRCDIR)/GetSteps.F90
GetFileName.o 		:   $(SRCDIR)/GetFileName.F90 	           ; $(F90) -c $(SRCDIR)/GetFileName.F90
Input_file.o 		:   $(SRCDIR)/Input_file.F90 	           ; $(F90) -c $(SRCDIR)/Input_file.F90
Error_msg.o               :   $(SRCDIR)/Error_msg.F90                ; $(F90) -c $(SRCDIR)/Error_msg.F90
MTD_Unbais.o    	        :   $(SRCDIR)/MTD_Unbais.F90	           ; $(F90) $(FCFLAGS) -c $(SRCDIR)/MTD_Unbais.F90
MTD_Potential.o 	        :   $(SRCDIR)/MTD_Potential.F90	           ; $(F90) $(FCFLAGS) -c $(SRCDIR)/MTD_Potential.F90
US_Prob.o   		:   $(SRCDIR)/US_Prob.F90	           ; $(F90) $(FCFLAGS) -c $(SRCDIR)/US_Prob.F90
US_MTD.o    		:   $(SRCDIR)/US_MTD.F90	           ; $(F90) -c $(SRCDIR)/US_MTD.F90
US_TEMP.o   		:   $(SRCDIR)/US_TEMP.F90	           ; $(F90) $(FCFLAGS) -c $(SRCDIR)/US_TEMP.F90
Mean_Force.o 		:   $(SRCDIR)/Mean_Force.F90	           ; $(F90) $(FCFLAGS) -c $(SRCDIR)/Mean_Force.F90
Stat_Error_Estimate.o	:   $(SRCDIR)/Stat_Error_Estimate.F90	   ; $(F90) $(FCFLAGS) -c $(SRCDIR)/Stat_Error_Estimate.F90
bspline_kinds_module.o    :   $(BSPDIR)/src/bspline_kinds_module.f90 ; $(F90) -c $(BSPDIR)/src/bspline_kinds_module.f90
bspline_sub_module.o      :   $(BSPDIR)/src/bspline_sub_module.f90   ; $(F90) -c $(BSPDIR)/src/bspline_sub_module.f90
bspline_oo_module.o	:   $(BSPDIR)/src/bspline_oo_module.f90    ; $(F90) -c $(BSPDIR)/src/bspline_oo_module.f90
bspline_module.o	        :   $(BSPDIR)/src/bspline_module.f90	   ; $(F90) -c $(BSPDIR)/src/bspline_module.f90
pyplot_module.o		:   $(PYPDIR)/src/pyplot_module.f90        ; $(F90) -c $(PYPDIR)/src/pyplot_module.f90
B_Spline.o		:   $(SRCDIR)/B_Spline.F90                 ; $(F90) $(FCFLAGS) -c $(SRCDIR)/B_Spline.F90

.PHONY: install clean distclean bspline

install	: $(TARGET)
	@mkdir -p $(TOPDIR)/bin
	@mkdir -p $(TOPDIR)/lib
	@mv tass.x $(TOPDIR)/bin
	@mv *.o *.mod $(TOPDIR)/lib

bspline : $(EXE) $(EXE2)
	@mkdir -p $(TOPDIR)/bin
	@mkdir -p $(TOPDIR)/lib
	@mv 1d_bspline.x $(TOPDIR)/bin
	@mv 2d_bspline.x $(TOPDIR)/bin
	@mv *.o *.mod $(TOPDIR)/lib

clean   :
	@rm -rf $(TOPDIR)/bin/Probability_analysis.x
	@rm -rf $(TOPDIR)/lib/*.o
	@rm -rf $(TOPDIR)/lib/*.mod
	@rm *.o *.mod

distclean :
	@echo "Removing executables and libraries from $(TOPDIR)"
	@echo "rm -rf bin/* "
	@echo "rm -rf lib/* "
	@rm -rf $(TOPDIR)/bin
	@rm -rf $(TOPDIR)/lib
