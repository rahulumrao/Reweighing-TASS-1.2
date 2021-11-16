#makefile
test.x: bspline_kinds_module.o  bspline_module.o  bspline_oo_module.o  bspline_sub_module.o
	gfortran -o test.x test.F90 bspline_kinds_module.o  bspline_module.o  bspline_oo_module.o  bspline_sub_module.o

clean   :
	rm *.x

