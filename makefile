# The path to SHTOOLS in $(INCLUDE) and $(LIBS) may have to be changed
# if you install SHTOOLS in a different location.
# 
# The -flto option in $(FCFLAGS) will have to be removed if using a
# version of gfortran prior to 4.6.

FC = ncargf90
FCFLAGS = -zero -save -openmp -debug -O0 -g -traceback 
FCFLAGS += -liomp5 
FILES = earth3d_main.f90 \
		earthd_main.f90 \
		earthg_main.f90 \
		earths_main.f90
LIBS = -L/opt/intel/composer_xe_2013/lib/intel64

earth.x: $(OBJECTS)
	$(FC) $(FCFLAGS) -o earth.x $(FILES)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.x
	rm -f *.log
	rm -f gmeta
	rm -f *.dat
