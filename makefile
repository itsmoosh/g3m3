# The path to SHTOOLS in $(INCLUDE) and $(LIBS) may have to be changed
# if you install SHTOOLS in a different location.
# 
# The -flto option in $(FCFLAGS) will have to be removed if using a
# version of gfortran prior to 4.6.

FC = ncargf90
FCFLAGS = -zero -save -openmp 
FCFLAGS += -liomp5 
FILES = saturn3d_main.f90 \
	    saturns_main.f90 \
		saturng_main.f90 \
		saturnd_main.f90
LIBS = -L/opt/intel/composer_xe_2013/lib/intel64

saturn.x: $(OBJECTS)
	$(FC) $(FCFLAGS) -o saturn.x $(FILES)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.x
	rm -f fluid*
	rm -f gmeta
	rm -f *.dat
