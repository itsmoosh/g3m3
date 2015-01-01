# The path to SHTOOLS in $(INCLUDE) and $(LIBS) may have to be changed
# if you install SHTOOLS in a different location.
# 
# The -flto option in $(FCFLAGS) will have to be removed if using a
# version of gfortran prior to 4.6.

FC = ncargf90
FCFLAGS = -zero -save -openmp 
FCFLAGS += -liomp5 
FILES = saturn3d_main.f saturns_main.f \
		saturng_main.f saturnd_main.f
#OBJECTS = sat_titans_mvgd.o sat_titang_mvgd.o \
		  sat_titand_mvgd.o sat_titan3d_mvgd.o
#INCLUDE = -I/usr/local/ncarg/include \
		  -I/opt/intel/fce/10.0.023/include
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
