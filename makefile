# The path to SHTOOLS in $(INCLUDE) and $(LIBS) may have to be changed
# if you install SHTOOLS in a different location.

FC = ncargf90
FCFLAGS = -zero -save -openmp
FCFLAGS += -liomp5 -g -traceback 
FILES = code/*.f90
LIBS = -L/opt/intel/composer_xe_2013/lib/intel64
TESTP = code/test.f90
#	For testing, edit the following line to refer to the file to test.
TESTF = code/test.f90

multifluid: $(OBJECTS)
	$(FC) $(FCFLAGS) -o master.x $(FILES)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.x
	rm -f *.log
	rm -f gmeta
	rm -f *.dat

test:
	$(FC) $(FCFLAGS) -o test.x $(TESTP) $(TESTF)
