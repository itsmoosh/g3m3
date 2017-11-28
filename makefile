# The path to SHTOOLS in $(INCLUDE) and $(LIBS) may have to be changed
# if you install SHTOOLS in a different location.
# To debug, run make debug and then use idb debug.x

FC = ncargf90
DEBUGGER = ifort
FCFLAGS = -zero -save -openmp -liomp5
DBFLAGS = -g -traceback
DBLIBS += -L/usr/local/pkgs/ncarg/lib -L/usr/local/lib -lncarg -lncarg_gks -lncarg_c -lXpm -lX11 -lXext -lcairo -lfontconfig -lpixman-1 -lfreetype -lexpat -lpng -lz -lpthread -lbz2 -lXrender
FILES = code/m_*.f90 code/3d_*.f90 code/s_*.f90 code/g_*.f90 code/c_*.f90 code/main.f90
#	LIBS is usually added to $PATH
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
	rm -f fluid*
	rm nohup.out

test:
	$(DEBUGGER) $(FCFLAGS) $(DBFLAGS) -o test.x $(TESTF) $(TESTP) $(DBLIBS)

debug:
	$(DEBUGGER) $(FCFLAGS) $(DBFLAGS) -o debug.x $(FILES) $(DBLIBS) 
