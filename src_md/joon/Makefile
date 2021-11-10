MAKEFILE = Makefile
exe = problema_dm
fcomp = gfortran #ifort # /opt/intel/compiler70/ia32/bin/ifc  
flags =  -O3 
# Remote compilation
OBJS = ziggurat.o globals.o mymodule.o problema_dm.o
#OBJS = ziggurat.o mymodule.o ising_2.o

.SUFFIXES:            # this deletes the default suffixes 
.SUFFIXES: .f90 .o    # this defines the extensions I want 

.f90.o:  
	$(fcomp) -c $(flags) $<  
        

$(exe):  $(OBJS) Makefile 
	$(fcomp) $(flags) -o $(exe) $(OBJS)


clean:
	rm ./*.o ./*.mod	
