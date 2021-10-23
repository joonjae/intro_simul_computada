MAKEFILE = Makefile
exe = problema_ising
#exe = ising_2
fcomp = gfortran #ifort # /opt/intel/compiler70/ia32/bin/ifc  
# Warning: the debugger doesn't get along with the optimization options
# So: not use -O3 WITH -g option
#flags =  -O3 
flags =  -O3 
# Remote compilation
OBJS = ziggurat.o mymodule.o problema_ising.o
#OBJS = ziggurat.o mymodule.o ising_2.o

.SUFFIXES:            # this deletes the default suffixes 
.SUFFIXES: .f90 .o    # this defines the extensions I want 

.f90.o:  
	$(fcomp) -c $(flags) $<  
        

$(exe):  $(OBJS) Makefile 
	$(fcomp) $(flags) -o $(exe) $(OBJS)


clean:
	rm ./*.o ./*.mod	

problema2b.o: problema_ising.f90 ziggurat.o mymodule.o
#problema2b.o: ising_2.f90 ziggurat.o mymodule.o
