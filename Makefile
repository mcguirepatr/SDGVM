FFLAGS= -fbounds-check -O3 -Wuninitialized -ftrapv -fimplicit-none -fno-automatic -fbacktrace 
#For debugging:
#FFLAGS= -g -fbounds-check -O0 -Wuninitialized -ftrapv -fimplicit-none -fno-automatic -fbacktrace -fcheck=all
NETCDFFLAGS= `nf-config --fflags`
NETCDFLIBS= `nf-config --flibs`
#HDF5LIBS= `h5cc -show`
HDF5LIBS= -lhdf5 
FF = gfortran

OBJ = states_convertSDGVM_func.o sdgvm0.o sdgvm1.o data.o growth.o parameter_adjustment.o hydrology.o phenology.o func.o doly.o soil.o nppcalc.o light.o sunshade.o weathergenerator.o metdos.o luna.o
#OBJ = tmp.o

sdgvm0:	$(OBJ)
	$(FF) $(FFLAGS)  -o sdgvm0 $(OBJ) $(NETCDFFLAGS) $(NETCDFLIBS) $(HDF5LIBS)

tmp:	$(OBJ)
	$(FF) $(FFLAGS)  -o tmp $(OBJ) $(NETCDFFLAGS) $(NETCDFLIBS) $(HDF5LIBS)

tmp.o : tmp.f90
	$(FF)  $(FFLAGS)  -c $< $(NETCDFFLAGS) $(NETCDFLIBS) $(HDF5LIBS)

states_convertSDGVM_func.o : states_convertSDGVM_func.f90
	$(FF)  $(FFLAGS)  -c $< $(NETCDFFLAGS) $(NETCDFLIBS) $(HDF5LIBS)

%.o: %.f
	$(FF)  $(FFLAGS)  -c $<



clean:
	rm -f *.o *.il *.OBJ *~ *.exe a.out *.MOD

