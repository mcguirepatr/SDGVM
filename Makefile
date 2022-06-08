FFLAGS= -fbounds-check -O3 -Wuninitialized -ftrapv -fimplicit-none -fno-automatic -fbacktrace 
NETCDFFLAGS= `nf-config --fflags`
NETCDFLIBS= `nf-config --flibs`
FF = gfortran

OBJ = states_convertSDGVM_func_v7a.o sdgvm0.o sdgvm1.o data.o growth.o parameter_adjustment.o hydrology.o phenology.o func.o doly.o soil.o nppcalc.o light.o sunshade.o weathergenerator.o metdos.o luna.o

sdgvm0:	$(OBJ)
	$(FF) $(FFLAGS)  -o sdgvm0 $(OBJ) $(NETCDFFLAGS) $(NETCDFLIBS)

states_convertSDGVM_func_v7a.o : states_convertSDGVM_func_v7a.f90
	$(FF)  $(FFLAGS)  -c $< $(NETCDFFLAGS) $(NETCDFLIBS)

%.o: %.f
	$(FF)  $(FFLAGS)  -c $<



clean:
	rm -f *.o *.il *.OBJ *~ *.exe a.out *.MOD

