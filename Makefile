#FFLAGS= -fbounds-check -O3 -Wuninitialized -ftrapv -fimplicit-none -fno-automatic -fbacktrace 
#For debugging:
FFLAGS= -g -fbounds-check -O0 -Wuninitialized -ftrapv -fimplicit-none -fno-automatic -fbacktrace -fcheck=all
FF = gfortran

OBJ = sdgvm0.o sdgvm1.o data.o growth.o parameter_adjustment.o hydrology.o phenology.o func.o doly.o soil.o nppcalc.o light.o sunshade.o weathergenerator.o metdos.o luna.o

sdgvm0:	$(OBJ)
	$(FF) $(FFLAGS)  -o sdgvm0 $(OBJ)

%.o: %.f
	$(FF)  $(FFLAGS)  -c $<



clean:
	rm -f *.o *.il *.OBJ *~ *.exe a.out *.MOD

