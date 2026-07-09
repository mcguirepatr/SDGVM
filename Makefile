OBJ = sdgvm0.o sdgvm1.o data.o growth.o parameter_adjustment.o hydrology.o phenology.o func.o doly.o soil.o nppcalc.o light.o sunshade.o weathergenerator.o metdos.o luna.o
FFLAGS= -fbounds-check -O3 -Wuninitialized -ftrapv -fimplicit-none -fno-automatic -fbacktrace 
#For debugging:
#FFLAGS= -g -fbounds-check -O0 -Wuninitialized -ftrapv -fimplicit-none -fno-automatic -fbacktrace -fcheck=all
FF = gfortran

# allow to work with macs
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SDKROOT := $(shell xcrun --show-sdk-path 2>/dev/null)
  ifneq ($(SDKROOT),)
    LDFLAGS += -Wl,-syslibroot,$(SDKROOT)
  endif
endif

# compile and link
sdgvm0:	$(OBJ)
	$(FF) $(FFLAGS) $(LDFLAGS)  -o sdgvm0 $(OBJ)

%.o: %.f
	$(FF) $(FFLAGS) $(LDFLAGS)  -c $<


clean:
	rm -f *.o *.il *.OBJ *~ *.exe a.out *.MOD

