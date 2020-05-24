FC := gfortran

FFLAGS := -O0 -g
LDFLAGS := 

EXEFILE = run

.PHONY: all clean

all: $(EXEFILE)
       
$(EXEFILE): main.o ao_integrals.o gaussian_integration.o utils_fcn.o
		 $(FC) $^ -o $@ $(LDFLAGS)

main.o: ao_integrals.o gaussian_integration.o utils_fcn.o
ao_integrals.o: gaussian_integration.o utils_fcn.o
gaussian_integration.o: utils_fcn.o

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.mod: ;

clean:
	rm -f *.o *.mod $(EXEFILE)
