FC := gfortran

FFLAGS := -O0 -g
LDFLAGS := 

EXEFILE = run

.PHONY: all clean

#VPATH := utils/ : derived_types/ : mmd_integrals 

all: $(EXEFILE)
       
$(EXEFILE): main.o ao_integral_module.o mmd_integrals.o utils.o mmd_primitives.o orbital_types.o gaussian_integral_utils.o special_functions.o constants.o
		 $(FC) $^ -o $@ $(LDFLAGS)

main.o: ao_integral_module.o utils.o mmd_primitives.o orbital_types.o gaussian_integral_utils.o special_functions.o constants.o

ao_integral_module.o: mmd_primitives.o

mmd_integrals.o: mmd_primitives.o

mmd_primitives.o: gaussian_integral_utils.o constants.o

gaussian_integral_utils.o: special_functions.o constants.o

special_functions.o: constants.o

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.mod: ;

clean:
	rm -f *.o *.mod $(EXEFILE)
