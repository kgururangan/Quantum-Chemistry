FC := gfortran

FFLAGS := -O3 -g -pedantic -Wall -Wextra  -std=f2018 -fdefault-real-8 -fdefault-integer-8
LDFLAGS := -llapack -lblas
LIBS := -Wl, --start-group /opt/intel/lib/intel64_mac  

EXEFILE = run

.PHONY: all clean

all: $(EXEFILE)
       
$(EXEFILE): system_types.o integral_types.o printing.o dgemm_module.o dgeev_module.o sort_module.o permutils.o einsum_module.o integrals.o cis_updates.o davidson_module.o main.o 
		 $(FC) $^ -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.mod: ;

clean:
	rm -f *.o *.mod $(EXEFILE)
