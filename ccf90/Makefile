FC := gfortran

FFLAGS := -O3 -g -pedantic -Wall -Wextra  -std=f2018 -fdefault-real-8 -fdefault-integer-8
LDFLAGS := -llapack -lblas
LIBS := -Wl, --start-group /opt/intel/lib/intel64_mac  

EXEFILE = run

.PHONY: all clean

all: $(EXEFILE)
       
$(EXEFILE): dgemm_module.o sort_module.o diis.o einsum_module.o permutils.o system_types.o integral_types.o integrals.o printing.o cc_energy.o hbar.o testing_module.o mp2.o ccsd_module.o main.o 
		 $(FC) $^ -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.mod: ;

clean:
	rm -f *.o *.mod $(EXEFILE)
