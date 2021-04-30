program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use mp2, only: calculate_mp2
        use ccsd_module, only: hbar_ccs_intermediates 
        use printing, only: print_calc_params, print_header
        use testing_module, only: test_einsum, write_hbar_ccs, test_update_t1b

        implicit none

        integer, parameter :: Norb = 14, Nfroz=0, Nelec=10
        type(e1int_t) :: fA, fB
        type(e2int_t) :: vA, vB, vC
        type(sys_t) :: sys

        call print_header()
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/CCpy_tests/H2O-2Re-DZ/onebody.inp',&
                           '/users/karthik/Desktop/CC_matlab_tests/CCpy_tests/H2O-2Re-DZ/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        !call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/onebody.inp',&
        !                   '/users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/twobody.inp',&
        !                   Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        call print_calc_params(sys)

        ! Testing einsum
        !call test_einsum(sys,fA,fB,vA,vB,vC)

        ! MP2 energy correction
        !call calculate_mp2(sys,fA,fB,vA,vB,vC,Emp2)
        !print*,'E(MP2) = ',Emp2

        ! Checking HBar CCS construction
        !call write_hbar_ccs(sys,fA,fB,vA,vB,vC)

        call test_update_t1b(sys,fA,fB,vA,vB,vC)

end program main

