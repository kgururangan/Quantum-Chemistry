program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use mp2, only: calculate_mp2
        use ccsd_module, only: ccsd 
        use printing, only: print_calc_params, print_header
        use cis, only: cis_mat
        use testing_module, only: test_eomccsd

        implicit none

        integer, parameter :: Norb = 14, Nfroz=0, Nelec=10
        type(e1int_t) :: fA, fB, H1A, H1B
        type(e2int_t) :: vA, vB, vC, H2A, H2B, H2C
        type(sys_t) :: sys
        real :: Emp2, Ecorr
        integer, parameter :: ndiis = 6, maxit = 100
        real, parameter :: shift = 0.0, tol = 1.0e-08
        real, allocatable :: t1a(:,:), t1b(:,:), t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:), &
                             c1a(:,:), c1b(:,:), w_cis(:)


        call print_header()
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/CCpy_tests/H2O-2Re-DZ/onebody.inp',&
                           '/users/karthik/Desktop/CC_matlab_tests/CCpy_tests/H2O-2Re-DZ/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        call print_calc_params(sys)

        !allocate(c1a(sys%Nunocc_a*sys%Nocc_a,5),c1b(sys%Nunocc_b*sys%Nocc_b,5),w_cis(5))
        !call cis_mat(sys,fA,fB,vA,vB,vC,5,c1a,c1b,w_cis)
        !print*,w_cis
        !deallocate(c1a,c1b,w_cis)


        !call calculate_mp2(sys,fA,fB,vA,vB,vC,Emp2)


        allocate(t1a(sys%Nunocc_a,sys%Nocc_a),t1b(sys%Nunocc_b,sys%Nocc_b),&
                 t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                 t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                 t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b))

        call test_eomccsd(sys,fA,fB,vA,vB,vC)

        deallocate(t1a,t1b,t2a,t2b,t2c)
end program main

