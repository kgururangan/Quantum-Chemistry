program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use mp2, only: calculate_mp2
        use ccsd_module, only: ccsd 
        use printing, only: print_calc_params, print_header
        use testing_module, only: test_einsum, write_hbar_ccs, test_update_t1b
        use diis, only: solvegauss

        implicit none

        integer, parameter :: Norb = 14, Nfroz=0, Nelec=10
        type(e1int_t) :: fA, fB
        type(e2int_t) :: vA, vB, vC
        type(sys_t) :: sys
        integer, parameter :: ndiis = 6, maxit = 1000
        real, parameter :: shift = 0.0, tol = 1.0e-08
        real, allocatable :: t1a(:,:), t1b(:,:), t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:)
        real :: A(10,10), b(10), x(10)
        integer :: i, j

        call print_header()
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/CCpy_tests/H2O-2Re-DZ/onebody.inp',&
                           '/users/karthik/Desktop/CC_matlab_tests/CCpy_tests/H2O-2Re-DZ/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        !call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/onebody.inp',&
        !                   '/users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/twobody.inp',&
        !                   Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        call print_calc_params(sys)


        allocate(t1a(sys%Nunocc_a,sys%Nocc_a),t1b(sys%Nunocc_b,sys%Nocc_b),&
                 t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                 t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                 t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b))

        call ccsd(sys,fA,fB,vA,vB,vC,ndiis,maxit,shift,tol,t1a,t1b,t2a,t2b,t2c)

       ! do i = 1,10
       !   do j = 1,10
       !      A(i,j) = float(i) + float(j)
       !   end do
       !   b(i) = float(i)
       !end do

       !call solvegauss(A,b,10,x)

       !print*,x






end program main

