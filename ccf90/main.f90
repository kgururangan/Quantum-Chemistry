program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use mp2, only: calculate_mp2
        use ccsd_module, only: ccsd 
        use printing, only: print_calc_params, print_header
        use cis, only: cis_mat
        use eomccsd, only: davidson_eomccsd
        use hbar, only: hbar_ccsd
        use leftccsd, only: left_ccsd

        implicit none

        integer, parameter :: Norb = 26, Nfroz=0, Nelec=6
        type(e1int_t) :: fA, fB, H1A, H1B
        type(e2int_t) :: vA, vB, vC, H2A, H2B, H2C
        type(sys_t) :: sys
        real :: Emp2, Ecorr
        integer, parameter :: ndiis = 6, maxit = 200, nroot = 15, init_guess = 1
        real, parameter :: shift = 1.0, tol = 1.0e-08
        real, allocatable :: t1a(:,:), t1b(:,:), t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:), &
                             c1a(:,:), c1b(:,:), w_cis(:), Rvec(:,:), omega(:), Lvec(:,:)
        integer :: n1a, n1b, n2a, n2b, n2c, ndim, iroot


        call print_header()
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/chplus-1.0-olsen/onebody.inp',&
                           '/users/karthik/Desktop/CC_matlab_tests/chplus-1.0-olsen/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        call print_calc_params(sys)

        !call test_leftccsd(sys,fA,fB,vA,vB,vC)

        n1a = sys%Nunocc_a * sys%Nocc_a
        n1b = sys%Nunocc_b * sys%Nocc_b
        n2a = sys%Nunocc_a**2 * sys%Nocc_a**2
        n2b = sys%Nunocc_a * sys%Nunocc_b * sys%Nocc_a * sys%Nocc_b
        n2c = sys%Nunocc_b**2 * sys%Nocc_b**2
        ndim = n1a + n1b + n2a + n2b + n2c

        allocate(t1a(sys%Nunocc_a,sys%Nocc_a),t1b(sys%Nunocc_b,sys%Nocc_b),&
                 t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                 t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                 t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b), &
                 Rvec(ndim,nroot),Lvec(ndim,nroot+1),omega(nroot))

        call ccsd(sys,fA,fB,vA,vB,vC,ndiis,maxit,shift,tol,t1a,t1b,t2a,t2b,t2c,Ecorr)
        call hbar_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,H1A,H1B,H2A,H2B,H2C)
        call davidson_eomccsd(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                               t1a,t1b,t2a,t2b,t2c,ndim,nroot,Rvec,omega,tol,maxit,init_guess)

        do iroot = 0,nroot
           call left_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Lvec(:,iroot+1),Rvec(:,iroot),&
                   H1A,H1B,H2A,H2B,H2C,tol,ndiis,maxit,shift,iroot,omega)
        end do
        !deallocate(t1a,t1b,t2a,t2b,t2c)
end program main

