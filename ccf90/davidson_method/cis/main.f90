program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use printing, only: print_calc_params, print_header
        use davidson_module, only: davidson_cis

        implicit none

        integer, parameter :: Norb = 76, Nfroz=4, Nelec=28
        type(e1int_t) :: fA, fB
        type(e2int_t) :: vA, vB, vC
        type(sys_t) :: sys
        integer :: ndim, maxit, nroot
        real, allocatable :: Cmat(:,:), omega(:)
        real :: tol
        

        call print_header()
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/onebody.inp',&
                           '/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        call print_calc_params(sys)

        ndim = sys%Nocc_a * sys%Nunocc_a + sys%Nocc_b * sys%Nunocc_b
        maxit = 300
        nroot = 30
        tol = 1e-08
        allocate(Cmat(ndim,nroot),omega(nroot))
        call davidson_cis(sys,fA,fB,vA,vB,vC,ndim,nroot,Cmat,omega,tol,maxit)


end program main

