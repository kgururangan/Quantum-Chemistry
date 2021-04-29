program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use printing, only: print_calc_params, print_header
        use davidson_module, only: davidson_cis

        implicit none

        integer, parameter :: Norb = 8, Nfroz=0, Nelec=8
        type(e1int_t) :: fA, fB
        type(e2int_t) :: vA, vB, vC
        type(sys_t) :: sys
        integer :: ndim, maxit, nroot
        real, allocatable :: Cmat(:,:), omega(:)
        real :: tol
        

        call print_header()
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/h8-mbs-0.0001/onebody.inp',&
                           '/users/karthik/Desktop/CC_matlab_tests/h8-mbs-0.0001/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB)
        call print_calc_params(sys)

        ndim = sys%Nocc_a * sys%Nunocc_a + sys%Nocc_b * sys%Nunocc_b
        maxit = 300
        nroot = 10
        tol = 1e-06
        allocate(Cmat(ndim,nroot),omega(nroot))
        call davidson_cis(sys,fA,fB,vA,vB,vC,ndim,nroot,Cmat,omega,tol,maxit)


end program main

