program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use mp2, only: calculate_mp2
        use ccsd_module, only: ccsd 
        use printing, only: print_calc_params, print_header, print_L_amps, print_R_amps
        use cis, only: cis_mat
        use eomccsd, only: davidson_eomccsd
        use hbar, only: hbar_ccsd
        use leftccsd, only: left_ccsd
        use crcc23_module, only: crcc23
        use creomcc23_module, only: creomcc23

        implicit none

        integer, parameter :: io_file = 201
        integer, parameter :: Norb = 26, Nfroz=0, Nelec=6
        type(e1int_t) :: fA, fB, H1A, H1B
        type(e2int_t) :: vA, vB, vC, H2A, H2B, H2C
        type(sys_t) :: sys
        real :: Emp2, Ecorr
        integer, parameter :: ndiis = 6, nroot = 20, init_guess = 2
        real :: shift_cc, shift_lcc, tol
        real, allocatable :: t1a(:,:), t1b(:,:), t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:), &
                             Rvec(:,:), omega(:), Lvec(:,:), left_resid(:)
        integer :: n1a, n1b, n2a, n2b, n2c, ndim, iroot, maxit_cc, maxit_eom
        real :: E23A, E23B, E23C, E23D
        integer :: i, j, a, b

        maxit_cc = 2000
        maxit_eom = 50
        shift_cc = 0.0
        shift_lcc = 0.8
        tol = 1.0e-08

        open(unit=io_file,file='OUTPUT.log',status='replace',form='formatted')


        call print_header(io_file)
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests//chplus-1.0-olsen/onebody.inp',&
                           '/users/karthik/Desktop/CC_matlab_tests/chplus-1.0-olsen/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB,io_file)
        call print_calc_params(sys,io_file)

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

        

        call ccsd(sys,fA,fB,vA,vB,vC,ndiis,maxit_cc,shift_cc,tol,t1a,t1b,t2a,t2b,t2c,Ecorr,io_file)
        call hbar_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,H1A,H1B,H2A,H2B,H2C)
        call davidson_eomccsd(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,&
                               t1a,t1b,t2a,t2b,t2c,ndim,nroot,Rvec,omega,tol,maxit_eom,init_guess,io_file)

        allocate(left_resid(nroot+1))
        do iroot = 0,nroot
           call left_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Lvec(:,iroot+1),Rvec,&
                   H1A,H1B,H2A,H2B,H2C,tol,ndiis,maxit_cc,shift_lcc,iroot,omega,left_resid(iroot),io_file)
        end do
        if (left_resid(1) <= tol) then
                write(io_file,fmt=*) 'Ground state converged'
        else
                write(io_file,fmt=*) 'Ground state not converged'
        end if
        do iroot = 1,nroot
           if (left_resid(iroot+1) <= tol) then
                   write(io_file,fmt=*) 'Root - ',iroot,'e = ',omega(iroot),'CONVERGED'
           else
                   write(io_file,fmt=*) 'Root - ',iroot,'e = ',omega(iroot),'NOT CONVERGED'
           end if
        end do
        deallocate(left_resid)

          
        !call print_L_amps(sys,Lvec(:,6),1.0e-04)
        !call print_R_amps(sys,Rvec(:,5),1.0e-02)

        !call crcc23(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Lvec(:,1),H1A,H1B,H2A,H2B,H2C,E23A,E23B,E23C,E23D,io_file)
        do iroot = 1,nroot
                call creomcc23(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Lvec(:,iroot+1),Rvec(:,iroot),&
                        omega(iroot),H1A,H1B,H2A,H2B,H2C,E23A,E23B,E23C,E23D,io_file)
        end do

        deallocate(t1a,t1b,t2a,t2b,t2c,Rvec,Lvec,omega)


        close(io_file)

end program main

