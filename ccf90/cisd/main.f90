program main

        use integral_types, only: e1int_t, e2int_t
        use integrals, only: get_integrals
        use system_types, only: sys_t
        use printing, only: print_calc_params, print_header
        use cisd, only: build_H_slater

        implicit none

        integer, parameter :: io_file = 201
        integer, parameter :: Norb = 20, Nfroz=0, Nelec=4
        type(e1int_t) :: fA, fB, zA, zB
        type(e2int_t) :: vA, vB, vC
        type(sys_t) :: sys

        open(unit=io_file,file='out.log',form='formatted',status='replace')

        call print_header(io_file)
        call get_integrals('/Users/karthik/Desktop/CC_matlab_tests/h4-dz-alpha-0.1/onebody.inp',&
                           '/users/karthik/Desktop/CC_matlab_tests/h4-dz-alpha-0.1/twobody.inp',&
                           Norb,Nfroz,Nelec,sys,vA,vB,vC,fA,fB,zA,zB,io_file)
        call print_calc_params(sys,io_file)

        call build_H_slater(sys,zA,zB,vA,vB,vC)

        close(io_file)


end program main
