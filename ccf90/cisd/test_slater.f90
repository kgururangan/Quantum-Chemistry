program main

        ! some notes:
        ! determinants are represented as an array of 64-bit integers
        ! a 64-bit integer is integer(kind=8)
        ! each integer(kind=8) has 64 slots, e.g. _ _ _ _ _ ..... _ _ _ (64 of them)
        ! each slot is occupied (1) or unoccupied (0)
        ! the actual integer number that results from that binary representation is what
        ! is stored in the array. E.g, if you have an alpha bitstring of 
        ! '01110101000' (Norb = 11), then the integer is 2^1 + 2^2 + 2^3 + 2^5 + 2^7 = 174
        ! and a beta bitstring of '11101010000' that corresponds to 2^0 + 2^1 + 2^2 + 2^4 + 2^6 = 87
        ! then the overall determinant is I = [174, 87]

        ! since each integer only has 64 bits, you can represent only 64 MO's per integer. Hence 
        ! if you have a system with 120 MO's, you will need two integers per spin and I will now look 
        ! like 
        ! I = [ I1(alpha), I1(beta);
        !       I2(alpha), I2(beta]

        use slater
        use const, only: i0_length

        implicit none

        integer(kind=8), allocatable :: I1(:,:), I2(:,:)
        type(excit_t) :: exc
        integer :: N_int, Nelec, Norb
        integer, allocatable :: D1(:), D2(:)

        print*,i0_length


        print*,'Test 1:'
        Nelec = 4
        Norb = 24
        allocate(D1(Nelec),D2(Nelec))
        N_int = floor(real(Norb/64))+1

        D1 = (/1, 2, 3, 4/)
        D2 = (/1, 2, 3, 6/)
        print*,'D1 = 1A  1B  2A  2B'
        print*,'D2 = 1A  1B  2A  3B'
        print*,'SHOULD BE 2B -> 3B'

        allocate(I1(N_int,2),I2(N_int,2))
        call create_det(D1,N_int,I1)
        call create_det(D2,N_int,I2)
        call get_excitation(I1,I2,N_int,exc)
        call print_excitation(exc)
        deallocate(I1,I2,D1,D2)

        print*,''

        print*,'Test 2:'
        Nelec = 4
        Norb = 24
        allocate(D1(Nelec),D2(Nelec))
        N_int = floor(real(Norb/64))+1

        D1 = (/1, 3, 6, 10/)
        D2 = (/1, 5, 6, 8/)
        print*,'D1 = 1A  2A  3B  5B'
        print*,'D2 = 1A  3B  3A  4B'
        print*,'SHOULD BE 2B -> 3B'

        allocate(I1(N_int,2),I2(N_int,2))
        call create_det(D1,N_int,I1)
        call create_det(D2,N_int,I2)
        call get_excitation(I1,I2,N_int,exc)
        call print_excitation(exc)
        deallocate(I1,I2,D1,D2)
       
        print*,'' 

        print*,'Test 3:'
        Nelec = 4
        Norb = 24
        allocate(D1(Nelec),D2(Nelec))
        N_int = floor(real(Norb/64))+1

        D1 = (/1, 2, 4, 5/)
        D2 = (/1, 2, 6, 7/)
        print*,'D1 = 1A  1B  2B  3A'
        print*,'D2 = 1A  1B  3B  4A'
        print*,'3A -> 4A   2B -> 3B'

        allocate(I1(N_int,2),I2(N_int,2))
        call create_det(D1,N_int,I1)
        call create_det(D2,N_int,I2)
        call get_excitation(I1,I2,N_int,exc)
        call print_excitation(exc)
        deallocate(I1,I2,D1,D2)

        print*,'' 

        print*,'Test 4:'
        Nelec = 4
        Norb = 24
        allocate(D1(Nelec),D2(Nelec))
        N_int = floor(real(Norb/64))+1

        D1 = (/2, 4, 9, 11/)
        D2 = (/2, 4, 11, 15/)
        print*,'D1 = 1B  2B  5A  6A'
        print*,'D2 = 1B  2B  6A  8A' 
        print*,'5A -> 8A'
        allocate(I1(N_int,2),I2(N_int,2))
        call create_det(D1,N_int,I1)
        call create_det(D2,N_int,I2)
        call get_excitation(I1,I2,N_int,exc)
        call print_excitation(exc)
        deallocate(I1,I2,D1,D2)

        print*,'' 

        print*,'Test 5:'
        Nelec = 8
        Norb = 8
        allocate(D1(Nelec),D2(Nelec))
        N_int = floor(real(Norb/64))+1

        D1 = (/2, 3, 4, 5, 6, 7, 8, 9/)
        D2 = (/1, 2, 4, 5, 6, 7, 8, 11/)
        print*,'D1 = 1B  2A  2B  3A  3B  4A  4B  5A'
        print*,'D2 = 1A  1B  2B  3A  3B  4A  4B  6A' 
        print*,'2A -> 1A    5A -> 6A'
        allocate(I1(N_int,2),I2(N_int,2))
        call create_det(D1,N_int,I1)
        call create_det(D2,N_int,I2)
        call get_excitation(I1,I2,N_int,exc)
        call print_excitation(exc)
        deallocate(I1,I2,D1,D2)


 end program main  
