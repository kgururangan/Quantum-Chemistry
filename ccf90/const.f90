module const

    ! This module holds all constant values in ccq.
    ! For example, integer and real types, numerical constants,
    ! file units, etc.

    ! 32-bit (4-byte) integer.
    integer, parameter :: int_32 = selected_int_kind(6)
    ! 64-bit (8-byte) integer.
    integer, parameter :: int_64 = selected_int_kind(15)

    ! Integer select for determinants
    integer, parameter :: i0 = int_32
    ! Number of bits in an integer of type i0.
    ! Note that pgi 10.3 has a bug and returns 32 if bit_size(int(0,i0)) is used.
    integer, parameter :: i0_length = bit_size(0_i0)
    ! Index of the last bit in an integer of type i0.
    ! (Bit indexing in fortran ranges from 0 to bit_size-1.)
    integer, parameter :: i0_end = i0_length - 1

    ! 32-bit (4-byte) real
    integer, parameter :: sp = selected_real_kind(6, 37)
    ! 64-bit (8-byte) real
    integer, parameter :: dp = selected_real_kind(15, 307)
    ! 128-bit (16-byte) real
    integer, parameter :: qp = selected_real_kind(33, 4931)
    ! Selected precision
    integer, parameter :: p = dp



    ! Maximum parameters
    integer, parameter :: max_p_space = 50000000
    integer, parameter :: line_len = 255
    integer, parameter :: c_len = 500

    ! File units
    integer, parameter :: part_ints_a_unit = 130
    integer, parameter :: part_ints_b_unit = 131
    integer, parameter :: part_ints_c_unit = 132

    integer, parameter :: log_unit = 140
    integer, parameter :: config_unit = 141
    integer, parameter :: tmp_unit = 142
    integer, parameter :: cnf_unit = 147
    integer, parameter :: t_unit = 389
    integer, parameter :: t_vecs_unit = 390
    integer, parameter :: l_unit = 422
    integer, parameter :: l_vecs_unit = 423
    integer, parameter :: hbar_unit = 488
    integer, parameter :: walk_unit = 489
    integer, parameter :: sym_unit = 490

    ! Adapt for Ilias old scheme
    integer, parameter :: ta = 29
    integer, parameter :: tb = 30
    integer, parameter :: tc = 31
    integer, parameter :: td = 32
    integer, parameter :: te = 33

end module const