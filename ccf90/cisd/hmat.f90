module hmat

        use integral_types, only: e1int_t, e2int_t
        use excitation_types, only: excit_t
        use system_types, only: sys_t
        use slater, only: get_excitation, create_det
        use sort_module, only: argsort_int
        use permutils, only: permsign

        implicit none

        contains

            function get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc) result(val)

                    real :: val
                    type(sys_t), intent(in) :: sys
                    integer, intent(in) :: det1(:), det2(:)
                    type(e1int_t), intent(in) :: zA, zB
                    type(e2int_t), intent(in) :: vA, vB, vC

                    type(excit_t), intent(out) :: exc
                    integer :: N_int, i2, j2, a2, b2, idx(sys%Nelec/2)
                    integer, allocatable :: Int1(:,:), Int2(:,:)
                    real :: sgn1, sgn2, spinperm

                    spinperm = 1.0

                    N_int = floor(real(sys%Norb/64))+1
                    allocate(Int1(N_int,2),Int2(N_int,2))

                    call create_det(det1,N_int,Int1)
                    idx = argsort_int(det1(1:sys%Nelec-1:2))
                    call permsign(sgn1,idx)
                    idx = argsort_int(det1(2:sys%Nelec:2))
                    call permsign(sgn2,idx)
                    spinperm = sgn1 * sgn2

                    call create_det(det2,N_int,Int2)
                    idx = argsort_int(det2(1:sys%Nelec-1:2))
                    call permsign(sgn1,idx)
                    idx = argsort_int(det2(2:sys%Nelec:2))
                    call permsign(sgn2,idx)
                    spinperm = spinperm * sgn1 * sgn2

                    call get_excitation(Int1,Int2,N_int,exc)
                    deallocate(Int1,Int2)

                    val = 0.0

                    if (exc%nexcit == -1) then
                            val = 0.0
                    else
                            if (exc%nexcit == 2) then
                                    if (exc%nexcit_alpha == 2) then
                                            i2 = (2*exc%from_a(1)-1)
                                            j2 = (2*exc%from_a(2)-1)
                                            a2 = (2*exc%to_a(1)-1)
                                            b2 = (2*exc%to_a(2)-1)
                                    else if (exc%nexcit_alpha == 1) then
                                            i2 = (2*exc%from_a(1)-1)
                                            j2 = (2*exc%from_b(1))
                                            a2 = (2*exc%to_a(1)-1)
                                            b2 = (2*exc%to_b(1))
                                    else
                                            i2 = 2*exc%from_b(1)
                                            j2 = 2*exc%from_b(2)
                                            a2 = 2*exc%to_b(1)
                                            b2 = 2*exc%to_b(2)
                                    end if

                                    val = HBar_2(sys,det1,zA%full,zB%full,&
                                            vA%full,vB%full,vC%full,i2,j2,a2,b2,exc%phase)
                            end if

                            if (exc%nexcit == 1) then
                                    if (exc%nexcit_alpha == 1) then
                                            i2 = (2*exc%from_a(1)-1)
                                            a2 = (2*exc%to_a(1)-1)
                                    else
                                            i2 = 2*exc%from_b(1)
                                            a2 = 2*exc%to_b(1)
                                    end if

                                    val = HBar_1(sys,det1,zA%full,zB%full,&
                                            vA%full,vB%full,vC%full,i2,a2,exc%phase)

                            end if

                            if (exc%nexcit == 0) then
                                    val = HBar_0(sys,det1,zA%full,zB%full,&
                                            vA%full,vB%full,vC%full)
                            end if

                            val = val * spinperm
                            exc%phase = exc%phase * spinperm

                    end if

            end function get_matrix_element


            function Hbar_0(sys, occ_list, zA, zB, vA, vB, vC) result(hmatel)

                    real :: hmatel
                    type(sys_t), intent(in) :: sys
                    integer, intent(in) :: occ_list(:)
                    real, intent(in) :: zA(:,:), zB(:,:),&
                            vA(:,:,:,:), vB(:,:,:,:), vC(:,:,:,:)

                    integer :: i, j
                    integer :: iel, jel
                    real :: Eref_ele

                    hmatel = 0.0
                    Eref_ele = sys%Escf - sys%Vnuc

                    do iel = 1, sys%Nelec
                       i = occ_list(iel)
                       hmatel = hmatel + get_z(zA, zB, i, i)

                       do jel = 1, sys%Nelec
                          j = occ_list(jel)
                          hmatel = hmatel + 0.5*get_v(vA, vB, vC, i, j, i, j)
                       end do
                    end do

                    hmatel = hmatel - Eref_ele

            end function HBar_0

            function HBar_1(sys, occ_list, zA, zB, vA, vB, vC, i, a, phase) result(hmatel)

                    real :: hmatel
                    type(sys_t), intent(in) :: sys
                    integer, intent(in) :: occ_list(:)
                    integer, intent(in) :: i, a
                    real, intent(in) :: phase
                    real, intent(in) :: zA(:,:), zB(:,:),&
                            vA(:,:,:,:), vB(:,:,:,:), vC(:,:,:,:)

                    integer :: iel

                    hmatel = get_z(zA, zB, i, a)

                    do iel = 1, sys%Nelec
                       hmatel = hmatel + get_v(vA, vB, vC, i, occ_list(iel), a, occ_list(iel))
                    end do

                    hmatel = hmatel * phase

            end function HBar_1

            function HBar_2(sys, occ_list, zA, zB, vA, vB, vC, i, j, a, b, phase) result(hmatel)

                    real :: hmatel
                    type(sys_t), intent(in) :: sys
                    integer, intent(in) :: occ_list(:)
                    real, intent(in) :: zA(:,:), zB(:,:),&
                            vA(:,:,:,:), vB(:,:,:,:), vC(:,:,:,:)
                    real, intent(in) :: phase
                    integer, intent(in) :: i, j, a, b

                    hmatel = phase * get_v(vA, vB, vC, i, j, a, b)

            end function HBar_2
                
            function get_z(z_a, z_b, i, a) result(z_int)

                ! Get one-body matrix element <i|z|a>

                ! In:
                !    ints: system's integrals
                !    i: ith spin-orbital
                !    a: ath spin-orbital

                ! Out:
                !    z_int: one-body operator matrix element

                real :: z_int
                real, intent(in) :: z_a(:,:), z_b(:,:)
                integer, intent(in) :: i, a

                ! Spatial orbitals
                integer :: i_sp, a_sp

                i_sp = int((i + 1) / 2)
                a_sp = int((a + 1) / 2)

                ! Choose spin case
                if (mod(i, 2) == 0) then
                    z_int = z_b(a_sp, i_sp)
                else
                    z_int = z_a(a_sp, i_sp)
                endif

            end function get_z


            function get_v(v_aa, v_ab, v_bb, i, j, a, b) result(v_int)

                ! Get two-body matrix element <ij|v|ab>

                ! In:
                !    ints: system's integrals
                !    i: ith spin-orbital
                !    j: jth spin-orbital
                !    a: ath spin-orbital
                !    b: bth spin-orbital

                ! Out:
                !    v_int: two-body operator matrix element

                real :: v_int
                real, intent(in) :: v_aa(:,:,:,:), v_ab(:,:,:,:), v_bb(:,:,:,:)
                integer, intent(in) :: i, j, a, b
                integer :: dod
                integer :: i_sp, j_sp, a_sp, b_sp

                ! Spatial orbitals
                i_sp = int((i + 1) / 2)
                j_sp = int((j + 1) / 2)
                a_sp = int((a + 1) / 2)
                b_sp = int((b + 1) / 2)

                ! Total spin
                dod = mod(i,2) + mod(j,2) + mod(a,2) + mod(b,2)

                v_int = 0.0

                ! All spins are the same
                if (dod == 4) then
                    v_int = v_aa(a_sp, b_sp, i_sp, j_sp)
                else if (dod == 0) then
                    v_int = v_bb(a_sp, b_sp, i_sp, j_sp)
                else if (mod(i, 2) == mod(a, 2)) then
                    ! Bra and ket indices match spin
                    v_int = v_ab(a_sp, b_sp, i_sp, j_sp)
                else if (mod(i, 2) == mod(b, 2)) then
                    ! Bra and ket indices are flipped
                    v_int = -v_ab(b_sp, a_sp, i_sp, j_sp)
                endif

            end function get_v

end module hmat

