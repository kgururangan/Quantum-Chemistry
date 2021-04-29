module cis_updates

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use einsum_module, only: einsum
        
        implicit none

        contains

                subroutine HC_1A(sys,fA,fB,vA,vB,vC,c1a,c1b,sigma_rs)

                ! Calculates <Phi_i^a | H_N C1 | Phi >
                !
                ! In: 
                !     sys - system derived type 
                !     fA, fB - Fock aa and bb onebody integral types
                !     vA, vB, vC - Antisymmetrized aa, ab, and bb twobody integral types
                !     c1a, c1b - CI amplitudes
                ! Out:
                !     sigma_rs - 1D array of alpha singles projection


                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: c1a(sys%Nunocc_a,sys%Nocc_a), c1b(sys%Nunocc_b,sys%Nocc_b)
                        real :: sigma(sys%Nunocc_a,sys%Nocc_a), Q(sys%Nunocc_a,sys%Nocc_a)
                        real, intent(out) :: sigma_rs(sys%Nunocc_a*sys%Nocc_a)

                        sigma = fA%uo
                        call einsum('ae,ei->ai',fA%uu,c1a,Q)
                        sigma = sigma + Q
                        call einsum('mi,am->ai',fA%oo,c1a,Q)
                        sigma = sigma - Q
                        call einsum('amie,em->ai',vA%uoou,c1a,Q)
                        sigma = sigma + Q
                        call einsum('amie,em->ai',vB%uoou,c1b,Q)
                        sigma = sigma + Q
                        sigma_rs = reshape(sigma,(/sys%Nocc_a*sys%Nunocc_a/)) 

                 end subroutine HC_1A

                 subroutine HC_1B(sys,fA,fB,vA,vB,vC,c1a,c1b,sigma_rs)

                ! Calculates <Phi_i~^a~ | H_N C1 | Phi >
                !
                ! In: 
                !     sys - system derived type 
                !     fA, fB - Fock aa and bb onebody integral types
                !     vA, vB, vC - Antisymmetrized aa, ab, and bb twobody integral types
                !     c1a, c1b - CI amplitudes
                ! Out:
                !     sigma_rs - 1D array of beta singles projection

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: c1a(sys%Nunocc_a,sys%Nocc_a), c1b(sys%Nunocc_b,sys%Nocc_b)
                        real :: sigma(sys%Nunocc_b,sys%Nocc_b), Q(sys%Nunocc_b,sys%Nocc_b)
                        real, intent(out) :: sigma_rs(sys%Nunocc_b*sys%Nocc_b)

                        sigma = fB%uo
                        call einsum('ae,ei->ai',fB%uu,c1b,Q)
                        sigma = sigma + Q
                        call einsum('mi,am->ai',fB%oo,c1b,Q)
                        sigma = sigma - Q
                        call einsum('amie,em->ai',vC%uoou,c1b,Q)
                        sigma = sigma + Q
                        call einsum('maei,em->ai',vB%ouuo,c1a,Q)
                        sigma = sigma + Q
                        sigma_rs = reshape(sigma,(/sys%Nocc_b*sys%Nunocc_b/))

                  end subroutine HC_1B

                  subroutine HC_matmat(sys,fA,fB,vA,vB,vC,cvec,ndim,crsz,sigma)

                ! Calculates <Phi_i^a | H_N C1 | Phi >
                ! In: 
                !     sys - system derived type 
                !     fA, fB - Fock aa and bb onebody integral types
                !     vA, vB, vC - Antisymmetrized aa, ab, and bb twobody integral types
                !     cvec - 2D array of Davidson search vectors (c1a and c1b)
                !     ndim - size of total singles projection
                !     crsz - number of guess vectors
                ! Out:
                !     sigma_rs - 2D array of total singles projection for each guess vector

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        integer, intent(in) :: ndim, crsz
                        real, intent(in) :: cvec(ndim,crsz)
                        real, intent(out) :: sigma(ndim,crsz)
                        real :: c1a(sys%Nunocc_a,sys%Nocc_a), c1b(sys%Nunocc_b,sys%Nocc_b), &
                                sigma_1a(sys%Nunocc_a*sys%Nocc_a), sigma_1b(sys%Nunocc_b*sys%Nocc_b)
                        integer :: j, n1a, n1b

                        n1a = sys%Nunocc_a * sys%Nocc_a
                        n1b = sys%Nunocc_b * sys%Nocc_b

                        do j = 1,crsz

                                c1a = reshape(cvec(1:n1a,j),shape(c1a))
                                c1b = reshape(cvec(n1a+1:n1a+n1b,j),shape(c1b))

                                call HC_1A(sys,fA,fB,vA,vB,vC,c1a,c1b,sigma_1a)
                                call HC_1B(sys,fA,fB,vA,vB,vC,c1a,c1b,sigma_1b)

                                sigma(1:n1a,j) = sigma_1a
                                sigma(n1a+1:n1a+n1b,j) = sigma_1b

                        end do

                   end subroutine HC_matmat


                        


end module cis_updates
