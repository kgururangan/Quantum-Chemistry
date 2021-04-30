module cc_energy

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use einsum_module, only: einsum
        use dgemm_module, only: gemm
        
        implicit none

        contains

                subroutine calc_cc_energy(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Ecorr)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, intent(out) :: Ecorr
                        integer :: m, n, e, f
                        real :: E1A, E1B, E1A1A, E1A1B, E1B1B, E2A, E2B, E2C

                        E1A = 0.0
                        do m = 1,sys%Nocc_a
                           do e = 1,sys%Nunocc_a
                              E1A = E1A + fA%ou(m,e)*t1a(e,m)
                           end do
                        end do

                        E1B = 0.0
                        do m = 1,sys%Nocc_b
                           do e = 1,sys%Nunocc_b
                              E1B = E1B + fB%ou(m,e)*t1b(e,m)
                           end do
                        end do

                        E2A = 0.0
                        E1A1A = 0.0
                        do m = 1,sys%Nocc_a
                           do n = 1,sys%Nocc_a
                              do e = 1,sys%Nunocc_a
                                 do f = 1,sys%Nunocc_a
                                    E2A = E2A + 0.25*vA%oouu(m,n,e,f)*t2a(e,f,m,n)
                                    E1A1A = E1A1A + 0.5*vA%oouu(m,n,e,f)*t1a(e,m)*t1a(f,n)
                                 end do
                              end do
                           end do
                        end do

                        E2B = 0.0
                        E1A1B = 0.0
                        do m = 1,sys%Nocc_a
                           do n = 1,sys%Nocc_b
                              do e = 1,sys%Nunocc_a
                                 do f = 1,sys%Nunocc_b
                                    E2B = E2B + vB%oouu(m,n,e,f)*t2b(e,f,m,n)
                                    E1A1B = E1A1B + vB%oouu(m,n,e,f)*t1a(e,m)*t1b(f,n)
                                 end do
                              end do
                           end do
                        end do

                        E2C = 0.0
                        E1B1B = 0.0
                        do m = 1,sys%Nocc_b
                           do n = 1,sys%Nocc_b
                              do e = 1,sys%Nunocc_b
                                 do f = 1,sys%Nunocc_b
                                    E2C = E2C + 0.25*vC%oouu(m,n,e,f)*t2c(e,f,m,n)
                                    E1B1B = E1B1B + 0.5*vC%oouu(m,n,e,f)*t1b(e,m)*t1b(f,n)
                                 end do
                              end do
                           end do
                        end do

                        Ecorr = E1A + E1B + E1A1A + E1A1B + E1B1B + E2A + E2B + E2C

                end subroutine calc_cc_energy

end module cc_energy
                        
                        


