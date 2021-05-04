module testing_module

        use einsum_module, only: einsum
        use integral_types, only: e1int_t, e2int_t
        use system_types, only: sys_t

        implicit none

        contains

                subroutine test_einsum(sys,fA,fB,vA,vB,vC)

                   type(e1int_t), intent(in) :: fA, fB
                   type(e2int_t), intent(in) :: vA, vB, vC
                   type(sys_t), intent(in) :: sys
                   real, parameter :: TOL = 1.0e-14
                   real :: val, xsum
                   integer :: i, j, a, b, m, n, e, f
                   real :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                           t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                           t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                           t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                   real, allocatable :: Z1(:,:), Z2(:,:,:,:) 

                   write(*,'(a/)') 'TESTING EINSTEIN SUMMATION ROUTINE'
                   call get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)

                   associate(nua=>sys%Nunocc_a, nub=>sys%Nunocc_b, noa=>sys%Nocc_a, nob=>sys%Nocc_b)

                           print*,'Case 1 : vA(amie)t2a(ebmj) -> Z(abij)'
                           allocate(Z2(nua,nua,noa,noa))
                           call einsum('amie,ebmj->abij',vA%uoou,t2a,Z2)
                           xsum = 0.0
                           do a = 1,nua
                              do b = 1,nua
                                 do i = 1,noa
                                    do j = 1,noa
                                       val = 0.0
                                       do m = 1,noa
                                          do e = 1,nua
                                             val = val + vA%uoou(a,m,i,e)*t2a(e,b,m,j)
                                          end do
                                       end do
                                       xsum = xsum + Z2(a,b,i,j) - val
                                    end do
                                 end do
                              end do
                           end do
                           deallocate(Z2)
                           call check_pass(xsum,TOL)


                           print*,'Case 2 : 0.5*vC(abef)t2c(efij) -> Z(abij)'
                           allocate(Z2(nub,nub,nob,nob))
                           call einsum('abef,efij->abij',0.5*vC%uuuu,t2c,Z2)
                           xsum = 0.0
                           do a = 1,nua
                              do b = 1,nub
                                 do i = 1,nob
                                    do j = 1,nob
                                       val = 0.0
                                       do e = 1,nub
                                          do f = e+1,nub
                                             val = val + vC%uuuu(a,b,e,f)*t2c(e,f,i,j)
                                          end do
                                       end do
                                       xsum = xsum + Z2(a,b,i,j) - val
                                    end do
                                 end do
                              end do
                           end do
                           deallocate(Z2)
                           call check_pass(xsum,TOL)

                           print*,'Case 3 : 0.5*vA(anef)*t2a(efin) -> Z(ai)'
                           allocate(Z1(nua,noa))
                           call einsum('anef,efin->ai',0.5*vA%uouu,t2a,Z1)
                           xsum = 0.0
                           do a = 1,nua
                              do i = 1,noa
                                 val = 0.0
                                 do e = 1,nua
                                    do f = e+1,nua
                                       do n = 1,noa
                                          val = val + vA%uouu(a,n,e,f)*t2a(e,f,i,n)
                                       end do
                                    end do
                                 end do
                                 xsum = xsum + Z1(a,i) - val
                              end do
                           end do
                           deallocate(Z1)
                           call check_pass(xsum,TOL)

                           print*, 'Case 4 : vB(anef)*t2b(efin) - > Z(ai)'
                           allocate(Z1(nua,noa))
                           call einsum('anef,efin->ai',vB%uouu,t2b,Z1)
                           xsum = 0.0
                           do a = 1,nua
                              do i = 1,noa
                                 val = 0.0
                                 do e = 1,nua
                                    do f = 1,nub
                                       do n = 1,nob
                                          val = val + vB%uouu(a,n,e,f)*t2b(e,f,i,n)
                                       end do
                                    end do
                                 end do
                                 xsum = xsum + Z1(a,i) - val
                              end do
                           end do
                           deallocate(Z1)
                           call check_pass(xsum,TOL)

                   end associate

                end subroutine test_einsum

                subroutine check_pass(xsum,tol)

                        real, intent(in) :: xsum, tol

                        if (abs(xsum) <= tol) then
                                write(*,'(a)') 'PASSED'
                        else
                                write(*,'(a)') 'FAILED'
                                print*,'Error = ',xsum
                        end if

                end subroutine check_pass


                subroutine get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)

                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        type(sys_t), intent(in) :: sys
                        real, intent(out) :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: val, denom, denom2
                        integer :: i, j, a, b, m, n, e, f


                        ! 2ND-ORDER APPROXIMATE t1a : <Phi_i^a|(R0 V_N R0 V_N)_C|Phi>
                        do i = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_a

                              val = 0.0
                              denom2 = fA%oo(i,i) - fA%uu(a,a)

                              do e = 1,sys%Nunocc_a
                                 do f = e+1,sys%Nunocc_a
                                    do n = 1,sys%Nocc_a
                                       denom = fA%oo(n,n)-fA%uu(e,e)-fA%uu(f,f)+fA%oo(i,i)
                                       val = val + vA%uouu(a,n,e,f)*vA%uuoo(e,f,i,n)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              do m = 1,sys%Nocc_a
                                 do n = m+1,sys%Nocc_a
                                    do f = 1,sys%Nunocc_a
                                       denom = fA%oo(m,m)+fA%oo(n,n)-fA%uu(f,f)-fA%uu(a,a)
                                       val = val - vA%ooou(m,n,i,f)*vA%uuoo(a,f,m,n)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              do e = 1,sys%Nunocc_a
                                 do f = 1,sys%Nunocc_b
                                    do n = 1,sys%Nocc_b
                                       denom = fB%oo(n,n)-fA%uu(e,e)-fB%uu(f,f)+fA%oo(i,i)
                                       val = val + vB%uouu(a,n,e,f)*vB%uuoo(e,f,i,n)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              do m = 1,sys%Nocc_a
                                 do n = 1,sys%Nocc_b
                                    do f = 1,sys%Nunocc_b
                                       denom = fA%oo(m,m)+fB%oo(n,n)-fB%uu(f,f)-fA%uu(a,a)
                                       val = val - vB%ooou(m,n,i,f)*vB%uuoo(a,f,m,n)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              t1a(a,i) = val

                           end do
                        end do

                        ! 2ND-ORDER APPROXIMATE t1b : <Phi_i~^a~|(R0 V_N R0 V_N)_C|Phi>
                        do i = 1,sys%Nocc_b
                           do a = 1,sys%Nunocc_b

                              val = 0.0
                              denom2 = fB%oo(i,i) - fB%uu(a,a)

                              do e = 1,sys%Nunocc_b
                                 do f = e+1,sys%Nunocc_b
                                    do n = 1,sys%Nocc_b
                                       denom = fB%oo(n,n)-fB%uu(e,e)-fB%uu(f,f)+fB%oo(i,i)
                                       val = val + vC%uouu(a,n,e,f)*vC%uuoo(e,f,i,n)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              do m = 1,sys%Nocc_b
                                 do n = m+1,sys%Nocc_b
                                    do f = 1,sys%Nunocc_b
                                       denom = fB%oo(m,m)+fB%oo(n,n)-fB%uu(f,f)-fB%uu(a,a)
                                       val = val - vC%ooou(m,n,i,f)*vC%uuoo(a,f,m,n)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              do e = 1,sys%Nunocc_b
                                 do f = 1,sys%Nunocc_a
                                    do n = 1,sys%Nocc_a
                                       denom = fA%oo(n,n)-fB%uu(e,e)-fA%uu(f,f)+fB%oo(i,i)
                                       val = val + vB%ouuu(n,a,f,e)*vB%uuoo(f,e,n,i)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              do m = 1,sys%Nocc_b
                                 do n = 1,sys%Nocc_a
                                    do f = 1,sys%Nunocc_a
                                       denom = fB%oo(m,m)+fA%oo(n,n)-fA%uu(f,f)-fB%uu(a,a)
                                       val = val - vB%oouo(n,m,f,i)*vB%uuoo(f,a,n,m)/(denom*denom2)
                                    end do
                                 end do
                              end do

                              t1b(a,i) = val

                           end do
                        end do

                        ! 1ST-ORDER APPROXIMATE t2a : <Phi_ij^ab | R0 V_N | Phi>
                        do i = 1,sys%Nocc_a
                           do j = i+1,sys%Nocc_a
                              do a = 1,sys%Nunocc_a
                                 do b = a+1,sys%Nunocc_a
                                        denom = fA%oo(i,i)+fA%oo(j,j)-fA%uu(a,a)-fA%uu(b,b)
                                        t2a(a,b,i,j) = vA%uuoo(a,b,i,j)/denom
                                        t2a(b,a,i,j) = -t2a(a,b,i,j)
                                        t2a(a,b,j,i) = -t2a(a,b,i,j)
                                        t2a(b,a,j,i) = t2a(a,b,i,j)
                                 end do
                              end do 
                           end do
                        end do

                        ! 1ST-ORDER APPROXIMATE t2b : <Phi_ij~^ab~ | R0 V_N | Phi>
                        do i = 1,sys%Nocc_a
                           do j = 1,sys%Nocc_b
                              do a = 1,sys%Nunocc_a
                                 do b = 1,sys%Nunocc_b
                                        denom = fA%oo(i,i)+fB%oo(j,j)-fA%uu(a,a)-fB%uu(b,b)
                                        t2b(a,b,i,j) = vB%uuoo(a,b,i,j)/denom
                                 end do
                              end do 
                           end do
                        end do

                        ! 1ST-ORDER APPROXIMATE t2c : <Phi_i~j~^a~b~ | R0 V_N | Phi>
                        do i = 1,sys%Nocc_b
                           do j = i+1,sys%Nocc_b
                              do a = 1,sys%Nunocc_b
                                 do b = a+1,sys%Nunocc_b
                                        denom = fB%oo(i,i)+fB%oo(j,j)-fB%uu(a,a)-fB%uu(b,b)
                                        t2c(a,b,i,j) = vC%uuoo(a,b,i,j)/denom
                                        t2c(b,a,i,j) = -t2c(a,b,i,j)
                                        t2c(a,b,j,i) = -t2c(a,b,i,j)
                                        t2c(b,a,j,i) = t2c(a,b,i,j)
                                 end do
                              end do 
                           end do
                        end do
                        
                end subroutine get_t_approx

                subroutine test_update_t1a(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: update_t1a

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, parameter :: shift = 0.0
                        integer :: a, i

                        call get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)
                        call update_t1a(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,shift)

                        do a = 1,sys%Nunocc_a
                           do i = 1,sys%Nocc_a
                              if (abs(t1a(a,i)) >= 1.0e-07) then
                                 print*,'t1a(',a,i,') = ',t1a(a,i)
                              end if
                           end do
                        end do

                end subroutine test_update_t1a

                subroutine test_update_t1b(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: update_t1b

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, parameter :: shift = 0.0
                        integer :: a, i

                        call get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)
                        call update_t1b(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,shift)

                        do a = 1,sys%Nunocc_b
                           do i = 1,sys%Nocc_b
                              if (abs(t1a(a,i)) >= 1.0e-07) then
                                 print*,'t1b(',a,i,') = ',t1b(a,i)
                              end if
                           end do
                        end do

                end subroutine test_update_t1b

                subroutine test_update_t2a(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: update_t2a
                        use ccsd_module, only: hbar_ccs_intermediates

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        type(e1int_t) :: H1A, H1B
                        type(e2int_t) :: H2A, H2B, H2C
                        real :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, parameter :: shift = 0.0
                        integer :: a, b, i, j

                        call get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)
                        call hbar_ccs_intermediates(sys,fA,fB,vA,vB,vC,t1a,t1b,H1A,H1B,H2A,H2B,H2C)
                        call update_t2a(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        do a = 1,sys%Nunocc_a
                           do b = a+1,sys%Nunocc_a
                              do i = 1,sys%Nocc_a
                                 do j = i+1,sys%Nocc_a
                                    if (abs(t2a(a,b,i,j)) >= 1e-07) then
                                       write(*,'(1x,a5,i0,i0,i0,i0,a5,f12.8)') 't2a(',a,b,i,j,') = ',t2a(a,b,i,j)
                                    end if
                                 end do
                              end do
                           end do
                        end do

                end subroutine test_update_t2a

                subroutine test_update_t2b(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: update_t2b
                        use ccsd_module, only: hbar_ccs_intermediates

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        type(e1int_t) :: H1A, H1B
                        type(e2int_t) :: H2A, H2B, H2C
                        real :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, parameter :: shift = 0.0
                        integer :: a, b, i, j

                        call get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)
                        call hbar_ccs_intermediates(sys,fA,fB,vA,vB,vC,t1a,t1b,H1A,H1B,H2A,H2B,H2C)
                        call update_t2b(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        do a = 1,sys%Nunocc_a
                           do b = 1,sys%Nunocc_b
                              do i = 1,sys%Nocc_a
                                 do j = 1,sys%Nocc_b
                                    if (abs(t2b(a,b,i,j)) >= 1e-07) then
                                       write(*,'(1x,a5,i0,i0,i0,i0,a5,f12.8)') 't2b(',a,b,i,j,') = ',t2b(a,b,i,j)
                                    end if
                                 end do
                              end do
                           end do
                        end do

                end subroutine test_update_t2b

                subroutine test_update_t2c(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: update_t2c
                        use ccsd_module, only: hbar_ccs_intermediates

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        type(e1int_t) :: H1A, H1B
                        type(e2int_t) :: H2A, H2B, H2C
                        real :: t1a(sys%Nunocc_a, sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_a), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, parameter :: shift = 0.0
                        integer :: a, b, i, j

                        call get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)
                        call hbar_ccs_intermediates(sys,fA,fB,vA,vB,vC,t1a,t1b,H1A,H1B,H2A,H2B,H2C)
                        call update_t2c(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C,t1a,t1b,t2a,t2b,t2c,shift)

                        do a = 1,sys%Nunocc_b
                           do b = a+1,sys%Nunocc_b
                              do i = 1,sys%Nocc_b
                                 do j = i+1,sys%Nocc_b
                                    if (abs(t2c(a,b,i,j)) >= 1e-07) then
                                       write(*,'(1x,a5,i0,i0,i0,i0,a5,f12.8)') 't2c(',a,b,i,j,') = ',t2c(a,b,i,j)
                                    end if
                                 end do
                              end do
                           end do
                        end do

                end subroutine test_update_t2c

                subroutine write_hbar_ccs(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: hbar_ccs_intermediates

                        integer, parameter :: io_h1a_ov = 1, io_h1a_oo = 2, io_h1a_vv = 3, & 
                                io_h1b_ov = 4, io_h1b_oo = 5, io_h1b_vv = 6, &
                                io_h2a_oooo = 7, io_h2a_vvvv = 8, io_h2a_vooo = 9, io_h2a_vvov = 10, io_h2a_voov = 11, &
                                io_h2b_oooo = 12, io_h2b_vvvv = 13, &
                                io_h2b_voov = 14, io_h2b_vovo = 15, io_h2b_ovov = 16, io_h2b_ovvo = 17, &
                                io_h2b_vooo = 18, io_h2b_ovoo = 19, io_h2b_vvov = 20, io_h2b_vvvo = 21, &
                                io_h2c_oooo = 22, io_h2c_vvvv = 23, io_h2c_vooo = 24, io_h2c_vvvo = 25, io_h2c_voov = 26
                        type(e1int_t), intent(in) :: fA, fB
                        type(e1int_t) :: h1A, h1B
                        type(e2int_t), intent(in) :: vA, vB, vC
                        type(e2int_t) :: h2A, h2B, h2C
                        type(sys_t), intent(in) :: sys
                        real, allocatable :: t1a(:,:), t1b(:,:), t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:)
                        integer :: a, b, i, j, m, n, e, f

                        ! HBAR CCS CONSTRUCTION
                        allocate(t1a(sys%Nunocc_a,sys%Nocc_a),t1b(sys%Nunocc_b,sys%Nocc_b),&
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a),&
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b),&
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b))

                        call get_t_approx(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)

                        call hbar_ccs_intermediates(sys,fA,fB,vA,vB,vC,t1a,t1b,h1a,h1b,h2a,h2b,h2c)

                        open(unit=io_h1a_ov,file='h1a_ov',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do i = 1,sys%Nocc_a
                              write(io_h1a_ov,fmt=*) h1a%ou(i,a)
                           end do
                        end do
                        close(io_h1a_ov)

                        open(unit=io_h1b_ov,file='h1b_ov',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_b
                           do i = 1,sys%Nocc_b
                              write(io_h1b_ov,fmt=*) h1b%ou(i,a)
                           end do
                        end do
                        close(io_h1b_ov)

                        open(io_h1a_vv,file='h1a_vv',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do b = 1,sys%Nunocc_a
                              write(io_h1a_vv,fmt=*) h1a%uu(a,b)
                           end do
                        end do
                        close(io_h1a_vv)

                        open(io_h1b_vv,file='h1b_vv',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_b
                           do b = 1,sys%Nunocc_b
                              write(io_h1b_vv,fmt=*) h1b%uu(a,b)
                           end do
                        end do
                        close(io_h1b_vv)

                        open(io_h1a_oo,file='h1a_oo',status='replace',form='formatted')
                        do i = 1,sys%Nocc_a
                           do j = 1,sys%Nocc_a
                              write(io_h1a_oo,fmt=*) h1a%oo(i,j)
                           end do
                        end do
                        close(io_h1a_oo)

                        open(io_h1b_oo,file='h1b_oo',status='replace',form='formatted')
                        do i = 1,sys%Nocc_b
                           do j = 1,sys%Nocc_b
                              write(io_h1b_oo,fmt=*) h1b%oo(i,j)
                           end do
                        end do
                        close(io_h1b_oo)

                        open(io_h2a_oooo,file='h2a_oooo',status='replace',form='formatted')
                        do m = 1,sys%Nocc_a
                           do n = 1,sys%Nocc_a
                              do i = 1,sys%Nocc_a
                                 do j = 1,sys%Nocc_a
                                    write(io_h2a_oooo,fmt=*) h2a%oooo(m,n,i,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2a_oooo)

                        open(io_h2a_vvvv,file='h2a_vvvv',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do b = 1,sys%Nunocc_a
                              do e = 1,sys%Nunocc_a
                                 do f = 1,sys%Nunocc_a
                                    write(io_h2a_vvvv,fmt=*) h2a%uuuu(a,b,e,f)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2a_vvvv)

                        open(io_h2a_vooo,file='h2a_vooo',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do m = 1,sys%Nocc_a
                              do i = 1,sys%Nocc_a
                                 do j = 1,sys%Nocc_a
                                    write(io_h2a_vooo,fmt=*) h2a%uooo(a,m,i,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2a_vooo)

                        open(io_h2a_vvov,file='h2a_vvov',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do b = 1,sys%Nunocc_a
                              do i = 1,sys%Nocc_a
                                 do e = 1,sys%Nunocc_a
                                    write(io_h2a_vvov,fmt=*) h2a%uuou(a,b,i,e)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2a_vvov)

                        open(io_h2a_voov,file='h2a_voov',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do m = 1,sys%Nocc_a
                              do i = 1,sys%Nocc_a
                                 do e = 1,sys%Nunocc_a
                                    write(io_h2a_voov,fmt=*) h2a%uoou(a,m,i,e)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2a_voov)

                        
                        open(io_h2b_oooo,file='h2b_oooo',status='replace',form='formatted')
                        do m = 1,sys%Nocc_a
                           do n = 1,sys%Nocc_b
                              do i = 1,sys%Nocc_a
                                 do j = 1,sys%Nocc_b
                                    write(io_h2b_oooo,fmt=*) h2b%oooo(m,n,i,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_oooo)

                        open(io_h2b_vvvv,file='h2b_vvvv',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do b = 1,sys%Nunocc_b
                              do e = 1,sys%Nunocc_a
                                 do f = 1,sys%Nunocc_b
                                    write(io_h2b_vvvv,fmt=*) h2b%uuuu(a,b,e,f)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_vvvv)

                        open(io_h2b_voov,file='h2b_voov',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do m = 1,sys%Nocc_b
                              do i = 1,sys%Nocc_a
                                 do e = 1,sys%Nunocc_b
                                    write(io_h2b_voov,fmt=*) h2b%uoou(a,m,i,e)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_voov)

                        open(io_h2b_vovo,file='h2b_vovo',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do m = 1,sys%Nocc_b
                              do e = 1,sys%Nunocc_a
                                 do i = 1,sys%Nocc_b
                                    write(io_h2b_vovo,fmt=*) h2b%uouo(a,m,e,i)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_vovo)

                        open(io_h2b_ovov,file='h2b_ovov',status='replace',form='formatted')
                        do m = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_b
                              do i = 1,sys%Nocc_a
                                 do e = 1,sys%Nunocc_b
                                    write(io_h2b_ovov,fmt=*) h2b%ouou(m,a,i,e)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_ovov)

                        open(io_h2b_ovvo,file='h2b_ovvo',status='replace',form='formatted')
                        do m = 1,sys%Nocc_a
                           do a = 1,sys%Nunocc_b
                              do e = 1,sys%Nunocc_a
                                 do i = 1,sys%Nocc_b
                                    write(io_h2b_ovvo,fmt=*) h2b%ouuo(m,a,e,i)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_ovvo)

                        open(io_h2b_ovoo,file='h2b_ovoo',status='replace',form='formatted')
                        do m = 1,sys%Nocc_a
                           do b = 1,sys%Nunocc_b
                              do i = 1,sys%Nocc_a
                                 do j = 1,sys%Nocc_b
                                    write(io_h2b_ovoo,fmt=*) h2b%ouoo(m,b,i,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_ovoo)

                        open(io_h2b_vooo,file='h2b_vooo',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do m = 1,sys%Nocc_b
                              do i = 1,sys%Nocc_a
                                 do j = 1,sys%Nocc_b
                                    write(io_h2b_vooo,fmt=*) h2b%uooo(a,m,i,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_vooo)

                        open(io_h2b_vvvo,file='h2b_vvvo',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do b = 1,sys%Nunocc_b
                              do e = 1,sys%Nunocc_a
                                 do j = 1,sys%Nocc_b
                                    write(io_h2b_vvvo,fmt=*) h2b%uuuo(a,b,e,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_vvvo)

                        open(io_h2b_vvov,file='h2b_vvov',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_a
                           do b = 1,sys%Nunocc_b
                              do i = 1,sys%Nocc_a
                                 do e = 1,sys%Nunocc_b
                                    write(io_h2b_vvov,fmt=*) h2b%uuou(a,b,i,e)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2b_vvov)

                        open(io_h2c_oooo,file='h2c_oooo',status='replace',form='formatted')
                        do m = 1,sys%Nocc_b
                           do n = 1,sys%Nocc_b
                              do i = 1,sys%Nocc_b
                                 do j = 1,sys%Nocc_b
                                    write(io_h2c_oooo,fmt=*) h2c%oooo(m,n,i,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2c_oooo)

                        open(io_h2c_vvvv,file='h2c_vvvv',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_b
                           do b = 1,sys%Nunocc_b
                              do e = 1,sys%Nunocc_b
                                 do f = 1,sys%Nunocc_b
                                    write(io_h2c_vvvv,fmt=*) h2c%uuuu(a,b,e,f)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2c_vvvv)

                        open(io_h2c_voov,file='h2c_voov',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_b
                           do m = 1,sys%Nocc_b
                              do i = 1,sys%Nocc_b
                                 do e = 1,sys%Nunocc_b
                                    write(io_h2c_voov,fmt=*) h2c%uoou(a,m,i,e)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2c_voov)

                        open(io_h2c_vooo,file='h2c_vooo',status='replace',form='formatted')
                        do m = 1,sys%Nocc_b
                           do a = 1,sys%Nunocc_b
                              do i = 1,sys%Nocc_b
                                 do j = 1,sys%Nocc_b
                                    write(io_h2c_vooo,fmt=*) h2c%uooo(a,m,i,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2c_vooo)

                        open(io_h2c_vvvo,file='h2c_vvvo',status='replace',form='formatted')
                        do a = 1,sys%Nunocc_b
                           do b = 1,sys%Nunocc_b
                              do e = 1,sys%Nunocc_b
                                 do j = 1,sys%Nocc_b
                                    write(io_h2c_vvvo,fmt=*) h2c%uuuo(a,b,e,j)
                                 end do
                              end do
                           end do
                        end do
                        close(io_h2c_vvvo)

                        deallocate(t1a,t1b,t2a,t2b,t2c)

                end subroutine write_hbar_ccs

                subroutine test_hbar_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c)
                        use hbar, only: hbar_ccsd

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                            t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        type(e1int_t) :: H1A, H1B
                        type(e2int_t) :: H2A, H2B, H2C
                        integer :: a, b, i, j, m, n, e, f
                        real, parameter :: tol = 1.0e-07

                        call hbar_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,H1A,H1B,H2A,H2B,H2C)

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)


                        print*,'TESTING H2B_ovov'
                        do m = 1,noa
                           do a = 1,nua
                              do i = 1,noa
                                 do e = 1,nua
                                    if (abs(H2B%ouou(m,a,i,e)) >= tol) then
                                            print*,'h2B_ovov(',m,a,i,e,') = ',H2B%ouou(m,a,i,e)
                                    end if
                                 end do
                              end do
                           end do
                        end do


                        end associate

                end subroutine test_hbar_ccsd
                
                subroutine test_eomccsd(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: ccsd
                        use hbar, only: hbar_ccsd
                        use eomccsd, only: HR_1A, HR_1B, HR_2A, HR_2B, HR_2C

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: r1a(sys%Nunocc_a,sys%Nocc_a), r1b(sys%Nunocc_b,sys%Nocc_b), &
                                r2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                r2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                r2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: HR1A(sys%Nunocc_a,sys%Nocc_a), HR1B(sys%Nunocc_b,sys%Nocc_b), &
                                HR2A(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                HR2B(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                HR2C(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        type(e1int_t) :: H1A, H1B
                        type(e2int_t) :: H2A, H2B, H2C
                        integer :: a, b, i, j
                        real, parameter :: tol = 1.0e-07, shift = 0.0, tol_cc = 1.0e-08
                        real :: Ecorr
                        integer, parameter :: ndiis = 5, maxit = 100

                        call ccsd(sys,fA,fB,vA,vB,vC,ndiis,maxit,shift,tol_cc,t1a,t1b,t2a,t2b,t2c,Ecorr)
                        call hbar_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,H1A,H1B,H2A,H2B,H2C)

                        call get_t_approx(sys,fA,fB,vA,vB,vC,r1a,r1b,r2a,r2b,r2c)

                        call HR_1A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR1A)
                        call HR_1B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR1B)
                        call HR_2A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR2A)
                        call HR_2B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR2B)
                        call HR_2C(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,r1a,r1b,r2a,r2b,r2c,HR2C)


                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)


                                print*,'TESTING R2B UPDATE'
                                do a =1,nua
                                   do b = 1,nub
                                      do i = 1,noa
                                         do j = 1,nob
                                            if (abs(HR2B(a,b,i,j)) >= tol) then
                                                    print*,'X2B(',a,b,i,j,') = ',HR2B(a,b,i,j)
                                            end if
                                         end do
                                      end do
                                   end do
                                end do



                        end associate

                end subroutine test_eomccsd

                subroutine test_leftccsd(sys,fA,fB,vA,vB,vC)

                        use ccsd_module, only: ccsd
                        use hbar, only: hbar_ccsd
                        use leftccsd, only: LH_1A, LH_1B, LH_2A, LH_2B, LH_2C

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real :: t1a(sys%Nunocc_a,sys%Nocc_a), t1b(sys%Nunocc_b,sys%Nocc_b), &
                                t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: l1a(sys%Nunocc_a,sys%Nocc_a), l1b(sys%Nunocc_b,sys%Nocc_b), &
                                l2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                l2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                l2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real :: LH1A(sys%Nunocc_a,sys%Nocc_a), LH1B(sys%Nunocc_b,sys%Nocc_b), &
                                LH2A(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a), &
                                LH2B(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b), &
                                LH2C(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        type(e1int_t) :: H1A, H1B
                        type(e2int_t) :: H2A, H2B, H2C
                        integer :: a, b, i, j
                        real, parameter :: tol = 1.0e-07, shift = 0.0, tol_cc = 1.0e-08
                        real :: Ecorr
                        integer, parameter :: ndiis = 5, maxit = 100, iroot = 0

                        call ccsd(sys,fA,fB,vA,vB,vC,ndiis,maxit,shift,tol_cc,t1a,t1b,t2a,t2b,t2c,Ecorr)
                        call hbar_ccsd(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,H1A,H1B,H2A,H2B,H2C)

                        call get_t_approx(sys,fA,fB,vA,vB,vC,l1a,l1b,l2a,l2b,l2c)

                        call LH_1A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH1A)
                        call LH_1B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH1B)
                        call LH_2A(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH2A)
                        call LH_2B(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH2B)
                        call LH_2C(sys,fA,fB,vA,vB,vC,H1A,H1B,H2A,H2B,H2C, &
                                t1a,t1b,t2a,t2b,t2c,l1a,l1b,l2a,l2b,l2c,iroot,LH2C)


                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)


                                print*,'TESTING L2B UPDATE'
                                do a =1,nua
                                   do b =1,nub
                                      do i = 1,noa
                                         do j = 1,nob
                                            if (abs(LH2B(a,b,i,j)) >= tol) then
                                                    print*,'X2B(',a,b,i,j,') = ',LH2B(a,b,i,j)
                                            end if
                                         end do
                                      end do
                                   end do
                                end do



                        end associate

                end subroutine test_leftccsd

 end module testing_module
