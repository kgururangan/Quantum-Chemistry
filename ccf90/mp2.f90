module mp2

        use integral_types, only: e1int_t, e2int_t
        use system_types, only: sys_t

        implicit none

        contains


                subroutine calculate_mp2(sys,fA,fB,vA,vB,vC,Emp2)
                        
                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(out) :: Emp2
                        real :: tmp
                        integer :: i, j, a, b

                        write(*,fmt=*) ''
                        write(*,fmt=*) 'MP2 ROUTINE'
                        write(*,fmt=*) ''

                        do i = 1,sys%Nocc_a
                           do j = i+1,sys%Nocc_a
                              do a = 1,sys%Nunocc_a
                                 do b = a+1,sys%Nunocc_a
                                    tmp = vA%oouu(i,j,a,b)*vA%uuoo(a,b,i,j)
                                    tmp = -1.0*tmp/(fA%uu(a,a)+fA%uu(b,b)-fA%oo(i,i)-fA%oo(j,j))
                                    Emp2 = Emp2 + tmp
                                 end do
                               end do
                            end do
                         end do

                         do i = 1,sys%Nocc_a
                            do j = 1,sys%Nocc_b
                               do a = 1,sys%Nunocc_a
                                  do b = 1,sys%Nunocc_b
                                     tmp = vB%oouu(i,j,a,b)*vB%uuoo(a,b,i,j)
                                     tmp = -1.0*tmp/(fA%uu(a,a)+fB%uu(b,b)-fA%oo(i,i)-fB%oo(j,j))
                                     Emp2 = Emp2 + tmp
                                  end do
                               end do
                            end do
                         end do

                         do i = 1,sys%Nocc_b
                            do j = i+1,sys%Nocc_b
                               do a = 1,sys%Nunocc_b
                                  do b = a+1,sys%Nunocc_b
                                     tmp = vC%oouu(i,j,a,b)*vC%uuoo(a,b,i,j)
                                     tmp = -1.0*tmp/(fB%uu(a,a)+fB%uu(b,b)-fB%oo(i,i)-fB%oo(j,j))
                                     Emp2 = Emp2 + tmp
                                  end do
                               end do
                            end do
                         end do

                         write(*,fmt=*) 'MP2 CORRELATION ENERGY (HARTREE) = ',Emp2
                         write(*,fmt=*) 'MP2 ENERGY (HARTREE) = ',Emp2+sys%Escf

                end subroutine calculate_mp2

 end module mp2
