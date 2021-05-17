module cis
       
        use excitation_types, only: excit_t 
        use integral_types, only: e1int_t, e2int_t
        use system_types, only: sys_t
        use slater, only: create_det, get_excitation, print_excitation
        use hmat, only: get_matrix_element
        use sort_module, only: argsort_int, argsort
        use permutils, only: permsign
        use dgeev_module, only: eig

        ! perhaps there is a phase associated with using split dets
        ! e.g.
        ! aaaaa | bbbbb vs abababababa

        implicit none

        contains

                subroutine build_H_slater(sys,fA,fB,zA,zB,vA,vB,vC)

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: zA, zB, fA, fB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, allocatable :: omega(:), CIvec(:,:)

                        real, allocatable :: wi(:)

                        type(excit_t) :: exc
                        integer, allocatable :: HF(:), det1(:), det2(:), &
                                                Int1(:,:), Int2(:,:), idx(:)
                        integer :: n1a, n1b, i, j, a, b, N_int, i2, j2, a2, b2, idet, jdet
                        real, allocatable :: Hmat(:,:), H1A1A(:,:), H1B1A(:,:), H1A1B(:,:), H1B1B(:,:), Hmat2(:,:)
                        real :: val, xsum
                        integer :: cnt1, cnt2


                        n1a = sys%Nunocc_a * sys%Nocc_a
                        n1b = sys%Nunocc_b * sys%Nocc_b

                        allocate(Hmat(n1a+n1b,n1a+n1b),Hmat2(n1a+n1b,n1a+n1b))
                        Hmat = 0.0
                        Hmat2 = 0.0

                      cnt1 = 0
                      cnt2 = 0
                      do i = 1,sys%Nocc_a
                         do a = 1,sys%Nunocc_a
                            cnt1 = cnt1 + 1
                            cnt2 = 0
                            do j = 1,sys%Nocc_a
                               do b = 1,sys%Nunocc_a
                                  cnt2 = cnt2 + 1
                                  val = vA%uoou(a,j,i,b)
                                  if (i == j) then
                                     val = val + fA%uu(a,b)
                                  end if
                                  if (a == b) then
                                     val = val - fA%oo(j,i)
                                  end if                                  
                                  Hmat2(cnt1,cnt2) = val
                               end do
                            end do
                         end do
                      end do
                        
                      cnt1 = 0
                      do i = 1,sys%Nocc_a
                         do a = 1,sys%Nunocc_a
                            cnt1 = cnt1 + 1
                            cnt2 = 0
                            do j = 1,sys%Nocc_b
                               do b = 1,sys%Nunocc_b
                                  cnt2 = cnt2 + 1
                                  val = vB%uoou(a,j,i,b)
                                  Hmat2(cnt1,n1a+cnt2) = val
                               end do
                            end do
                         end do
                      end do

                      cnt1 = 0
                      do i = 1,sys%Nocc_b
                         do a = 1,sys%Nunocc_b
                            cnt1 = cnt1 + 1
                            cnt2 = 0
                            do j = 1,sys%Nocc_a
                               do b = 1,sys%Nunocc_a
                                  cnt2 = cnt2 + 1
                                  val = vB%ouuo(j,a,b,i)
                                  Hmat2(n1a+cnt1,cnt2) = val
                               end do
                            end do
                         end do
                      end do

                      cnt1 = 0
                      do i = 1,sys%Nocc_b
                         do a = 1,sys%Nunocc_b
                            cnt1 = cnt1 + 1
                            cnt2 = 0
                            do j = 1,sys%Nocc_b
                               do b = 1,sys%Nunocc_b
                                  cnt2 = cnt2 + 1
                                  val = vC%uoou(a,j,i,b)
                                  if (i == j) then
                                     val = val + fB%uu(a,b)
                                  end if
                                  if (a == b) then
                                     val = val - fB%oo(j,i)
                                  end if
                                  Hmat2(n1a+cnt1,n1a+cnt2) = val
                               end do
                            end do
                         end do
                      end do



                        allocate(HF(sys%Nelec),det1(sys%Nelec),det2(sys%Nelec))
                        do i = 1,sys%Nelec
                           HF(i) = i
                        end do


                        ! Slater rule construction
                        allocate(H1A1A(n1a,n1a),H1A1B(n1a,n1b),H1B1A(n1b,n1a),H1B1B(n1b,n1b))

                        xsum = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do a = 1,2*sys%Nunocc_a-1,2

                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec

                              do j = 1,2*sys%Nocc_a-1,2
                                 do b = 1,2*sys%Nunocc_a-1,2

                                    jdet = jdet + 1
                                    det2 = HF
                                    det2(j) = b+sys%Nelec

                                    H1A1A(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)

                                    xsum = xsum + H1A1A(idet,jdet) - Hmat2(idet,jdet)


                                 end do
                              end do
                           end do
                        end do
                        print*,xsum


                        xsum = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do a = 1,2*sys%Nunocc_a-1,2

                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec

                              do j = 2,2*sys%Nocc_b,2
                                 do b = 2,2*sys%Nunocc_b,2

                                    jdet = jdet + 1
                                    det2 = HF
                                    det2(j) = b+sys%Nelec

                                    H1A1B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)

                                    xsum = xsum + (H1A1B(idet,jdet) - Hmat2(idet,n1a+jdet))


                                 end do
                              end do
                           end do
                        end do
                        print*,xsum


                        xsum = 0.0
                        idet = 0
                        do i = 2,2*sys%Nocc_b,2
                           do a = 2,2*sys%Nunocc_b,2

                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec

                              do j = 1,2*sys%Nocc_a-1,2
                                 do b = 1,2*sys%Nunocc_a-1,2

                                    jdet = jdet + 1
                                    det2 = HF
                                    det2(j) = b+sys%Nelec

                                    H1B1A(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)

                                    xsum = xsum + (H1B1A(idet,jdet) - Hmat2(n1a+idet,jdet))


                                 end do
                              end do
                           end do
                        end do
                        print*,xsum


                        xsum = 0.0
                        idet = 0
                        do i = 2,2*sys%Nocc_b,2
                           do a = 2,2*sys%Nunocc_b,2

                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec

                              do j = 2,2*sys%Nocc_b,2
                                 do b = 2,2*sys%Nunocc_b,2

                                    jdet = jdet + 1
                                    det2 = HF
                                    det2(j) = b+sys%Nelec

                                    H1B1B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)

                                    xsum = xsum + H1B1B(idet,jdet) - Hmat2(n1a+idet,n1a+jdet)


                                 end do
                              end do
                           end do
                        end do

                        print*,xsum


                        Hmat(1:n1a,1:n1a) = H1A1A
                        Hmat(n1a+1:n1a+n1b,1:n1a) = H1B1A
                        Hmat(1:n1a,n1a+1:n1a+n1b) = H1A1B
                        Hmat(n1a+1:n1a+n1b,n1a+1:n1a+n1b) = H1B1B

                        allocate(omega(n1a+n1b),wi(n1a+n1b),CIvec(n1a+n1b,n1a+n1b))
                        call eig(Hmat,CIvec,omega,wi)

                        allocate(idx(n1a+n1b))
                        idx = argsort(omega)
                        omega = omega(idx)
                        deallocate(idx)

                        do i = 1,20
                           print*,'Root ',i,' = ',omega(i)
                        end do

                        deallocate(Hmat,H1A1A,H1A1B,H1B1A,H1B1B,CIvec,omega,wi)
                        deallocate(det1,det2,HF)

                        allocate(omega(n1a+n1b),wi(n1a+n1b),CIvec(n1a+n1b,n1a+n1b))
                        call eig(Hmat2,CIvec,omega,wi)

                        allocate(idx(n1a+n1b))
                        idx = argsort(omega)
                        omega = omega(idx)
                        deallocate(idx)
                        print*,''
                        do i = 1,20
                           print*,'Root ',i,' = ',omega(i)
                        end do

                        deallocate(Hmat2,CIvec,omega,wi)


                end subroutine build_H_slater
                                    
                                    
                                    
                            




 end module cis
