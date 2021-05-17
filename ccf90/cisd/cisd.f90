module cisd
       
        use excitation_types, only: excit_t 
        use integral_types, only: e1int_t, e2int_t
        use system_types, only: sys_t
        use slater, only: create_det, get_excitation, print_excitation
        use hmat, only: get_matrix_element
        use sort_module, only: argsort_int, argsort
        use dgeev_module, only: eig

        ! something to look into would be the convention in hmat.f90 wher
        ! get_matrix_element is finding the excitation from D1 to D2 where
        ! <D1|H|D2>. In fact, we should be finding the excitation from D2 to D1 


        implicit none

        contains

                subroutine build_H_slater(sys,zA,zB,vA,vB,vC)

                        integer, parameter :: io_file = 205

                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: zA, zB
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, allocatable :: omega(:), CIvec(:,:)

                        real, allocatable :: wi(:)

                        type(excit_t) :: exc
                        integer, allocatable :: HF(:), det1(:), det2(:), idx(:)
                        integer :: n1a, n1b, n2a, n2b, n2c, i, j, a, b, c, d, k, l, idet, jdet, pos(5), ndim
                        real, allocatable :: Hmat(:,:), H1A1A(:,:), H1B1A(:,:), H1A1B(:,:), H1B1B(:,:), &
                                             H2A2A(:,:), H2A2B(:,:), H2A1A(:,:), H1A2A(:,:), &
                                             H2B2B(:,:), H2B2A(:,:), H2B2C(:,:), H2B1A(:,:), H2B1B(:,:), H1A2B(:,:), H1B2B(:,:), &
                                             H2C2C(:,:), H2C2B(:,:), H2C1B(:,:), H1B2C(:,:), H02A(:), H02B(:), H02C(:),&
                                             H2A0(:), H2B0(:), H2C0(:)

                        open(unit=io_file,file='tmpdebug.log',form='formatted',status='replace')

                        n1a = sys%Nunocc_a * sys%Nocc_a
                        n1b = sys%Nunocc_b * sys%Nocc_b
                        n2b = sys%Nunocc_a * sys%Nunocc_b * sys%Nocc_a * sys%Nocc_b

                        n2a = 0
                        do i = 1,sys%Nocc_a
                           do j = i+1,sys%Nocc_a
                              do a = 1,sys%Nunocc_a
                                 do b = a+1,sys%Nunocc_a
                                    n2a = n2a + 1
                                 end do
                              end do
                           end do
                        end do

                        n2c = 0
                        do i = 1,sys%Nocc_b
                           do j = i+1,sys%Nocc_b
                              do a = 1,sys%Nunocc_b
                                 do b = a+1,sys%Nunocc_b
                                    n2c = n2c + 1
                                 end do
                              end do
                           end do
                        end do

                        pos(1) = n1a+1
                        pos(2) = n1a+n1b+1
                        pos(3) = n1a+n1b+n2a+1
                        pos(4) = n1a+n1b+n2a+n2b+1
                        pos(5) = n1a+n2b+n2a+n2b+n2c+1

                        ndim = n1a+n1b+n2a+n2b+n2c

                        allocate(H02A(n2a), H02B(n2b), H02C(n2c), H2A0(n2a), H2B0(n2b), H2C0(n2c),&
                        H1A1A(n1a,n1a),&
                        H1A1B(n1a,n1b),H1B1A(n1b,n1a),&
                        H1B1B(n1b,n1b),&
                        H2A2A(n2a,n2a), H2A1A(n2a,n1a), H1A2A(n1a,n2a), H2A2B(n2a,n2b),&
                        H2B2B(n2b,n2b), H2B1A(n2b,n1a), H2B1B(n2b,n1b), H1A2B(n1a,n2b), &
                        H1B2B(n1b,n2b), H2B2A(n2b,n2a), H2B2C(n2b,n2c), &
                        H2C2C(n2c,n2c), H2C1B(n2c,n1b), H1B2C(n1b,n2c), H2C2B(n2c,n2b))

                        allocate(HF(sys%Nelec),det1(sys%Nelec),det2(sys%Nelec))
                        do i = 1,sys%Nelec
                           HF(i) = i
                        end do

                        ! H02A and H2A0
                        H02A = 0.0
                        H2A0 = 0.0
                        jdet = 0
                        do i = 1,sys%Nocc_a
                           do j = i+1,sys%Nocc_a
                              do a = 1,sys%Nunocc_a
                                 do b = a+1,sys%Nunocc_a
                                    jdet = jdet + 1
                                    H02A(jdet) = vA%uuoo(a,b,i,j)  
                                    H2A0(jdet) = vA%oouu(i,j,a,b)
                                 end do
                              end do
                           end do
                        end do              
                               
                        ! H02B and H2B0
                        H02B = 0.0
                        H2B0 = 0.0
                        jdet = 0
                        do i = 1,sys%Nocc_a
                           do j = 1,sys%Nocc_b
                              do a = 1,sys%Nunocc_a
                                 do b = 1,sys%Nunocc_b
                                    jdet = jdet + 1
                                    H02B(jdet) = vB%uuoo(a,b,i,j)  
                                    H2B0(jdet) = vB%oouu(i,j,a,b)
                                 end do
                              end do
                           end do
                        end do       

                        ! H02C and H2C0
                        H02C = 0.0
                        H2C0 = 0.0
                        jdet = 0
                        do i = 1,sys%Nocc_b
                           do j = i+1,sys%Nocc_b
                              do a = 1,sys%Nunocc_b
                                 do b = a+1,sys%Nunocc_b
                                    jdet = jdet + 1
                                    H02C(jdet) = vC%uuoo(a,b,i,j)  
                                    H2C0(jdet) = vC%oouu(i,j,a,b)
                                 end do
                              end do
                           end do
                        end do              


                        ! H1A1A
                        H1A1A = 0.0
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
                                 end do
                              end do
                           end do
                        end do

                        ! H1A1B
                        H1A1B = 0.0
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
                                 end do
                              end do
                           end do
                        end do

                        ! H1B1A
                        H1B1A = 0.0
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
                                 end do
                              end do
                           end do
                        end do

                        ! H1B1B
                        H1B1B = 0.0
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
                                 end do
                              end do
                           end do
                        end do

                        ! H2A1A
                        H2A1A = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = i+2,2*sys%Nocc_a-1,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = a+2,2*sys%Nunocc_a-1,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 1,2*sys%Nocc_a-1,2
                                       do c = 1,2*sys%Nunocc_a-1,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(k) = c+sys%Nelec
                                          H2A1A(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)

                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do


                        ! H2A2A
                        H2A2A = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = i+2,2*sys%Nocc_a-1,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = a+2,2*sys%Nunocc_a-1,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 1,2*sys%Nocc_a-1,2
                                       do l = k+2,2*sys%Nocc_a-1,2
                                          do c = 1,2*sys%Nunocc_a-1,2
                                             do d = c+2,2*sys%Nunocc_a-1,2
                                                  jdet = jdet + 1
                                                  det2 = HF
                                                  det2(k) = c+sys%Nelec
                                                  det2(l) = d+sys%Nelec
                                                  H2A2A(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H2A2B
                        H2A2B = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = i+2,2*sys%Nocc_a-1,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = a+2,2*sys%Nunocc_a-1,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 1,2*sys%Nocc_a-1,2
                                       do l = 2,2*sys%Nocc_b,2
                                          do c = 1,2*sys%Nunocc_a-1,2
                                             do d = 2,2*sys%Nunocc_b,2
                                                  jdet = jdet + 1
                                                  det2 = HF
                                                  det2(k) = c+sys%Nelec
                                                  det2(l) = d+sys%Nelec
                                                  H2A2B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H1A2A
                        H1A2A = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do a = 1,2*sys%Nunocc_a-1,2
                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec
                              do j = 1,2*sys%Nocc_a-1,2
                                 do k = j+2,2*sys%Nocc_a-1,2
                                    do b = 1,2*sys%Nunocc_a-1,2
                                       do c = b+2,2*sys%Nunocc_a-1,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(j) = b+sys%Nelec
                                          det2(k) = c+sys%Nelec
                                          H1A2A(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do


                        ! H2B2B
                        H2B2B = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = 2,2*sys%Nocc_b,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = 2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 1,2*sys%Nocc_a-1,2
                                       do l = 2,2*sys%Nocc_b,2
                                          do c = 1,2*sys%Nunocc_a-1,2
                                             do d = 2,2*sys%Nunocc_b,2
                                                  jdet = jdet + 1
                                                  det2 = HF
                                                  det2(k) = c+sys%Nelec
                                                  det2(l) = d+sys%Nelec
                                                  H2B2B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H2B2A
                        H2B2A = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = 2,2*sys%Nocc_b,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = 2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 1,2*sys%Nocc_a-1,2
                                       do l = k+2,2*sys%Nocc_a-1,2
                                          do c = 1,2*sys%Nunocc_a-1,2
                                             do d = c+2,2*sys%Nunocc_a-1,2
                                                  jdet = jdet + 1
                                                  det2 = HF
                                                  det2(k) = c+sys%Nelec
                                                  det2(l) = d+sys%Nelec
                                                  H2B2A(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H2B2C
                        H2B2C = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = 2,2*sys%Nocc_b,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = 2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 2,2*sys%Nocc_b,2
                                       do l = k+2,2*sys%Nocc_b,2
                                          do c = 2,2*sys%Nunocc_b,2
                                             do d = c+2,2*sys%Nunocc_b,2
                                                  jdet = jdet + 1
                                                  det2 = HF
                                                  det2(k) = c+sys%Nelec
                                                  det2(l) = d+sys%Nelec
                                                  H2B2C(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do


                        ! H2B1A
                        H2B1A = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = 2,2*sys%Nocc_b,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = 2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 1,2*sys%Nocc_a-1,2
                                       do c = 1,2*sys%Nunocc_a-1,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(k) = c+sys%Nelec
                                          H2B1A(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H2B1B
                        H2B1B = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do j = 2,2*sys%Nocc_b,2
                              do a = 1,2*sys%Nunocc_a-1,2
                                 do b = 2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 2,2*sys%Nocc_b,2
                                       do c = 2,2*sys%Nunocc_b,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(k) = c+sys%Nelec
                                          H2B1B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H1A2B
                        H1A2B = 0.0
                        idet = 0
                        do i = 1,2*sys%Nocc_a-1,2
                           do a = 1,2*sys%Nunocc_a-1,2
                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec
                              do j = 1,2*sys%Nocc_a-1,2
                                 do k = 2,2*sys%Nocc_b,2
                                    do b = 1,2*sys%Nunocc_a-1,2
                                       do c = 2,2*sys%Nunocc_b,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(j) = b+sys%Nelec
                                          det2(k) = c+sys%Nelec
                                          H1A2B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H1B2B
                        H1B2B = 0.0
                        idet = 0
                        do i = 2,2*sys%Nocc_b,2
                           do a = 2,2*sys%Nunocc_b,2
                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec
                              do j = 1,2*sys%Nocc_a-1,2
                                 do k = 2,2*sys%Nocc_b,2
                                    do b = 1,2*sys%Nunocc_a-1,2
                                       do c = 2,2*sys%Nunocc_b,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(j) = b+sys%Nelec
                                          det2(k) = c+sys%Nelec
                                          H1B2B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H2C1B
                        H2C1B = 0.0
                        idet = 0
                        do i = 2,2*sys%Nocc_b,2
                           do j = i+2,2*sys%Nocc_b,2
                              do a = 2,2*sys%Nunocc_b,2
                                 do b = a+2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 2,2*sys%Nocc_b,2
                                       do c = 2,2*sys%Nunocc_b,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(k) = c+sys%Nelec
                                          H2C1B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H2C2C
                        H2C2C = 0.0
                        idet = 0
                        do i = 2,2*sys%Nocc_b,2
                           do j = i+2,2*sys%Nocc_b,2
                              do a = 2,2*sys%Nunocc_b,2
                                 do b = a+2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 2,2*sys%Nocc_b,2
                                       do l = k+2,2*sys%Nocc_b,2
                                          do c = 2,2*sys%Nunocc_b,2
                                             do d = c+2,2*sys%Nunocc_b,2
                                                  jdet = jdet + 1
                                                  det2 = HF
                                                  det2(k) = c+sys%Nelec
                                                  det2(l) = d+sys%Nelec
                                                  H2C2C(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H2C2B
                        H2C2B = 0.0
                        idet = 0
                        do i = 2,2*sys%Nocc_b,2
                           do j = i+2,2*sys%Nocc_b,2
                              do a = 2,2*sys%Nunocc_b,2
                                 do b = a+2,2*sys%Nunocc_b,2
                                    idet = idet + 1
                                    jdet = 0
                                    det1 = HF
                                    det1(i) = a+sys%Nelec
                                    det1(j) = b+sys%Nelec
                                    do k = 1,2*sys%Nocc_a-1,2
                                       do l = 2,2*sys%Nocc_b,2
                                          do c = 1,2*sys%Nunocc_a-1,2
                                             do d = 2,2*sys%Nunocc_b,2
                                                  jdet = jdet + 1
                                                  det2 = HF
                                                  det2(k) = c+sys%Nelec
                                                  det2(l) = d+sys%Nelec
                                                  H2C2B(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do

                        ! H1B2C
                        H1B2C = 0.0
                        idet = 0
                        do i = 2,2*sys%Nocc_b,2
                           do a = 2,2*sys%Nunocc_b,2
                              idet = idet + 1
                              jdet = 0
                              det1 = HF
                              det1(i) = a+sys%Nelec
                              do j = 2,2*sys%Nocc_b,2
                                 do k = j+2,2*sys%Nocc_b,2
                                    do b = 2,2*sys%Nunocc_b,2
                                       do c = b+2,2*sys%Nunocc_b,2
                                          jdet = jdet + 1
                                          det2 = HF
                                          det2(j) = b+sys%Nelec
                                          det2(k) = c+sys%Nelec
                                          H1B2C(idet,jdet) = get_matrix_element(sys,det1,det2,zA,zB,vA,vB,vC,exc)
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do


                        allocate(Hmat(n1a+n1b+n2a+n2b+n2c+1,n1a+n1b+n2a+n2b+n2c+1))

                        Hmat = 0.0

                        Hmat(1,pos(2)+1:pos(3)) = H02A
                        Hmat(1,pos(3)+1:pos(4)) = H02B
                        Hmat(1,pos(4)+1:pos(5)) = H02C

                        Hmat(pos(2)+1:pos(3),1) = H2A0
                        Hmat(pos(3)+1:pos(4),1) = H2B0
                        Hmat(pos(4)+1:pos(5),1) = H2C0

                        Hmat(2:pos(1),2:pos(1)) = H1A1A
                        Hmat(2:pos(1),pos(1)+1:pos(2)) = H1A1B
                        Hmat(2:pos(1),pos(2)+1:pos(3)) = H1A2A
                        Hmat(2:pos(1),pos(3)+1:pos(4)) = H1A2B

                        Hmat(pos(1)+1:pos(2),2:pos(1)) = H1B1A
                        Hmat(pos(1)+1:pos(2),pos(1)+1:pos(2)) = H1B1B
                        Hmat(pos(1)+1:pos(2),pos(3)+1:pos(4)) = H1B2B
                        Hmat(pos(1)+1:pos(2),pos(4)+1:pos(5)) = H1B2C

                        Hmat(pos(2)+1:pos(3),2:pos(1)) = H2A1A
                        Hmat(pos(2)+1:pos(3),pos(2)+1:pos(3)) = H2A2A
                        Hmat(pos(2)+1:pos(3),pos(3)+1:pos(4)) = H2A2B

                        Hmat(pos(3)+1:pos(4),2:pos(1)) = H2B1A
                        Hmat(pos(3)+1:pos(4),pos(1)+1:pos(2)) = H2B1B
                        Hmat(pos(3)+1:pos(4),pos(2)+1:pos(3)) = H2B2A
                        Hmat(pos(3)+1:pos(4),pos(3)+1:pos(4)) = H2B2B
                        Hmat(pos(3)+1:pos(4),pos(4)+1:pos(5)) = H2B2C

                        Hmat(pos(4)+1:pos(5),pos(1)+1:pos(2)) = H2C1B
                        Hmat(pos(4)+1:pos(5),pos(3)+1:pos(4)) = H2C2B
                        Hmat(pos(4)+1:pos(5),pos(4)+1:pos(5)) = H2C2C

                        allocate(omega(n1a+n1b+n2a+n2b+n2c+1),wi(n1a+n1b+n2a+n2b+n2c+1),&
                                CIvec(n1a+n1b+n2a+n2b+n2c+1,n1a+n1b+n2a+n2b+n2c+1))
                        call eig(Hmat,CIvec,omega,wi)

                        allocate(idx(n1a+n1b+n2a+n2b+n2c+1))
                        idx = argsort(omega)
                        omega = omega(idx)
                        deallocate(idx)

                        do i = 0,ndim
                           print*,'Root ',i,'E = ',omega(i+1)+sys%Escf,'VEE = ',omega(i+1)-omega(1)
                        end do

                        deallocate(Hmat,H1A1A,H1A1B,H1B1A,H1B1B,&
                                H2A2A,H2A1A,H1A2A,H2B1B,H2B2B,H1A2B,H1B2B,H2B1A,H2C2C,H2C1B,H1B2C,&
                                H2A2B, H2B2A, H2B2C, H2C2B, &
                                CIvec,omega,wi)
                        deallocate(det1,det2,HF)

                        close(io_file)



                end subroutine build_H_slater
                                    
                                    
                                    
                            




 end module cisd
