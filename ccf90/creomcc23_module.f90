module creomcc23_module

        use dgemm_module, only: dot
        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use permutils, only: reorder_stripe
        use einsum_module, only: einsum

        implicit none

        contains
                subroutine creomcc23(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,L,R,omega,H1A,H1B,H2A,H2B,H2C,E23A,E23B,E23C,E23D,io)

                        use cc_energy, only: calc_cc_energy

                        integer, intent(in) :: io
                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(:,:), t1b(:,:), t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:),&
                                            L(:), R(:), omega
                        real, intent(out) :: E23A, E23B, E23C, E23D
                        real, allocatable :: D3A_V(:,:,:), D3A_O(:,:,:), D3B_V(:,:,:), D3B_O(:,:,:), &
                                             D3C_V(:,:,:), D3C_O(:,:,:), D3D_V(:,:,:), D3D_O(:,:,:), &
                                             l1a(:,:), l1b(:,:), l2a(:,:,:,:), l2b(:,:,:,:), l2c(:,:,:,:),&
                                             r1a(:,:), r1b(:,:), r2a(:,:,:,:), r2b(:,:,:,:), r2c(:,:,:,:)
                        real, allocatable :: I1_d1(:,:,:,:), I2_d1(:,:,:,:), I3_d1(:,:,:,:), I1_d2(:,:,:,:),&
                                             I2_d2(:,:,:,:), I3_d2(:,:,:,:),&
                                             I1_d3(:,:,:,:), I1_d4(:,:,:,:), I1_d5(:,:,:,:), I2_d5(:,:,:,:),&
                                             I3_d5(:,:,:,:), I1_d6(:,:,:,:), I2_d6(:,:,:,:), I3_d6(:,:,:,:),&
                                             I1_d9(:,:), I2_d9(:,:), I3_d9(:,:,:,:), Q(:,:,:,:), &
                                             x1(:,:,:), x2(:,:,:), x3(:,:,:), x4(:,:,:),&
                                             y1(:,:,:), y2(:,:,:), y3(:,:,:), y4(:,:,:),&
                                             xd1(:,:,:), xd2(:,:,:), xd3(:,:,:), xd4(:,:,:), xd5(:,:,:), xd6(:,:,:),&
                                             xd7(:,:,:), xd8(:,:,:), xd9(:,:,:)
                        integer :: a, b, c, i, j, k, n1a, n1b, n2a, n2b, n2c
                        real :: MM23, L3, denom, val, z, Ecorr, xsum

                        write(io,fmt=*) ''
                        write(io,fmt=*) '++++++++++++++++++++ CR-EOMCC(2,3) ROUTINE +++++++++++++++++++++'
                        write(io,fmt=*) ''
                        
                        ! initialize
                        E23A = 0.0
                        E23B = 0.0
                        E23C = 0.0
                        E23D = 0.0

                        xsum = 0.0

                        associate(noa=>sys%Nocc_a,nob=>sys%Nocc_b,nua=>sys%Nunocc_a,nub=>sys%Nunocc_b)

                                n1a = nua*noa
                                n1b = nub*nob
                                n2a = nua**2 * noa**2
                                n2b = nua*nua*noa*nob
                                n2c = nub**2 * nob**2

                                allocate(l1a(nua,noa),l1b(nub,nob),l2a(nua,nua,noa,noa),l2b(nua,nub,noa,nob),l2c(nub,nub,nob,nob),&
                                         r1a(nua,noa),r1b(nub,nob),r2a(nua,nua,noa,noa),r2b(nua,nub,noa,nob),r2c(nub,nub,nob,nob))

                                l1a = reshape(L(1:n1a),(/nua,noa/))
                                l1b = reshape(L(n1a+1:n1a+n1b),(/nub,nob/))
                                l2a = reshape(L(n1a+n1b+1:n1a+n1b+n2a),(/nua,nua,noa,noa/))
                                l2b = reshape(L(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b),(/nua,nub,noa,nob/))
                                l2c = reshape(L(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c),(/nub,nub,nob,nob/))

                                r1a = reshape(R(1:n1a),(/nua,noa/))
                                r1b = reshape(R(n1a+1:n1a+n1b),(/nub,nob/))
                                r2a = reshape(R(n1a+n1b+1:n1a+n1b+n2a),(/nua,nua,noa,noa/))
                                r2b = reshape(R(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b),(/nua,nub,noa,nob/))
                                r2c = reshape(R(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c),(/nub,nub,nob,nob/))

                                ! Be sure to enforce that LR = 1 explicitly, even if it seems like it's 
                                ! close to 1 after the left-EOMCC routine.
                                call enforce_LR(l1a,l1b,l2a,l2b,l2c,r1a,r1b,r2a,r2b,r2c)

                                ! get 3-body HBar triples diagonal
                                call triples_3body_diagonal(sys,vA,vB,vC,t2a,t2b,t2c,&
                                        D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O)

                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!!!! MM(2,3)A CORRECTION !!!!!
                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                ! Diagram 1 intermediates
                                allocate(I1_d1(nua,nua,noa,nua),I2_d1(nua,nua,noa,nua))
                                call einsum('amie,cm->acie',H2A%uoou,r1a,I1_d1)
                                call reorder_stripe(4,shape(I1_d1),size(I1_d1),'2134',I1_d1,I2_d1)
                                I1_d1 = -1.0*I1_d1 + I2_d1
                                deallocate(I2_d1)
                                allocate(I2_d1(noa,nua,noa,noa))
                                call einsum('nmjk,cm->ncjk',H2A%oooo,r1a,I2_d1)
                                I2_d1 = -1.0*I2_d1

                                ! Diagram 2 intermediates
                                allocate(I1_d2(nua,noa,noa,noa),I2_d2(nua,noa,noa,noa))
                                call einsum('amie,ej->amij',H2A%uoou,r1a,I1_d2)
                                call reorder_stripe(4,shape(I1_d2),size(I1_d2),'1243',I1_d2,I2_d2)
                                I1_d2= I1_d2 - I2_d2
                                deallocate(I2_d2)
                                allocate(I2_d2(nua,nua,noa,nua))
                                call einsum('abfe,fi->abie',H2A%uuuu,r1a,I2_d2)

                                ! Diagram 3 intermediates
                                allocate(I1_d3(nua,noa,noa,noa))
                                call einsum('amef,efij->amij',H2A%uouu,r2a,I1_d3)
                                I1_d3 = 0.5*I1_d3

                                ! Diagram 4 intermediates
                                allocate(I1_d4(nua,nua,noa,nua))
                                call einsum('mnie,abmn->abie',H2A%ooou,r2a,I1_d4)
                                I1_d4 = 0.5*I1_d4

                                ! Diagram 5 intermediates
                                allocate(I1_d5(nua,nua,noa,nua),I2_d5(nua,nua,noa,nua))
                                call einsum('bmfe,aeim->abif',H2A%uouu,r2a,I1_d5)
                                call reorder_stripe(4,shape(I1_d5),size(I1_d5),'2134',I1_d5,I2_d5)
                                I1_d5 = I1_d5 - I2_d5
                                deallocate(I2_d5)
                                allocate(I2_d5(nua,noa,noa,noa),I3_d5(nua,noa,noa,noa))
                                call einsum('nmje,cekm->cnkj',H2A%ooou,r2a,I2_d5)
                                call reorder_stripe(4,shape(I2_d5),size(I2_d5),'1243',I2_d5,I3_d5)
                                I2_d5 = I2_d5 - I3_d5
                                deallocate(I3_d5)

                                ! Diagram 6 intermediates
                                allocate(I1_d6(nua,nua,noa,nua),I2_d6(nua,nua,noa,nua))
                                call einsum('bmfe,aeim->abif',H2B%uouu,r2b,I1_d6)
                                call reorder_stripe(4,shape(I1_d6),size(I1_d6),'2134',I1_d6,I2_d6)
                                I1_d6 = I1_d6 - I2_d6
                                deallocate(I2_d6)
                                allocate(I2_d6(nua,noa,noa,noa),I3_d6(nua,noa,noa,noa))
                                call einsum('nmje,cekm->cnkj',H2B%ooou,r2b,I2_d6)
                                call reorder_stripe(4,shape(I2_d6),size(I2_d6),'1243',I2_d6,I3_d6)
                                I2_d6 = I2_d6 - I3_d6
                                deallocate(I3_d6)

                                ! Diagram 9 intermediates
                                allocate(I1_d9(nua,noa),I2_d9(nua,noa))
                                call einsum('mnef,fn->me',vA%oouu,r1a,I1_d9)
                                call einsum('mnef,fn->me',vB%oouu,r1b,I2_d9)
                                I1_d9 = I1_d9 + I2_d9
                                deallocate(I2_d9)
                                allocate(I3_d9(noa,nua,noa,noa))
                                call einsum('me,ecjk->mcjk',I1_d9,t2a,I3_d9)
                                deallocate(I1_d9)

                                allocate(xd1(nua,nua,nua),xd2(nua,nua,nua),xd3(nua,nua,nua),&
                                         xd4(nua,nua,nua),xd5(nua,nua,nua),&
                                         xd7(nua,nua,nua),xd9(nua,nua,nua),&
                                         x1(nua,nua,nua), x2(nua,nua,nua),&
                                         y1(nua,nua,nua),y2(nua,nua,nua),y3(nua,nua,nua))

                                do i = 1,noa
                                   do j = i+1,noa
                                      do k = j+1,noa

                                         ! Approximate L3A
                                         y3 = 0.0
                                         call einsum('eba,ec->abc',H2A%uouu(:,i,:,:),l2a(:,:,j,k),y1)
                                         call einsum('mc,abm->abc',H2A%ooou(j,k,:,:),l2a(:,:,i,:),y2)
                                         y3 = y3 + (y1 - y2)                                 
                                         call einsum('eba,ec->abc',H2A%uouu(:,j,:,:),l2a(:,:,i,k),y1)
                                         call einsum('mc,abm->abc',H2A%ooou(i,k,:,:),l2a(:,:,j,:),y2)
                                         y3 = y3 - (y1 - y2)
                                         call einsum('eba,ec->abc',H2A%uouu(:,k,:,:),l2a(:,:,j,i),y1)
                                         call einsum('mc,abm->abc',H2A%ooou(j,i,:,:),l2a(:,:,k,:),y2)
                                         y3 = y3 - (y1 - y2)

                                         ! Diagram 1
                                         xd1 = 0.0
                                         call einsum('abe,ec->abc',I1_d1(:,:,i,:),t2a(:,:,j,k),x1)
                                         call einsum('nc,abn->abc',I2_d1(:,:,j,k),t2a(:,:,i,:),x2)
                                         xd1 = xd1 + (x1 - x2)
                                         call einsum('abe,ec->abc',I1_d1(:,:,j,:),t2a(:,:,i,k),x1)
                                         call einsum('nc,abn->abc',I2_d1(:,:,i,k),t2a(:,:,j,:),x2)
                                         xd1 = xd1 - (x1 - x2)
                                         call einsum('abe,ec->abc',I1_d1(:,:,k,:),t2a(:,:,j,i),x1)
                                         call einsum('nc,abn->abc',I2_d1(:,:,j,i),t2a(:,:,k,:),x2)
                                         xd1 = xd1 - (x1 - x2)

                                         ! Diagram 2
                                         xd2 = 0.0
                                         call einsum('am,bcm->abc',I1_d2(:,:,i,j),t2a(:,:,:,k),x1)
                                         call einsum('cbe,ae->abc',I2_d2(:,:,k,:),t2a(:,:,i,j),x2)
                                         xd2 = xd2 + (x2 - x1)
                                         call einsum('am,bcm->abc',I1_d2(:,:,k,j),t2a(:,:,:,i),x1)
                                         call einsum('cbe,ae->abc',I2_d2(:,:,i,:),t2a(:,:,k,j),x2)
                                         xd2 = xd2 - (x2 - x1)
                                         call einsum('am,bcm->abc',I1_d2(:,:,i,k),t2a(:,:,:,j),x1)
                                         call einsum('cbe,ae->abc',I2_d2(:,:,j,:),t2a(:,:,i,k),x2)
                                         xd2 = xd2 - (x2 - x1)

                                         ! Diagram 3
                                         xd3 = 0.0
                                         call einsum('am,bcm->abc',I1_d3(:,:,i,j),t2a(:,:,:,k),x1)
                                         xd3 = xd3 - x1
                                         call einsum('am,bcm->abc',I1_d3(:,:,k,j),t2a(:,:,:,i),x1)
                                         xd3 = xd3 + x1
                                         call einsum('am,bcm->abc',I1_d3(:,:,i,k),t2a(:,:,:,j),x1)
                                         xd3 = xd3 + x1

                                         ! Diagram 4
                                         xd4 = 0.0
                                         call einsum('abe,ec->abc',I1_d4(:,:,i,:),t2a(:,:,j,k),x1)
                                         xd4 = xd4 + x1
                                         call einsum('abe,ec->abc',I1_d4(:,:,j,:),t2a(:,:,i,k),x1)
                                         xd4 = xd4 - x1
                                         call einsum('abe,ec->abc',I1_d4(:,:,k,:),t2a(:,:,j,i),x1)
                                         xd4 = xd4 - x1

                                         ! Diagram 5 and 6
                                         xd5 = 0.0
                                         call einsum('abf,fc->abc',I1_d5(:,:,i,:)+I1_d6(:,:,i,:),t2a(:,:,j,k),x1)
                                         call einsum('cn,abn->abc',I2_d5(:,:,k,j)+I2_d6(:,:,k,j),t2a(:,:,i,:),x2)
                                         xd5 = xd5 + (x1 - x2)
                                         call einsum('abf,fc->abc',I1_d5(:,:,j,:)+I1_d6(:,:,j,:),t2a(:,:,i,k),x1)
                                         call einsum('cn,abn->abc',I2_d5(:,:,k,i)+I2_d6(:,:,k,i),t2a(:,:,j,:),x2)
                                         xd5 = xd5 - (x1 - x2)
                                         call einsum('abf,fc->abc',I1_d5(:,:,k,:)+I1_d6(:,:,k,:),t2a(:,:,j,i),x1)
                                         call einsum('cn,abn->abc',I2_d5(:,:,i,j)+I2_d6(:,:,i,j),t2a(:,:,k,:),x2)
                                         xd5 = xd5 - (x1 - x2)

                                         ! Diagram 7 and 8
                                         xd7 = 0.0
                                         call einsum('cm,abm->abc',H2A%uooo(:,:,k,j),r2a(:,:,i,:),x1)
                                         call einsum('abe,ec->abc',H2A%uuou(:,:,i,:),r2a(:,:,j,k),x2)
                                         xd7 = xd7 + (x2 - x1)                                 
                                         call einsum('cm,abm->abc',H2A%uooo(:,:,i,j),r2a(:,:,k,:),x1)
                                         call einsum('abe,ec->abc',H2A%uuou(:,:,k,:),r2a(:,:,j,i),x2)
                                         xd7 = xd7 - (x2 - x1)
                                         call einsum('cm,abm->abc',H2A%uooo(:,:,k,i),r2a(:,:,j,:),x1)
                                         call einsum('abe,ec->abc',H2A%uuou(:,:,j,:),r2a(:,:,i,k),x2)
                                         xd7 = xd7 - (x2 - x1)

                                         ! Diagram 9
                                         xd9 = 0.0
                                         call einsum('mc,abm->abc',I3_d9(:,:,j,k),t2a(:,:,i,:),x1)
                                         xd9 = xd9 - x1
                                         call einsum('mc,abm->abc',I3_d9(:,:,i,k),t2a(:,:,j,:),x1)
                                         xd9 = xd9 + x1
                                         call einsum('mc,abm->abc',I3_d9(:,:,j,i),t2a(:,:,k,:),x1)
                                         xd9 = xd9 + x1

                                         do a = 1,nua
                                            do b = a+1,nua
                                               do c = b+1,nua

                                                  z =      l1a(a,i)*vA%oouu(j,k,b,c) + l2a(b,c,j,k)*H1A%ou(i,a)  ! (1)
                                                  z = z - (l1a(b,i)*vA%oouu(j,k,a,c) + l2a(a,c,j,k)*H1A%ou(i,b)) ! (ab)
                                                  z = z - (l1a(c,i)*vA%oouu(j,k,b,a) + l2a(b,a,j,k)*H1A%ou(i,c)) ! (ac)
                                                  z = z - (l1a(a,j)*vA%oouu(i,k,b,c) + l2a(b,c,i,k)*H1A%ou(j,a)) ! (ij)
                                                  z = z - (l1a(a,k)*vA%oouu(j,i,b,c) + l2a(b,c,j,i)*H1A%ou(k,a)) ! (ik)
                                                  z = z + (l1a(b,j)*vA%oouu(i,k,a,c) + l2a(a,c,i,k)*H1A%ou(j,b)) ! (ab)(ij)
                                                  z = z + (l1a(b,k)*vA%oouu(j,i,a,c) + l2a(a,c,j,i)*H1A%ou(k,b)) ! (ab)(ik)
                                                  z = z + (l1a(c,j)*vA%oouu(i,k,b,a) + l2a(b,a,i,k)*H1A%ou(j,c)) ! (ac)(ij)
                                                  z = z + (l1a(c,k)*vA%oouu(j,i,b,a) + l2a(b,a,j,i)*H1A%ou(k,c)) ! (ac)(ik)
                                                  
                                                  L3 = y3(a,b,c) - y3(c,b,a) - y3(a,c,b) + z

                                                  MM23 =        xd1(a,b,c) - xd1(c,b,a) - xd1(a,c,b)
                                                  MM23 = MM23 + xd2(a,b,c) - xd2(b,a,c) - xd2(c,b,a)
                                                  MM23 = MM23 + xd3(a,b,c) - xd3(b,a,c) - xd3(c,b,a)
                                                  MM23 = MM23 + xd4(a,b,c) - xd4(a,c,b) - xd4(c,b,a)
                                                  MM23 = MM23 + xd5(a,b,c) - xd5(c,b,a) - xd5(a,c,b)
                                                  MM23 = MM23 + xd7(a,b,c) - xd7(c,b,a) - xd7(a,c,b)
                                                  MM23 = MM23 + xd9(a,b,c) - xd9(c,b,a) - xd9(a,c,b)

                                                  val = MM23*L3

                                                  denom = fA%oo(i,i)+fA%oo(j,j)+fA%oo(k,k)&
                                                         -fA%uu(a,a)-fA%uu(b,b)-fA%uu(c,c)

                                                  E23A = E23A + 2.0 * val/(omega+denom)

                                                  denom = H1A%oo(i,i)+H1A%oo(j,j)+H1A%oo(k,k)&
                                                         -H1A%uu(a,a)-H1A%uu(b,b)-H1A%uu(c,c)

                                                  E23B = E23B + 2.0 * val/(omega+denom)

                                                  denom = denom &
                                                        -H2A%uoou(a,i,i,a) - H2A%uoou(b,i,i,b) - H2A%uoou(c,i,i,c)&
                                                        -H2A%uoou(a,j,j,a) - H2A%uoou(b,j,j,b) - H2A%uoou(c,j,j,c)&
                                                        -H2A%uoou(a,k,k,a) - H2A%uoou(b,k,k,b) - H2A%uoou(c,k,k,c)&
                                                        -H2A%oooo(j,i,j,i) - H2A%oooo(k,i,k,i) - H2A%oooo(k,j,k,j)&
                                                        -H2A%uuuu(b,a,b,a) - H2A%uuuu(c,a,c,a) - H2A%uuuu(c,b,c,b)

                                                  E23C = E23C + 2.0 * val/(omega+denom)

                                                  denom = denom &
                                                        +D3A_O(a,i,j)+D3A_O(a,i,k)+D3A_O(a,j,k)&
                                                        +D3A_O(b,i,j)+D3A_O(b,i,k)+D3A_O(b,j,k)&
                                                        +D3A_O(c,i,j)+D3A_O(c,i,k)+D3A_O(c,j,k)&
                                                        -D3A_V(a,i,b)-D3A_V(a,i,c)-D3A_V(b,i,c)&
                                                        -D3A_V(a,j,b)-D3A_V(a,j,c)-D3A_V(b,j,c)&
                                                        -D3A_V(a,k,b)-D3A_V(a,k,c)-D3A_V(b,k,c)

                                                  E23D  = E23D + 2.0 * val/(omega+denom)


                                               end do
                                            end do
                                         end do
                                        
                                      end do
                                   end do
                                end do

                                deallocate(x1,x2,y1,y2,y3,xd1,xd2,xd3,xd4,xd5,xd7,xd9)
                                deallocate(I1_d1,I2_d1,I1_d2,I2_d2,I1_d3,I1_d4,I1_d5,I2_d5,I1_d6,I2_d6,I3_d9)

                                !print*,E23A
                                !print*,E23B
                                !print*,E23C
                                !print*,E23D

                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!!!! MM(2,3)B CORRECTION !!!!!
                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                allocate(I1_d1(noa,nub,noa,nob),I2_d1(nua,nub,nua,nob),I3_d1(nua,nob,noa,nob))
                                call einsum('mcie,ek->mcik',H2B%ouou,r1b,I1_d1)
                                call einsum('acfe,ek->acfk',H2B%uuuu,r1b,I2_d1)
                                call einsum('amie,ek->amik',H2B%uoou,r1b,I3_d1)

                                allocate(xd1(nua,nua,nua),xd2(nua,nua,nub),xd3(nua,nua,nub),&
                                         xd4(nua,nua,nub),xd5(nua,nua,nub),&
                                         xd7(nua,nua,nub),xd9(nua,nua,nub),&
                                         x1(nua,nua,nub), x2(nua,nua,nub),x3(nua,nua,nub),x4(nua,nua,nub),&
                                         y1(nua,nua,nub),y2(nua,nua,nub),y3(nua,nua,nub),y4(nua,nua,nub))

                                do i = 1,noa
                                   do j = i+1,noa
                                      do k = 1,nob

                                         ! Diagram 1
                                         xd1 = 0.0
                                         call einsum('mc,abm->abc',I1_d1(:,:,i,k),t2a(:,:,:,j),x1)
                                         call einsum('mc,abm->abc',I1_d1(:,:,j,k),t2a(:,:,:,i),x2)
                                         xd1 = xd1 - x1 + x2
                                         call einsum('acf,fb->abc',I2_d1(:,:,:,k),t2a(:,:,i,j),x1)
                                         xd1 = xd1 + x1
                                         call einsum('am,bcm->abc',I3_d1(:,:,i,k),t2b(:,:,j,:),x1)
                                         call einsum('am,bcm->abc',I3_d1(:,:,j,k),t2b(:,:,i,:),x2)
                                         xd1 = xd1 - x1 + x2


                                         do a = 1,nua
                                            do b = a+1,nua
                                               do c = 1,nub

                                                  ! Approximate L3B

                                                  ! Diagram 1
                                                  MM23 = xd1(a,b,c)

                                               end do
                                            end do
                                         end do

                                      end do
                                   end do
                                end do


                                                  






                        end associate

                        !deallocate(I1_d1,I2_d1,I1_d2,I2_d2,I1_d3,xd1,xd2,xd3,xd4,xd5,xd7,xd9,y1,y2,y3)


                        deallocate(l1a,l1b,l2a,l2b,l2c,r1a,r1b,r2a,r2b,r2c)
                
                        call calc_cc_energy(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Ecorr)

                        write(io,fmt=*) 'EOMCCSD = ',omega+Ecorr+sys%Escf,'omega = ',omega 
                        write(io,fmt=*) 'CR-EOMCC(2,3)A = ',E23A+omega+Ecorr+sys%Escf,'EcorrA = ',E23A+omega,'deltaA = ',E23A
                        write(io,fmt=*) 'CR-EOMCC(2,3)B = ',E23B+omega+Ecorr+sys%Escf,'EcorrB = ',E23B+omega,'deltaB = ',E23B
                        write(io,fmt=*) 'CR-EOMCC(2,3)C = ',E23C+omega+Ecorr+sys%Escf,'EcorrC = ',E23C+omega,'deltaC = ',E23C
                        write(io,fmt=*) 'CR-EOMCC(2,3)D = ',E23D+omega+Ecorr+sys%Escf,'EcorrD = ',E23D+omega,'deltaD = ',E23D

                end subroutine creomcc23

                subroutine triples_3body_diagonal(sys,vA,vB,vC,t2a,t2b,t2c,&
                                D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O)

                        type(sys_t), intent(in) :: sys
                        type(e2int_t), intent(in) :: vA, vB, vC
                        real, intent(in) :: t2a(sys%Nunocc_a,sys%Nunocc_a,sys%Nocc_a,sys%Nocc_a),&
                                            t2b(sys%Nunocc_a,sys%Nunocc_b,sys%Nocc_a,sys%Nocc_b),&
                                            t2c(sys%Nunocc_b,sys%Nunocc_b,sys%Nocc_b,sys%Nocc_b)
                        real, allocatable, intent(out) :: D3A_V(:,:,:), D3A_O(:,:,:), D3B_V(:,:,:), D3B_O(:,:,:), &
                                                          D3C_V(:,:,:), D3C_O(:,:,:), D3D_V(:,:,:), D3D_O(:,:,:)
                        integer :: a, b, c, i, j, k

                        associate(noa=>sys%Nocc_a,nua=>sys%Nunocc_a,nob=>sys%Nocc_b,nub=>sys%Nunocc_b)

                                allocate(D3A_V(nua,noa,nua),D3A_O(nua,noa,noa),D3B_V(nua,noa,nub),D3B_O(nua,noa,nob),&
                                         D3C_V(nua,nob,nub),D3C_O(nub,noa,nob),D3D_V(nub,nob,nub),D3D_O(nub,nob,nob))

                                do b = 1,nua
                                   do i = 1,noa
                                      do a = 1,nua
                                         D3A_V(a,i,b) = -1.0*dot(vA%oouu(i,:,a,b),t2a(a,b,i,:))
                                      end do
                                   end do
                                end do

                                do j = 1,noa
                                   do i = 1,noa
                                      do a = 1,nua
                                         D3A_O(a,i,j) = dot(vA%oouu(i,j,a,:),t2a(a,:,i,j))
                                      end do
                                   end do
                                end do

                                do c = 1,nub
                                   do i = 1,noa
                                      do a = 1,nua
                                         D3B_V(a,i,c) = -1.0*dot(vB%oouu(i,:,a,c),t2b(a,c,i,:))
                                      end do
                                   end do
                                end do

                                do k = 1,nob
                                   do i = 1,noa
                                      do a = 1,nua
                                         D3B_O(a,i,k) = dot(vB%oouu(i,k,a,:),t2b(a,:,i,k))
                                      end do
                                   end do
                                end do

                                do c = 1,nub
                                   do k = 1,nob
                                      do a = 1,nua
                                         D3C_V(a,k,c) = -1.0*dot(vB%oouu(:,k,a,c),t2b(a,c,:,k))
                                      end do
                                   end do
                                end do

                                do k = 1,nob
                                   do i = 1,noa
                                      do c = 1,nub
                                         D3C_O(c,i,k) = dot(vB%oouu(i,k,:,c),t2b(:,c,i,k))
                                      end do
                                   end do
                                end do

                                 
                                do b = 1,nub
                                   do i = 1,nob
                                      do a = 1,nub
                                         D3D_V(a,i,b) = -1.0*dot(vC%oouu(i,:,a,b),t2c(a,b,i,:))
                                      end do
                                   end do
                                end do

                                do j = 1,nob
                                   do i = 1,nob
                                      do a = 1,nub
                                         D3D_O(a,i,j) = dot(vC%oouu(i,j,a,:),t2c(a,:,i,j))
                                      end do
                                   end do
                                end do

                        end associate

                end subroutine triples_3body_diagonal

                subroutine enforce_LR(l1a,l1b,l2a,l2b,l2c,r1a,r1b,r2a,r2b,r2c)

                        real, intent(inout) :: l1a(:,:), l1b(:,:), &
                                               l2a(:,:,:,:), l2b(:,:,:,:), l2c(:,:,:,:)
                        real, intent(in) :: r1a(:,:), r1b(:,:), &
                                            r2a(:,:,:,:), r2b(:,:,:,:), r2c(:,:,:,:)
                        real :: LR
                        integer :: a, b, i, j, nua, noa, nub, nob

                        nua = size(l1a,1)
                        noa = size(l1a,2)
                        nub = size(l1b,1)
                        nob = size(l1b,2)

                        LR = 0.0

                        do i = 1,noa
                           do a = 1,nua
                              LR = LR + l1a(a,i)*r1a(a,i)
                           end do
                        end do

                        do i = 1,nob
                           do a = 1,nub
                              LR = LR + l1b(a,i)*r1b(a,i)
                           end do
                        end do

                        do i = 1,noa
                           do j = 1,noa
                              do a = 1,nua
                                 do b = 1,nua
                                    LR = LR + 0.25*l2a(a,b,i,j)*r2a(a,b,i,j)
                                 end do
                              end do
                           end do
                        end do

                        do i = 1,noa
                           do j = 1,nob
                              do a = 1,nua
                                 do b = 1,nub
                                    LR = LR + l2b(a,b,i,j)*r2b(a,b,i,j)
                                 end do
                              end do
                           end do
                        end do

                        do i = 1,nob
                           do j = 1,nob
                              do a = 1,nub
                                 do b = 1,nub
                                    LR = LR + 0.25*l2c(a,b,i,j)*r2c(a,b,i,j)
                                 end do
                              end do
                           end do
                        end do

                        l1a = l1a/LR
                        l1b = l1b/LR
                        l2a = l2a/LR
                        l2b = l2b/LR
                        l2c = l2c/LR

                end subroutine enforce_LR

end module creomcc23_module

