module crcc23_module

        use dgemm_module, only: dot
        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use permutils, only: reorder_stripe
        use einsum_module, only: einsum

        implicit none

        contains

                subroutine crcc23(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,L,H1A,H1B,H2A,H2B,H2C,E23A,E23B,E23C,E23D,io)

                        use cc_energy, only: calc_cc_energy

                        integer, intent(in) :: io
                        type(sys_t), intent(in) :: sys
                        type(e1int_t), intent(in) :: fA, fB, H1A, H1B
                        type(e2int_t), intent(in) :: vA, vB, vC, H2A, H2B, H2C
                        real, intent(in) :: t1a(:,:), t1b(:,:), t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:), L(:)
                        real, intent(out) :: E23A, E23B, E23C, E23D
                        real, allocatable :: D3A_V(:,:,:), D3A_O(:,:,:), D3B_V(:,:,:), D3B_O(:,:,:), &
                                             D3C_V(:,:,:), D3C_O(:,:,:), D3D_V(:,:,:), D3D_O(:,:,:), &
                                             l1a(:,:), l1b(:,:), l2a(:,:,:,:), l2b(:,:,:,:), l2c(:,:,:,:)
                        real, allocatable :: I1(:,:,:,:), I2(:,:,:,:), I3(:,:,:,:), &
                                             x1(:,:,:), x2(:,:,:), x3(:,:,:), x4(:,:,:),&
                                             y1(:,:,:), y2(:,:,:), y3(:,:,:), y4(:,:,:)
                        integer :: a, b, c, i, j, k, n1a, n1b, n2a, n2b, n2c
                        real :: MM23, L3, denom, val, z, Ecorr

                        write(io,fmt=*) ''
                        write(io,fmt=*) '++++++++++++++++++++ CR-CC(2,3) ROUTINE +++++++++++++++++++++'
                        write(io,fmt=*) ''
                        
                        ! initialize
                        E23A = 0.0
                        E23B = 0.0
                        E23C = 0.0
                        E23D = 0.0

                        associate(noa=>sys%Nocc_a,nob=>sys%Nocc_b,nua=>sys%Nunocc_a,nub=>sys%Nunocc_b)

                                n1a = nua*noa
                                n1b = nub*nob
                                n2a = nua**2 * noa**2
                                n2b = nua*nua*noa*nob
                                n2c = nub**2 * nob**2

                                allocate(l1a(nua,noa),l1b(nub,nob),l2a(nua,nua,noa,noa),l2b(nua,nub,noa,nob),l2c(nub,nub,nob,nob))

                                l1a = reshape(L(1:n1a),(/nua,noa/))
                                l1b = reshape(L(n1a+1:n1a+n1b),(/nub,nob/))
                                l2a = reshape(L(n1a+n1b+1:n1a+n1b+n2a),(/nua,nua,noa,noa/))
                                l2b = reshape(L(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b),(/nua,nub,noa,nob/))
                                l2c = reshape(L(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c),(/nub,nub,nob,nob/))

                                ! get 3-body HBar triples diagonal
                                call triples_3body_diagonal(sys,vA,vB,vC,t2a,t2b,t2c,&
                                        D3A_V,D3A_O,D3B_V,D3B_O,D3C_V,D3C_O,D3D_V,D3D_O)

                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!!!! MM(2,3)A CORRECTION !!!!!
                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                allocate(I1(nua,nua,noa,nua),&
                                         x1(nua,nua,nua),x2(nua,nua,nua),x3(nua,nua,nua),&
                                         y1(nua,nua,nua),y2(nua,nua,nua),y3(nua,nua,nua))

                                call einsum('me,abim->abie',H1A%ou,t2a,I1)
                                I1 = H2A%uuou + I1

                                do i = 1,noa
                                   do j = i+1,noa
                                      do k = j+1,noa

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

                                         x3 = 0.0
                                         call einsum('cm,abm->abc',H2A%uooo(:,:,k,j),t2a(:,:,i,:),x1)
                                         call einsum('abe,ec->abc',I1(:,:,i,:),t2a(:,:,j,k),x2)
                                         x3 = x3 + (x2 - x1)                                 
                                         call einsum('cm,abm->abc',H2A%uooo(:,:,i,j),t2a(:,:,k,:),x1)
                                         call einsum('abe,ec->abc',I1(:,:,k,:),t2a(:,:,j,i),x2)
                                         x3 = x3 - (x2 - x1)
                                         call einsum('cm,abm->abc',H2A%uooo(:,:,k,i),t2a(:,:,j,:),x1)
                                         call einsum('abe,ec->abc',I1(:,:,j,:),t2a(:,:,i,k),x2)
                                         x3 = x3 - (x2 - x1)

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
                                                  
                                                  MM23 = x3(a,b,c) - x3(c,b,a) - x3(a,c,b)
                                                  L3 = y3(a,b,c) - y3(c,b,a) - y3(a,c,b) + z
                                                  val = MM23*L3

                                                  denom = fA%oo(i,i)+fA%oo(j,j)+fA%oo(k,k)&
                                                         -fA%uu(a,a)-fA%uu(b,b)-fA%uu(c,c)

                                                  E23A = E23A + val/denom

                                                  denom = H1A%oo(i,i)+H1A%oo(j,j)+H1A%oo(k,k)&
                                                         -H1A%uu(a,a)-H1A%uu(b,b)-H1A%uu(c,c)

                                                  E23B = E23B + val/denom

                                                  denom = denom &
                                                        -H2A%uoou(a,i,i,a) - H2A%uoou(b,i,i,b) - H2A%uoou(c,i,i,c)&
                                                        -H2A%uoou(a,j,j,a) - H2A%uoou(b,j,j,b) - H2A%uoou(c,j,j,c)&
                                                        -H2A%uoou(a,k,k,a) - H2A%uoou(b,k,k,b) - H2A%uoou(c,k,k,c)&
                                                        -H2A%oooo(j,i,j,i) - H2A%oooo(k,i,k,i) - H2A%oooo(k,j,k,j)&
                                                        -H2A%uuuu(b,a,b,a) - H2A%uuuu(c,a,c,a) - H2A%uuuu(c,b,c,b)

                                                  E23C = E23C + val/denom

                                                  denom = denom &
                                                        +D3A_O(a,i,j)+D3A_O(a,i,k)+D3A_O(a,j,k)&
                                                        +D3A_O(b,i,j)+D3A_O(b,i,k)+D3A_O(b,j,k)&
                                                        +D3A_O(c,i,j)+D3A_O(c,i,k)+D3A_O(c,j,k)&
                                                        -D3A_V(a,i,b)-D3A_V(a,i,c)-D3A_V(b,i,c)&
                                                        -D3A_V(a,j,b)-D3A_V(a,j,c)-D3A_V(b,j,c)&
                                                        -D3A_V(a,k,b)-D3A_V(a,k,c)-D3A_V(b,k,c)

                                                  E23D  = E23D + val/denom

                                               end do
                                            end do
                                         end do

                                      end do
                                   end do
                                end do

                                deallocate(x1,x2,x3,y1,y2,y3)


                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!!!! MM(2,3)B CORRECTION !!!!!
                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                allocate(I2(nua,nub,nua,nob),I3(nua,nub,noa,nub),&
                                         x1(nua,nua,nub),x2(nua,nua,nub),x3(nua,nua,nub),x4(nua,nua,nub),&
                                         y1(nua,nua,nub),y2(nua,nua,nub),y3(nua,nua,nub),y4(nua,nua,nub))

                                call einsum('me,bcmk->bcek',H1A%ou,t2b,I2)
                                I2 = H2B%uuuo + I2
                                call einsum('me,bcjm->bcje',H1B%ou,t2b,I3)
                                I3 = H2B%uuou + I3

                                do i = 1,noa
                                   do j = i+1,noa
                                      do k = 1,nob

                                         y4 = 0.0
                                         call einsum('eba,ec->abc',H2A%uouu(:,i,:,:),l2b(:,:,j,k),y1)
                                         call einsum('mc,abm->abc',H2B%ooou(j,k,:,:),l2a(:,:,i,:),y2)
                                         y4 = y4 + (y1 - y2)
                                         call einsum('eba,ec->abc',H2A%uouu(:,j,:,:),l2b(:,:,i,k),y1)
                                         call einsum('mc,abm->abc',H2B%ooou(i,k,:,:),l2a(:,:,j,:),y2)
                                         y4 = y4 - (y1 - y2)

                                         y3 = 0.0
                                         call einsum('ebc,ae->abc',H2B%uouu(:,k,:,:),l2a(:,:,i,j),y1)
                                         call einsum('ma,bcm->abc',H2A%ooou(j,i,:,:),l2b(:,:,:,k),y2)
                                         y3 = y3 + (y1 - y2)

                                         call einsum('ebc,ae->abc',H2B%ouuu(j,:,:,:),l2b(:,:,i,k),y1)
                                         call einsum('bm,acm->abc',H2B%oouo(j,k,:,:),l2b(:,:,i,:),y2)
                                         y3 = y3 + (y1 - y2)
                                         call einsum('ebc,ae->abc',H2B%ouuu(i,:,:,:),l2b(:,:,j,k),y1)
                                         call einsum('bm,acm->abc',H2B%oouo(i,k,:,:),l2b(:,:,j,:),y2)
                                         y3 = y3 - (y1 - y2)                                

                                         x4 = 0.0
                                         call einsum('abe,ec->abc',I1(:,:,i,:),t2b(:,:,j,k),x1)
                                         call einsum('mc,abm->abc',H2B%ouoo(:,:,j,k),t2a(:,:,i,:),x2)                          
                                         x4 = x4 + (x1 - x2)
                                         call einsum('abe,ec->abc',I1(:,:,j,:),t2b(:,:,i,k),x1)
                                         call einsum('mc,abm->abc',H2B%ouoo(:,:,i,k),t2a(:,:,j,:),x2)                          
                                         x4 = x4 - (x1 - x2)

                                         x3 = 0.0
                                         call einsum('bce,ae->abc',I2(:,:,:,k),t2a(:,:,i,j),x1)
                                         call einsum('am,bcm->abc',H2A%uooo(:,:,i,j),t2b(:,:,:,k),x2)
                                         x3 = x3 + (x1 - x2)
                                         call einsum('bce,ae->abc',I3(:,:,j,:),t2b(:,:,i,k),x1)
                                         call einsum('bm,acm->abc',H2B%uooo(:,:,j,k),t2b(:,:,i,:),x2)
                                         x3 = x3 + (x1 - x2)
                                         call einsum('bce,ae->abc',I3(:,:,i,:),t2b(:,:,j,k),x1)
                                         call einsum('bm,acm->abc',H2B%uooo(:,:,i,k),t2b(:,:,j,:),x2) 
                                         x3 = x3 - (x1 - x2)


                                         do a = 1,nua
                                            do b = a+1,nua
                                               do c = 1,nub

                                                  z = l1b(c,k)*vA%oouu(i,j,a,b)
                                                  z = z + H1B%ou(k,c)*l2a(a,b,i,j)
                                                  z = z + l1a(b,j)*vB%oouu(i,k,a,c)
                                                  z = z - l1a(a,j)*vB%oouu(i,k,b,c)
                                                  z = z - l1a(b,i)*vB%oouu(j,k,a,c)
                                                  z = z + l1a(a,i)*vB%oouu(j,k,b,c)
                                                  z = z + H1A%ou(j,b)*l2b(a,c,i,k)
                                                  z = z - H1A%ou(j,a)*l2b(b,c,i,k)
                                                  z = z - H1A%ou(i,b)*l2b(a,c,j,k)
                                                  z = z + H1A%ou(i,a)*l2b(b,c,j,k)

                                                  ! calculate A(c/ab) x3(abc)
                                                  MM23 = x3(a,b,c) - x3(b,a,c) + x4(a,b,c)
                                                  L3 = y3(a,b,c) - y3(b,a,c) + y4(a,b,c) + z
                                                  val = MM23*L3

                                                  denom = fA%oo(i,i)+fA%oo(j,j)+fB%oo(k,k)&
                                                          -fA%uu(a,a)-fA%uu(b,b)-fB%uu(c,c)

                                                  E23A = E23A + val/denom

                                                  denom = H1A%oo(i,i)+H1A%oo(j,j)+H1B%oo(k,k)&
                                                          -H1A%uu(a,a)-H1A%uu(b,b)-H1B%uu(c,c)
                                                  
                                                  E23B = E23B + val/denom

                                                  denom = denom &
                                                         -H2A%uoou(a,i,i,a)-H2A%uoou(b,i,i,b)+H2B%ouou(i,c,i,c)&
                                                         -H2A%uoou(a,j,j,a)-H2A%uoou(b,j,j,b)+H2B%ouou(j,c,j,c)&
                                                         +H2B%uouo(a,k,a,k)+H2B%uouo(b,k,b,k)-H2C%uoou(c,k,k,c)&
                                                         -H2A%oooo(j,i,j,i)-H2B%oooo(i,k,i,k)-H2B%oooo(j,k,j,k)&
                                                         -H2A%uuuu(b,a,b,a)-H2B%uuuu(a,c,a,c)-H2B%uuuu(b,c,b,c)
            
                                                  E23C = E23C + val/denom 

                                                  denom = denom &
                                                        +D3A_O(a,i,j)+D3B_O(a,i,k)+D3B_O(a,j,k)&
                                                        +D3A_O(b,i,j)+D3B_O(b,i,k)+D3B_O(b,j,k)&
                                                        +D3C_O(c,i,k)+D3C_O(c,j,k)&
                                                        -D3A_V(a,i,b)-D3B_V(a,i,c)-D3B_V(b,i,c)&
                                                        -D3A_V(a,j,b)-D3B_V(a,j,c)-D3B_V(b,j,c)&
                                                        -D3C_V(a,k,c)-D3C_V(b,k,c)

                                                  E23D = E23D + val/denom

                                               end do
                                            end do
                                         end do

                                      end do
                                   end do
                                end do
                               
                                deallocate(I1,I2,I3,x1,x2,x3,x4,y1,y2,y3,y4)


                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!!!! MM(2,3)D CORRECTION !!!!!
                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                allocate(I1(nub,nub,nob,nub),&
                                         x1(nub,nub,nub),x2(nub,nub,nub),x3(nub,nub,nub),&
                                         y1(nub,nub,nub),y2(nub,nub,nub),y3(nub,nub,nub))

                                call einsum('me,abim->abie',H1B%ou,t2c,I1)
                                I1 = H2C%uuou + I1

                                do i = 1,nob
                                   do j = i+1,nob
                                      do k = j+1,nob

                                         y3 = 0.0
                                         call einsum('eba,ec->abc',H2C%uouu(:,i,:,:),l2c(:,:,j,k),y1)
                                         call einsum('mc,abm->abc',H2C%ooou(j,k,:,:),l2c(:,:,i,:),y2)
                                         y3 = y3 + (y1 - y2)                                 
                                         call einsum('eba,ec->abc',H2C%uouu(:,j,:,:),l2c(:,:,i,k),y1)
                                         call einsum('mc,abm->abc',H2C%ooou(i,k,:,:),l2c(:,:,j,:),y2)
                                         y3 = y3 - (y1 - y2)
                                         call einsum('eba,ec->abc',H2C%uouu(:,k,:,:),l2c(:,:,j,i),y1)
                                         call einsum('mc,abm->abc',H2C%ooou(j,i,:,:),l2c(:,:,k,:),y2)
                                         y3 = y3 - (y1 - y2)

                                         x3 = 0.0
                                         call einsum('cm,abm->abc',H2C%uooo(:,:,k,j),t2c(:,:,i,:),x1)
                                         call einsum('abe,ec->abc',I1(:,:,i,:),t2c(:,:,j,k),x2)
                                         x3 = x3 + (x2 - x1)                                 
                                         call einsum('cm,abm->abc',H2C%uooo(:,:,i,j),t2c(:,:,k,:),x1)
                                         call einsum('abe,ec->abc',I1(:,:,k,:),t2c(:,:,j,i),x2)
                                         x3 = x3 - (x2 - x1)
                                         call einsum('cm,abm->abc',H2C%uooo(:,:,k,i),t2c(:,:,j,:),x1)
                                         call einsum('abe,ec->abc',I1(:,:,j,:),t2c(:,:,i,k),x2)
                                         x3 = x3 - (x2 - x1)

                                         do a = 1,nub
                                            do b = a+1,nub
                                               do c = b+1,nub

                                                  z =      l1b(a,i)*vC%oouu(j,k,b,c) + l2c(b,c,j,k)*H1B%ou(i,a)  ! (1)
                                                  z = z - (l1b(b,i)*vC%oouu(j,k,a,c) + l2c(a,c,j,k)*H1B%ou(i,b)) ! (ab)
                                                  z = z - (l1b(c,i)*vC%oouu(j,k,b,a) + l2c(b,a,j,k)*H1B%ou(i,c)) ! (ac)
                                                  z = z - (l1b(a,j)*vC%oouu(i,k,b,c) + l2c(b,c,i,k)*H1B%ou(j,a)) ! (ij)
                                                  z = z - (l1b(a,k)*vC%oouu(j,i,b,c) + l2c(b,c,j,i)*H1B%ou(k,a)) ! (ik)
                                                  z = z + (l1b(b,j)*vC%oouu(i,k,a,c) + l2c(a,c,i,k)*H1B%ou(j,b)) ! (ab)(ij)
                                                  z = z + (l1b(b,k)*vC%oouu(j,i,a,c) + l2c(a,c,j,i)*H1B%ou(k,b)) ! (ab)(ik)
                                                  z = z + (l1b(c,j)*vC%oouu(i,k,b,a) + l2c(b,a,i,k)*H1B%ou(j,c)) ! (ac)(ij)
                                                  z = z + (l1b(c,k)*vC%oouu(j,i,b,a) + l2c(b,a,j,i)*H1B%ou(k,c)) ! (ac)(ik)
                                                  
                                                  MM23 = x3(a,b,c) - x3(c,b,a) - x3(a,c,b)
                                                  L3 = y3(a,b,c) - y3(c,b,a) - y3(a,c,b) + z
                                                  val = MM23*L3

                                                  denom = fB%oo(i,i)+fB%oo(j,j)+fB%oo(k,k)&
                                                         -fB%uu(a,a)-fB%uu(b,b)-fB%uu(c,c)

                                                  E23A = E23A + val/denom

                                                  denom = H1B%oo(i,i)+H1B%oo(j,j)+H1B%oo(k,k)&
                                                         -H1B%uu(a,a)-H1B%uu(b,b)-H1B%uu(c,c)

                                                  E23B = E23B + val/denom

                                                  denom = denom &
                                                        -H2C%uoou(a,i,i,a) - H2C%uoou(b,i,i,b) - H2C%uoou(c,i,i,c)&
                                                        -H2C%uoou(a,j,j,a) - H2C%uoou(b,j,j,b) - H2C%uoou(c,j,j,c)&
                                                        -H2C%uoou(a,k,k,a) - H2C%uoou(b,k,k,b) - H2C%uoou(c,k,k,c)&
                                                        -H2C%oooo(j,i,j,i) - H2C%oooo(k,i,k,i) - H2C%oooo(k,j,k,j)&
                                                        -H2C%uuuu(b,a,b,a) - H2C%uuuu(c,a,c,a) - H2C%uuuu(c,b,c,b)

                                                  E23C = E23C + val/denom

                                                  denom = denom &
                                                        +D3D_O(a,i,j)+D3D_O(a,i,k)+D3D_O(a,j,k)&
                                                        +D3D_O(b,i,j)+D3D_O(b,i,k)+D3D_O(b,j,k)&
                                                        +D3D_O(c,i,j)+D3D_O(c,i,k)+D3D_O(c,j,k)&
                                                        -D3D_V(a,i,b)-D3D_V(a,i,c)-D3D_V(b,i,c)&
                                                        -D3D_V(a,j,b)-D3D_V(a,j,c)-D3D_V(b,j,c)&
                                                        -D3D_V(a,k,b)-D3D_V(a,k,c)-D3D_V(b,k,c)

                                                  E23D  = E23D + val/denom

                                               end do
                                            end do
                                         end do

                                      end do
                                   end do
                                end do

                                deallocate(x1,x2,x3,y1,y2,y3)

                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                !!!!! MM(2,3)C CORRECTION !!!!!
                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                allocate(I2(nua,nub,noa,nub),I3(nua,nub,nua,nob),&
                                         x1(nua,nub,nub),x2(nua,nub,nub),x3(nua,nub,nub),x4(nua,nub,nub),&
                                         y1(nua,nub,nub),y2(nua,nub,nub),y3(nua,nub,nub),y4(nua,nub,nub))

                                call einsum('me,abim->abie',H1B%ou,t2b,I2)
                                I2 = H2B%uuou + I2
                                call einsum('me,abmj->abej',H1A%ou,t2b,I3)
                                I3 = H2B%uuuo + I3

                                do i = 1,noa
                                   do j = 1,nob
                                      do k = j+1,nob

                                         y4 = 0.0
                                         call einsum('ebc,ae->abc',H2C%uouu(:,k,:,:),l2b(:,:,i,j),y1)
                                         call einsum('am,bcm->abc',H2B%oouo(i,j,:,:),l2c(:,:,:,k),y2)
                                         y4 = y4 + (y1 - y2)                                         
                                         call einsum('ebc,ae->abc',H2C%uouu(:,j,:,:),l2b(:,:,i,k),y1)
                                         call einsum('am,bcm->abc',H2B%oouo(i,k,:,:),l2c(:,:,:,j),y2)
                                         y4 = y4 - (y1 - y2)

                                         y3 = 0.0
                                         call einsum('eab,ec->abc',H2B%ouuu(i,:,:,:),l2c(:,:,j,k),y1)
                                         call einsum('mc,abm->abc',H2C%ooou(j,k,:,:),l2b(:,:,i,:),y2)
                                         y3 = y3 + (y1 - y2)
                                         call einsum('eab,ec->abc',H2B%uouu(:,j,:,:),l2b(:,:,i,k),y1)
                                         call einsum('mb,acm->abc',H2B%ooou(i,j,:,:),l2b(:,:,:,k),y2)
                                         y3 = y3 + (y1 - y2)
                                         call einsum('eab,ec->abc',H2B%uouu(:,k,:,:),l2b(:,:,i,j),y1)
                                         call einsum('mb,acm->abc',H2B%ooou(i,k,:,:),l2b(:,:,:,j),y2)
                                         y3 = y3 - (y1 - y2)

                
                                         x4 = 0.0
                                         call einsum('bce,ae->abc',I1(:,:,j,:),t2b(:,:,i,k),x1)
                                         call einsum('am,bcm->abc',H2B%uooo(:,:,i,j),t2c(:,:,:,k),x2)
                                         x4 = x4 + (x1 - x2)
                                         call einsum('bce,ae->abc',I1(:,:,k,:),t2b(:,:,i,j),x1)
                                         call einsum('am,bcm->abc',H2B%uooo(:,:,i,k),t2c(:,:,:,j),x2)
                                         x4 = x4 - (x1 - x2)

                                         x3 = 0.0
                                         call einsum('abe,ec->abc',I2(:,:,i,:),t2c(:,:,j,k),x1)
                                         call einsum('cm,abm->abc',H2C%uooo(:,:,k,j),t2b(:,:,i,:),x2)
                                         x3 = x3 + (x1 - x2)
                                         call einsum('abe,ec->abc',I3(:,:,:,j),t2b(:,:,i,k),x1)
                                         call einsum('mb,acm->abc',H2B%ouoo(:,:,i,j),t2b(:,:,:,k),x2)
                                         x3 = x3 + (x1 - x2)
                                         call einsum('abe,ec->abc',I3(:,:,:,k),t2b(:,:,i,j),x1)
                                         call einsum('mb,acm->abc',H2B%ouoo(:,:,i,k),t2b(:,:,:,j),x2)
                                         x3 = x3 - (x1 - x2)

                                         do a = 1,nua
                                            do b = 1,nub
                                               do c = b+1,nub

                                                  z = l1a(a,i)*vC%oouu(j,k,b,c)
                                                  z = z + H1A%ou(i,a)*l2c(b,c,j,k)
                                                  z = z + l1b(b,j)*vB%oouu(i,k,a,c)
                                                  z = z - l1b(c,j)*vB%oouu(i,k,a,b)
                                                  z = z - l1b(b,k)*vB%oouu(i,j,a,c)
                                                  z = z + l1b(c,k)*vB%oouu(i,j,a,b)
                                                  z = z + H1B%ou(j,b)*l2b(a,c,i,k)
                                                  z = z - H1B%ou(k,b)*l2b(a,c,i,j)
                                                  z = z - H1B%ou(j,c)*l2b(a,b,i,k)
                                                  z = z + H1B%ou(k,c)*l2b(a,b,i,j)

                                                  MM23 = x3(a,b,c) - x3(a,c,b) + x4(a,b,c)
                                                  L3 = y3(a,b,c) - y3(a,c,b) + y4(a,b,c) + z
                                                  val = MM23*L3

                                                  denom = fA%oo(i,i)+fB%oo(j,j)+fB%oo(k,k)&
                                                          -fA%uu(a,a)-fB%uu(b,b)-fB%uu(c,c)

                                                  E23A = E23A + val/denom

                                                  denom = H1A%oo(i,i)+H1B%oo(j,j)+H1B%oo(k,k)&
                                                          -H1A%uu(a,a)-H1B%uu(b,b)-H1B%uu(c,c)
                                                  
                                                  E23B = E23B + val/denom

                                                  denom = denom &
                                                        -H2A%uoou(a,i,i,a)+H2B%ouou(i,b,i,b)+H2B%ouou(i,c,i,c)&
                                                        +H2B%uouo(a,j,a,j)-H2C%uoou(b,j,j,b)-H2C%uoou(c,j,j,c)&
                                                        +H2B%uouo(a,k,a,k)-H2C%uoou(b,k,k,b)-H2C%uoou(c,k,k,c)&
                                                        -H2B%oooo(i,j,i,j)-H2B%oooo(i,k,i,k)-H2C%oooo(k,j,k,j)&
                                                        -H2B%uuuu(a,b,a,b)-H2B%uuuu(a,c,a,c)-H2C%uuuu(c,b,c,b)
            
                                                  E23C = E23C + val/denom 

                                                  denom  = denom &
                                                        +D3B_O(a,i,j)+D3B_O(a,i,k)&
                                                        +D3C_O(b,i,j)+D3C_O(b,i,k)+D3D_O(b,j,k)&
                                                        +D3C_O(c,i,j)+D3C_O(c,i,k)+D3D_O(c,j,k)&
                                                        -D3B_V(a,i,b)-D3B_V(a,i,c)&
                                                        -D3C_V(a,j,b)-D3C_V(a,j,c)-D3D_V(b,j,c)&
                                                        -D3C_V(a,k,b)-D3C_V(a,k,c)-D3D_V(b,k,c)

                                                  E23D = E23D + val/denom

                                               end do
                                            end do
                                         end do

                                      end do
                                   end do
                                end do
                               
                                deallocate(I1,I2,I3,x1,x2,x3,x4,y1,y2,y3,y4)

                        end associate

                        deallocate(l1a,l1b,l2a,l2b,l2c)

                        call calc_cc_energy(sys,fA,fB,vA,vB,vC,t1a,t1b,t2a,t2b,t2c,Ecorr)

                        write(io,fmt=*) 'CR-CC(2,3)A = ',E23A+Ecorr+sys%Escf,'EcorrA = ',E23A+Ecorr,'deltaA = ',E23A
                        write(io,fmt=*) 'CR-CC(2,3)B = ',E23B+Ecorr+sys%Escf,'EcorrB = ',E23B+Ecorr,'deltaB = ',E23B
                        write(io,fmt=*) 'CR-CC(2,3)C = ',E23C+Ecorr+sys%Escf,'EcorrC = ',E23C+Ecorr,'deltaC = ',E23C
                        write(io,fmt=*) 'CR-CC(2,3)D = ',E23D+Ecorr+sys%Escf,'EcorrD = ',E23D+Ecorr,'deltaD = ',E23D

                end subroutine crcc23

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

end module crcc23_module
