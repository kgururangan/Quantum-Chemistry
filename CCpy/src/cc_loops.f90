module cc_loops

  implicit none
  
  contains

      subroutine update_t1a(t1a,X1A,fA_oo,fA_vv,shift,noa,nua)

              implicit none

              integer, intent(in) :: noa, nua
              real(8), intent(in) :: fA_oo(1:noa,1:noa), fA_vv(1:nua,1:nua), &
                                  X1A(1:nua,1:noa), shift               
              real(8), intent(inout) :: t1a(1:nua,1:noa)
              !f2py intent(in,out) :: t1a(0:nua-1,0:noa-1)
              integer :: i, a
              real(8) :: denom

              do i = 1,noa
                do a = 1,nua
                  denom = fA_oo(i,i) - fA_vv(a,a)
                  t1a(a,i) = t1a(a,i) + X1A(a,i)/(denom-shift)
                end do
              end do

      end subroutine update_t1a

      subroutine update_t1b(t1b,X1B,fB_oo,fB_vv,shift,nob,nub)

              implicit none

              integer, intent(in) :: nob, nub
              real(8), intent(in) :: fB_oo(1:nob,1:nob), fB_vv(1:nub,1:nub), &
                                  X1B(1:nub,1:nob), shift               
              real(8), intent(inout) :: t1b(1:nub,1:nob)
              !f2py intent(in,out) :: t1b(0:nub-1,0:nob-1)
              integer :: i, a
              real(8) :: denom

              do i = 1,nob
                do a = 1,nub
                  denom = fB_oo(i,i) - fB_vv(a,a)
                  t1b(a,i) = t1b(a,i) + X1B(a,i)/(denom-shift)
                end do
              end do

      end subroutine update_t1b

      subroutine update_t2a(t2a,X2A,fA_oo,fA_vv,shift,noa,nua)

              implicit none

              integer, intent(in) :: noa, nua
              real(8), intent(in) :: fA_oo(1:noa,1:noa), fA_vv(1:nua,1:nua), &
                                  X2A(1:nua,1:nua,1:noa,1:noa), shift               
              real(8), intent(inout) :: t2a(1:nua,1:nua,1:noa,1:noa)
              !f2py intent(in,out) :: t2a(0:nua-1,0:nua-1,0:noa-1,0:noa-1)
              integer :: i, j, a, b
              real(8) :: denom

              do i = 1,noa
                do j = i+1,noa
                  do a = 1,nua
                    do b = a+1,nua
                      denom = fA_oo(i,i) + fA_oo(j,j) - fA_vv(a,a) - fA_vv(b,b)
                      t2a(b,a,j,i) = t2a(b,a,j,i) + X2A(b,a,j,i)/(denom-shift)
                      t2a(a,b,j,i) = -t2a(b,a,j,i)
                      t2a(b,a,i,j) = -t2a(b,a,j,i)
                      t2a(a,b,i,j) = t2a(b,a,j,i)
                    end do
                  end do
                end do
              end do

      end subroutine update_t2a

      subroutine update_t2b(t2b,X2B,fA_oo,fA_vv,fB_oo,fB_vv,shift,noa,nua,nob,nub)

              implicit none

              integer, intent(in) :: noa, nua, nob, nub
              real(8), intent(in) :: fA_oo(1:noa,1:noa), fA_vv(1:nua,1:nua), &
                                  fB_oo(1:nob,1:nob), fB_vv(1:nub,1:nub), &
                                  X2B(1:nua,1:nub,1:noa,1:nob), shift               
              real(8), intent(inout) :: t2b(1:nua,1:nub,1:noa,1:nob)
              !f2py intent(in,out) :: t2b(0:nua-1,0:nub-1,0:noa-1,0:nob-1)
              integer :: i, j, a, b
              real(8) :: denom

              do j = 1,nob
                do i = 1,noa
                  do b = 1,nub
                    do a = 1,nua
                      denom = fA_oo(i,i) + fB_oo(j,j) - fA_vv(a,a) - fB_vv(b,b)
                      t2b(a,b,i,j) = t2b(a,b,i,j) + X2B(a,b,i,j)/(denom-shift)
                    end do
                  end do
                end do
              end do

      end subroutine update_t2b

      subroutine update_t2c(t2c,X2C,fB_oo,fB_vv,shift,nob,nub)

              implicit none

              integer, intent(in) :: nob, nub
              real(8), intent(in) :: fB_oo(1:nob,1:nob), fB_vv(1:nub,1:nub), &
                                  X2C(1:nub,1:nub,1:nob,1:nob), shift               
              real(8), intent(inout) :: t2c(1:nub,1:nub,1:nob,1:nob)
              !f2py intent(in,out) :: t2c(0:nub-1,0:nub-1,0:nob-1,0:nob-1)
              integer :: i, j, a, b
              real(8) :: denom

              do i = 1,nob
                do j = i+1,nob
                  do a = 1,nub
                    do b = a+1,nub
                      denom = fB_oo(i,i) + fB_oo(j,j) - fB_vv(a,a) - fB_vv(b,b)
                      t2c(b,a,j,i) = t2c(b,a,j,i) + X2C(b,a,j,i)/(denom-shift)
                      t2c(a,b,j,i) = -t2c(b,a,j,i)
                      t2c(b,a,i,j) = -t2c(b,a,j,i)
                      t2c(a,b,i,j) = t2c(b,a,j,i)
                    end do
                  end do
                end do
              end do

      end subroutine update_t2c

      subroutine update_t3a(t3a,X3A,fA_oo,fA_vv,shift,noa,nua)

              implicit none

              integer, intent(in) :: noa, nua
              real(8), intent(in) :: fA_oo(1:noa,1:noa), fA_vv(1:nua,1:nua), &
                                  X3A(1:nua,1:nua,1:nua,1:noa,1:noa,1:noa), shift               
              real(8), intent(inout) :: t3a(1:nua,1:nua,1:nua,1:noa,1:noa,1:noa)
              !f2py intent(in,out) :: t3a(0:nua-1,0:nua-1,0:nua-1,0:noa-1,0:noa-1,0:noa-1)
              integer :: i, j, k, a, b, c, ii, jj, kk, aa, bb, cc
              real(8) :: denom

              do ii = 1,noa
                  do jj = ii+1,noa
                      do kk = jj+1,noa
                          do aa = 1,nua
                              do bb = aa+1,nua
                                  do cc = bb+1,nua

                                      A = cc; B = bb; C = aa;
                                      I = kk; J = jj; K = ii;
                                      
                                      denom = fA_oo(I,I)+fA_oo(J,J)+fA_oo(K,K)-fA_vv(A,A)-fA_vv(B,B)-fA_vv(C,C)

                                      t3a(A,B,C,I,J,K) = t3a(A,B,C,I,J,K) + X3A(A,B,C,I,J,K)/(denom-shift)                            
                                      t3a(A,B,C,K,I,J) = t3a(A,B,C,I,J,K)
                                      t3a(A,B,C,J,K,I) = t3a(A,B,C,I,J,K)
                                      t3a(A,B,C,I,K,J) = -t3a(A,B,C,I,J,K)
                                      t3a(A,B,C,J,I,K) = -t3a(A,B,C,I,J,K)
                                      t3a(A,B,C,K,J,I) = -t3a(A,B,C,I,J,K)
                                      
                                      t3a(B,A,C,I,J,K) = -t3a(A,B,C,I,J,K)
                                      t3a(B,A,C,K,I,J) = -t3a(A,B,C,I,J,K)
                                      t3a(B,A,C,J,K,I) = -t3a(A,B,C,I,J,K)
                                      t3a(B,A,C,I,K,J) = t3a(A,B,C,I,J,K)
                                      t3a(B,A,C,J,I,K) = t3a(A,B,C,I,J,K)
                                      t3a(B,A,C,K,J,I) = t3a(A,B,C,I,J,K)
                                      
                                      t3a(A,C,B,I,J,K) = -t3a(A,B,C,I,J,K)
                                      t3a(A,C,B,K,I,J) = -t3a(A,B,C,I,J,K)
                                      t3a(A,C,B,J,K,I) = -t3a(A,B,C,I,J,K)
                                      t3a(A,C,B,I,K,J) = t3a(A,B,C,I,J,K)
                                      t3a(A,C,B,J,I,K) = t3a(A,B,C,I,J,K)
                                      t3a(A,C,B,K,J,I) = t3a(A,B,C,I,J,K)
                                      
                                      t3a(C,B,A,I,J,K) = -t3a(A,B,C,I,J,K)
                                      t3a(C,B,A,K,I,J) = -t3a(A,B,C,I,J,K)
                                      t3a(C,B,A,J,K,I) = -t3a(A,B,C,I,J,K)
                                      t3a(C,B,A,I,K,J) = t3a(A,B,C,I,J,K)
                                      t3a(C,B,A,J,I,K) = t3a(A,B,C,I,J,K)
                                      t3a(C,B,A,K,J,I) = t3a(A,B,C,I,J,K)
                                      
                                      t3a(B,C,A,I,J,K) = t3a(A,B,C,I,J,K)
                                      t3a(B,C,A,K,I,J) = t3a(A,B,C,I,J,K)
                                      t3a(B,C,A,J,K,I) = t3a(A,B,C,I,J,K)
                                      t3a(B,C,A,I,K,J) = -t3a(A,B,C,I,J,K)
                                      t3a(B,C,A,J,I,K) = -t3a(A,B,C,I,J,K)
                                      t3a(B,C,A,K,J,I) = -t3a(A,B,C,I,J,K)
                                      
                                      t3a(C,A,B,I,J,K) = t3a(A,B,C,I,J,K)
                                      t3a(C,A,B,K,I,J) = t3a(A,B,C,I,J,K)
                                      t3a(C,A,B,J,K,I) = t3a(A,B,C,I,J,K)
                                      t3a(C,A,B,I,K,J) = -t3a(A,B,C,I,J,K)
                                      t3a(C,A,B,J,I,K) = -t3a(A,B,C,I,J,K)
                                      t3a(C,A,B,K,J,I) = -t3a(A,B,C,I,J,K)
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do

      end subroutine update_t3a

      subroutine update_t3b(t3b,X3B,fA_oo,fA_vv,fB_oo,fB_vv,shift,noa,nua,nob,nub)

              implicit none

              integer, intent(in) :: noa, nua, nob, nub
              real(8), intent(in) :: fA_oo(1:noa,1:noa), fA_vv(1:nua,1:nua), &
                                  fB_oo(1:nob,1:nob), fB_vv(1:nub,1:nub), &
                                  X3B(1:nua,1:nua,1:nub,1:noa,1:noa,1:nob), shift               
              real(8), intent(inout) :: t3b(1:nua,1:nua,1:nub,1:noa,1:noa,1:nob)
              !f2py intent(in,out) :: t3b(0:nua-1,0:nua-1,0:nub-1,0:noa-1,0:noa-1,0:nob-1)
              integer :: i, j, k, a, b, c, ii, jj, kk, aa, bb, cc
              real(8) :: denom

              do ii = 1,noa
                  do jj = ii+1,noa
                      do kk = 1,nob
                          do aa = 1,nua
                              do bb = aa+1,nua
                                  do cc = 1,nub
                  
                                      a = bb; b = aa; c = cc;
                                      i = jj; j = ii; k = kk;

                                      denom = fA_oo(i,i)+fA_oo(j,j)+fB_oo(k,k)-fA_vv(a,a)-fA_vv(b,b)-fB_vv(c,c)
                                      t3b(a,b,c,i,j,k) = t3b(a,b,c,i,j,k) + X3B(a,b,c,i,j,k)/(denom-shift)
                                      t3b(b,a,c,i,j,k) = -t3b(a,b,c,i,j,k)
                                      t3b(a,b,c,j,i,k) = -t3b(a,b,c,i,j,k)
                                      t3b(b,a,c,j,i,k) = t3b(a,b,c,i,j,k)

                                  end do
                              end do
                          end do
                      end do
                  end do
              end do

      end subroutine update_t3b

      subroutine update_t3c(t3c,X3C,fA_oo,fA_vv,fB_oo,fB_vv,shift,noa,nua,nob,nub)

              implicit none

              integer, intent(in) :: noa, nua, nob, nub
              real(8), intent(in) :: fA_oo(1:noa,1:noa), fA_vv(1:nua,1:nua), &
                                  fB_oo(1:nob,1:nob), fB_vv(1:nub,1:nub), &
                                  X3C(1:nua,1:nub,1:nub,1:noa,1:nob,1:nob), shift               
              real(8), intent(inout) :: t3c(1:nua,1:nub,1:nub,1:noa,1:nob,1:nob)
              !f2py intent(in,out) :: t3c(0:nua-1,0:nub-1,0:nub-1,0:noa-1,0:nob-1,0:nob-1)
              integer :: i, j, k, a, b, c, ii, jj, kk, aa, bb, cc
              real(8) :: denom

              do ii = 1,noa
                  do jj = 1,nob
                      do kk = jj+1,nob
                          do aa = 1,nua
                              do bb = 1,nub
                                  do cc = bb+1,nub
                  
                                      a = aa; b = cc; c = bb;
                                      i = ii; j = kk; k = jj;

                                      denom = fA_oo(i,i)+fB_oo(j,j)+fB_oo(k,k)-fA_vv(a,a)-fB_vv(b,b)-fB_vv(c,c)
                                      t3c(a,b,c,i,j,k) = t3c(a,b,c,i,j,k) + X3C(a,b,c,i,j,k)/(denom-shift)
                                      t3c(a,c,b,i,j,k) = -t3c(a,b,c,i,j,k)
                                      t3c(a,b,c,i,k,j) = -t3c(a,b,c,i,j,k)
                                      t3c(a,c,b,i,k,j) = t3c(a,b,c,i,j,k)

                                  end do
                              end do
                          end do
                      end do
                  end do
              end do

      end subroutine update_t3c

      subroutine update_t3d(t3d,X3D,fB_oo,fB_vv,shift,nob,nub)

              implicit none

              integer, intent(in) :: nob, nub
              real(8), intent(in) :: fB_oo(1:nob,1:nob), fB_vv(1:nub,1:nub), &
                                  X3D(1:nub,1:nub,1:nub,1:nob,1:nob,1:nob), shift               
              real(8), intent(inout) :: t3d(1:nub,1:nub,1:nub,1:nob,1:nob,1:nob)
              !f2py intent(in,out) :: t3d(0:nub-1,0:nub-1,0:nub-1,0:nob-1,0:nob-1,0:nob-1)
              integer :: i, j, k, a, b, c, ii, jj, kk, aa, bb, cc
              real(8) :: denom

              do ii = 1,nob
                  do jj = ii+1,nob
                      do kk = jj+1,nob
                          do aa = 1,nub
                              do bb = aa+1,nub
                                  do cc = bb+1,nub

                                      A = cc; B = bb; C = aa;
                                      I = kk; J = jj; K = ii;
                                      
                                      denom = fB_oo(I,I)+fB_oo(J,J)+fB_oo(K,K)-fB_vv(A,A)-fB_vv(B,B)-fB_vv(C,C)

                                      t3d(A,B,C,I,J,K) = t3d(A,B,C,I,J,K) + X3D(A,B,C,I,J,K)/(denom-shift)                            
                                      t3d(A,B,C,K,I,J) = t3d(A,B,C,I,J,K)
                                      t3d(A,B,C,J,K,I) = t3d(A,B,C,I,J,K)
                                      t3d(A,B,C,I,K,J) = -t3d(A,B,C,I,J,K)
                                      t3d(A,B,C,J,I,K) = -t3d(A,B,C,I,J,K)
                                      t3d(A,B,C,K,J,I) = -t3d(A,B,C,I,J,K)
                                      
                                      t3d(B,A,C,I,J,K) = -t3d(A,B,C,I,J,K)
                                      t3d(B,A,C,K,I,J) = -t3d(A,B,C,I,J,K)
                                      t3d(B,A,C,J,K,I) = -t3d(A,B,C,I,J,K)
                                      t3d(B,A,C,I,K,J) = t3d(A,B,C,I,J,K)
                                      t3d(B,A,C,J,I,K) = t3d(A,B,C,I,J,K)
                                      t3d(B,A,C,K,J,I) = t3d(A,B,C,I,J,K)
                                      
                                      t3d(A,C,B,I,J,K) = -t3d(A,B,C,I,J,K)
                                      t3d(A,C,B,K,I,J) = -t3d(A,B,C,I,J,K)
                                      t3d(A,C,B,J,K,I) = -t3d(A,B,C,I,J,K)
                                      t3d(A,C,B,I,K,J) = t3d(A,B,C,I,J,K)
                                      t3d(A,C,B,J,I,K) = t3d(A,B,C,I,J,K)
                                      t3d(A,C,B,K,J,I) = t3d(A,B,C,I,J,K)
                                      
                                      t3d(C,B,A,I,J,K) = -t3d(A,B,C,I,J,K)
                                      t3d(C,B,A,K,I,J) = -t3d(A,B,C,I,J,K)
                                      t3d(C,B,A,J,K,I) = -t3d(A,B,C,I,J,K)
                                      t3d(C,B,A,I,K,J) = t3d(A,B,C,I,J,K)
                                      t3d(C,B,A,J,I,K) = t3d(A,B,C,I,J,K)
                                      t3d(C,B,A,K,J,I) = t3d(A,B,C,I,J,K)
                                      
                                      t3d(B,C,A,I,J,K) = t3d(A,B,C,I,J,K)
                                      t3d(B,C,A,K,I,J) = t3d(A,B,C,I,J,K)
                                      t3d(B,C,A,J,K,I) = t3d(A,B,C,I,J,K)
                                      t3d(B,C,A,I,K,J) = -t3d(A,B,C,I,J,K)
                                      t3d(B,C,A,J,I,K) = -t3d(A,B,C,I,J,K)
                                      t3d(B,C,A,K,J,I) = -t3d(A,B,C,I,J,K)
                                      
                                      t3d(C,A,B,I,J,K) = t3d(A,B,C,I,J,K)
                                      t3d(C,A,B,K,I,J) = t3d(A,B,C,I,J,K)
                                      t3d(C,A,B,J,K,I) = t3d(A,B,C,I,J,K)
                                      t3d(C,A,B,I,K,J) = -t3d(A,B,C,I,J,K)
                                      t3d(C,A,B,J,I,K) = -t3d(A,B,C,I,J,K)
                                      t3d(C,A,B,K,J,I) = -t3d(A,B,C,I,J,K)
                                  end do
                              end do
                          end do
                      end do
                  end do
              end do

      end subroutine update_t3d

      subroutine update_R(r1a,r1b,r2a,r2b,r2c,omega,H1A_oo,H1A_vv,H1B_oo,H1B_vv,shift,noa,nua,nob,nub)

              implicit none

              integer, intent(in) :: noa, nua, nob, nub
              real(8), intent(in) :: H1A_oo(1:noa,1:noa), H1A_vv(1:nua,1:nua), &
                                  H1B_oo(1:nob,1:nob), H1B_vv(1:nub,1:nub), shift, &
                                  omega  
              real(8), intent(inout) :: r1a(1:nua,1:noa)
              !f2py intent(in,out) :: r1a(0:nua-1,0:noa-1) 
              real(8), intent(inout) :: r1b(1:nub,1:nob)
              !f2py intent(in,out) :: r1b(0:nub-1,0:nob-1)   
              real(8), intent(inout) :: r2a(1:nua,1:nua,1:noa,1:noa)
              !f2py intent(in,out) :: r2a(0:nua-1,0:nua-1,0:noa-1,0:noa-1)        
              real(8), intent(inout) :: r2b(1:nua,1:nub,1:noa,1:nob)
              !f2py intent(in,out) :: r2b(0:nua-1,0:nub-1,0:noa-1,0:nob-1)
              real(8), intent(inout) :: r2c(1:nub,1:nub,1:nob,1:nob)
              !f2py intent(in,out) :: r2c(0:nub-1,0:nub-1,0:nob-1,0:nob-1)
              integer :: i, j, a, b
              real(8) :: denom

              do i = 1,noa
                do a = 1,nua
                  denom = H1A_vv(a,a) - H1A_oo(i,i)
                  r1a(a,i) = r1a(a,i)/(omega-denom+shift)
                end do
              end do

              do i = 1,nob
                do a = 1,nub
                  denom = H1B_vv(a,a) - H1B_oo(i,i)
                  r1b(a,i) = r1b(a,i)/(omega-denom+shift)
                end do
              end do

              do i = 1,noa
                do j = i+1,noa
                  do a = 1,nua
                    do b = a+1,nua
                      denom = H1A_vv(a,a) + H1A_vv(b,b) - H1A_oo(i,i) - H1A_oo(j,j)
                      r2a(b,a,j,i) = r2a(b,a,j,i)/(omega-denom+shift)
                      r2a(a,b,j,i) = -r2a(b,a,j,i)
                      r2a(b,a,i,j) = -r2a(b,a,j,i)
                      r2a(a,b,i,j) = r2a(b,a,j,i)
                    end do
                  end do
                end do
              end do

              do j = 1,nob
                do i = 1,noa
                  do b = 1,nub
                    do a = 1,nua
                      denom = H1A_vv(a,a) + H1B_vv(b,b) - H1A_oo(i,i) - H1B_oo(j,j)
                      r2b(a,b,i,j) = r2b(a,b,i,j)/(omega-denom+shift)
                    end do
                  end do
                end do
              end do

              do i = 1,nob
                do j = i+1,nob
                  do a = 1,nub
                    do b = a+1,nub
                      denom = H1B_vv(a,a) + H1B_vv(b,b) - H1B_oo(i,i) - H1B_oo(j,j)
                      r2c(b,a,j,i) = r2c(b,a,j,i)/(omega-denom+shift)
                      r2c(a,b,j,i) = -r2c(b,a,j,i)
                      r2c(b,a,i,j) = -r2c(b,a,j,i)
                      r2c(a,b,i,j) = r2c(b,a,j,i)
                    end do
                  end do
                end do
              end do

      end subroutine update_R

      subroutine update_L(l1a,l1b,l2a,l2b,l2c,X1A,X1B,X2A,X2B,X2C,H1A_oo,H1A_vv,H1B_oo,H1B_vv,shift,noa,nua,nob,nub)

              implicit none

              integer, intent(in) :: noa, nua, nob, nub
              real(8), intent(in) :: H1A_oo(1:noa,1:noa), H1A_vv(1:nua,1:nua), &
                                  H1B_oo(1:nob,1:nob), H1B_vv(1:nub,1:nub), shift  
              real(8), intent(inout) :: l1a(1:nua,1:noa)
              !f2py intent(in,out) :: l1a(0:nua-1,0:noa-1) 
              real(8), intent(inout) :: l1b(1:nub,1:nob)
              !f2py intent(in,out) :: l1b(0:nub-1,0:nob-1)   
              real(8), intent(inout) :: l2a(1:nua,1:nua,1:noa,1:noa)
              !f2py intent(in,out) :: l2a(0:nua-1,0:nua-1,0:noa-1,0:noa-1)        
              real(8), intent(inout) :: l2b(1:nua,1:nub,1:noa,1:nob)
              !f2py intent(in,out) :: l2b(0:nua-1,0:nub-1,0:noa-1,0:nob-1)
              real(8), intent(inout) :: l2c(1:nub,1:nub,1:nob,1:nob)
              !f2py intent(in,out) :: l2c(0:nub-1,0:nub-1,0:nob-1,0:nob-1)
              real(8), intent(in) :: X1A(1:nua,1:noa)
              !f2py intent(in) :: X1A(0:nua-1,0:noa-1) 
              real(8), intent(in) :: X1B(1:nub,1:nob)
              !f2py intent(in) :: X1B(0:nub-1,0:nob-1)   
              real(8), intent(in) :: X2A(1:nua,1:nua,1:noa,1:noa)
              !f2py intent(in) :: X2A(0:nua-1,0:nua-1,0:noa-1,0:noa-1)        
              real(8), intent(in) :: X2B(1:nua,1:nub,1:noa,1:nob)
              !f2py intent(in) :: X2B(0:nua-1,0:nub-1,0:noa-1,0:nob-1)
              real(8), intent(in) :: X2C(1:nub,1:nub,1:nob,1:nob)
              !f2py intent(in) :: X2C(0:nub-1,0:nub-1,0:nob-1,0:nob-1)
              integer :: i, j, a, b
              real(8) :: denom

              do i = 1,noa
                do a = 1,nua
                  denom = H1A_vv(a,a) - H1A_oo(i,i)
                  l1a(a,i) = l1a(a,i) - X1A(a,i)/(denom+shift)
                end do
              end do

              do i = 1,nob
                do a = 1,nub
                  denom = H1B_vv(a,a) - H1B_oo(i,i)
                  l1b(a,i) = l1b(a,i) - X1B(a,i)/(denom+shift)
                end do
              end do

              do i = 1,noa
                do j = i+1,noa
                  do a = 1,nua
                    do b = a+1,nua
                      denom = H1A_vv(a,a) + H1A_vv(b,b) - H1A_oo(i,i) - H1A_oo(j,j)
                      l2a(b,a,j,i) = l2a(b,a,j,i) - X2A(b,a,j,i)/(denom+shift)
                      l2a(a,b,j,i) = -l2a(b,a,j,i)
                      l2a(b,a,i,j) = -l2a(b,a,j,i)
                      l2a(a,b,i,j) = l2a(b,a,j,i)
                    end do
                  end do
                end do
              end do

              do j = 1,nob
                do i = 1,noa
                  do b = 1,nub
                    do a = 1,nua
                      denom = H1A_vv(a,a) + H1B_vv(b,b) - H1A_oo(i,i) - H1B_oo(j,j)
                      l2b(a,b,i,j) = l2b(a,b,i,j) - X2B(a,b,i,j)/(denom+shift)
                    end do
                  end do
                end do
              end do

              do i = 1,nob
                do j = i+1,nob
                  do a = 1,nub
                    do b = a+1,nub
                      denom = H1B_vv(a,a) + H1B_vv(b,b) - H1B_oo(i,i) - H1B_oo(j,j)
                      l2c(b,a,j,i) = l2c(b,a,j,i) - X2C(b,a,j,i)/(denom+shift)
                      l2c(a,b,j,i) = -l2c(b,a,j,i)
                      l2c(b,a,i,j) = -l2c(b,a,j,i)
                      l2c(a,b,i,j) = l2c(b,a,j,i)
                    end do
                  end do
                end do
              end do

      end subroutine update_L


end module cc_loops
