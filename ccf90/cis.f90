module cis

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use dgeev_module, only: eigh
        use sort_module, only: argsort

        implicit none

        contains

              subroutine cis_mat(sys,fA,fB,vA,vB,vC,C,omega)
                
                      type(sys_t), intent(in) :: sys
                      type(e1int_t), intent(in) :: fA, fB
                      type(e2int_t), intent(in) :: vA, vB, vC
                      real, allocatable, intent(out) :: C(:,:), omega(:)
                      integer :: n1a, n1b, cnt1, cnt2, a, b, i, j
                      real, allocatable :: Hcis(:,:)
                      real :: val

                      n1a = sys%Nunocc_a*sys%Nocc_a
                      n1b = sys%Nunocc_b*sys%Nocc_b

                      allocate(Hcis(n1a+n1b,n1a+n1b),C(n1a+n1b,n1a+n1b),omega(n1a+n1b))

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
                                  Hcis(cnt1,cnt2) = val
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
                                  Hcis(cnt1,n1a+cnt2) = val
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
                                  Hcis(n1a+cnt1,cnt2) = val
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
                                  Hcis(n1a+cnt1,n1a+cnt2) = val
                               end do
                            end do
                         end do
                      end do


                      call eigh(Hcis,C,omega)
                      !Cvec1A = C(1:n1a,:)
                      !Cvec1B = C(n1a+1:n1a+n1b,:)
                      deallocate(Hcis)

              end subroutine cis_mat


end module cis
