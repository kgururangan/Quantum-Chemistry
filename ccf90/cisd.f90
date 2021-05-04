module cisd

        use system_types, only: sys_t
        use integral_types, only: e1int_t, e2int_t
        use dgeev_module, only: eigh
        use sort_module, only: argsort

        implicit none

        contains

              subroutine cisd_mat(sys,fA,fB,vA,vB,vC,nact_h,nact_p,C,omega)
                
                      type(sys_t), intent(in) :: sys
                      type(e1int_t), intent(in) :: fA, fB
                      type(e2int_t), intent(in) :: vA, vB, vC
                      real, allocatable, intent(out) :: C(:,:), omega(:)
                      integer, intent(in) :: nact_h, nact_p
                      integer :: n1a, n1b, n2a, n2b, n2c, cnt1, cnt2, a, b, i, j, k, l, c, d
                      real, allocatable :: HDD(:,:), HDS(:,:), HSD(:,:), HSS(:,:)
                      integer, allocatable :: act_h_alpha_rng(nact_h), act_h_beta_rng(nact_h), &
                                              act_p_alpha_rng(nact_p), act_p_beta_rng(nact_p)
                      real :: val

                      n1a = sys%Nunocc_a*sys%Nocc_a
                      n1b = sys%Nunocc_b*sys%Nocc_b
                      n2a = sys%Nunocc_a**2 * sys%Nocc_a**2
                      n2b = sys%Nunocc_a*sys%Nunocc_b*sys%Nocc_a*sys%Nocc_b
                      n2c = sys%Nunocc_b**2 * sys%Nocc_b**2
                      ndim = n1a + n1b + n2a + n2b + n2c

                      nact_h_alpha_rng = (/do i = sys%Nocc_a-nact_h+1,sys%Nocc_a/)
                      nact_h_beta_rng = (/do i = sys%Nocc_b-nact_h+1,sys%Nocc_b/)
                      nact_p_alpha_rng = (/do i = 1,nact_p/)
                      nact_p_beta_rng = (/do i = 1,nact_p/)


                      allocate(HSS(n1a+n1b,n1a+n1b),HSD(n1a+n1b,n2a+n2b+n2c),&
                               HDS(n2a+n2b+n2c,n1a+n1b),HDD(n2a+n2b+n2c,n2a+n2b+n2c),&
                               C(ndim,ndim),omega(ndim))

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
                                  HSS(cnt1,cnt2) = val
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
                                  HSS(cnt1,n1a+cnt2) = val
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
                                  HSS(n1a+cnt1,cnt2) = val
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
                                  HSS(n1a+cnt1,n1a+cnt2) = val
                               end do
                            end do
                         end do
                      end do

                      cnt1 = 0
                      do i = 1,sys%Nocc_a
                         do j = i+1,sys%Nocc_a
                            do a = 1,sys%Nunocc_a
                               do b = a+1,sys%Nunocc_a

                                  if ((isactive((/i,j/),act_h_alpha_rng)) .AND. &
                                          (isactive((/a,b/),act_p_alpha_rng))) then


                                        cnt1 = cnt1 + 1
                                        cnt2 = 0

                                        do k = 1,sys%Nocc_a
                                           do l = k+1,sys%Nocc_a
                                              do c = 1,sys%Nunocc_a
                                                 do d = c+1,sys%Nunocc_a

                
                                                    if ((isactive((/i,j/),act_h_alpha_rng)) .AND. &
                                                            (isactive((/a,b/),act_p_alpha_rng))) then
                                                         
                                                        cnt2 = cnt + 1
                                                        val = 0.0

                                                    else
                                                        cycle
                                                    end if

                                                 end do
                                              end do
                                           end do
                                        end do


                                     else
                                        cycle
                                     end if


                                end do
                             end do
                          end do
                       end do

                                              
                                              



                      !call eigh(Hcis,C,omega)
                      !deallocate(Hcis)

              end subroutine cisd_mat

              function isactive(arr,act_range) result(bool)

                      integer, intent(in) :: arr(:)
                      logical :: bool
                      integer :: ct, i

                      ct = 0
                      do i = 1,size(arr)
                         if ( (arr(i) >= act_range(1)) .AND. (arr(i) <= act_range(2)) ) then
                                 ct = ct + 1
                         end if
                      end do

                      if (ct >= 1)
                         bool = .TRUE.
                      else
                         bool = .FALSE.
                      end if

              end function isactive





              subroutine get_orb_diff(D1,D2,spinA,spinB,diff_oa,diff_ua,diff_ob,diff_ub,sgn)

                      use sort_module, only: argsort_int, get_intersection_sorted
                      use permutils, only: permsign, reorder_max_coincidence

                      integer, intent(in) :: D1(:), D2(:)
                      integer, allocatable, intent(out) :: diff_oa(:), diff_ua(:), diff_ob(:), diff_ub(:)
                      real, intent(out) :: sgn
                      character, intent(in) :: spinA(1), spinB(1)
                      integer, allocatable :: i1a(:), i1b(:), i2a(:), i2b(:), a1a(:), a1b(:), a2a(:), a2b(:)
                      integer :: exc1, exc2, p, q
                      real :: sgn_o, sgn_u

                      exc1 = size(D1)/2
                      exc2 = size(D2)/2

                      if (exc1 == 1) then
                              if (spinA == 'A') then
                                      allocate(i1a(1),a1a(1))
                                      i1a = D1(1)
                                      a1a = D1(2)

                              else
                                      allocate(i1b(1),a1b(1))
                                      i1b = D1(1)
                                      a1b = D1(2)
                              end if

                      else if (exc1 == 2) then
                              if (spinA == 'A') then
                                      allocate(i1a(2),a1a(2))
                                      i1a = D1(1:2)
                                      a1a = D1(3:4)
                              else if (spinA == 'B') then
                                      allocate(i1a(1),i1b(1),a1a(1),a1b(1))
                                      i1a = D1(1)
                                      i1b = D1(2)
                                      a1a = D1(3)
                                      a1b = D1(4)
                              else
                                      allocate(i1b(2),a1b(2))
                                      i1b = D1(1:2)
                                      a1b = D1(3:4)
                              end if
                       else
                               print*,'EXCITATION RANK NOT SUPPORTED'
                       end if

                      if (exc2 == 1) then
                              if (spinB == 'A') then
                                      allocate(i2a(1),a2a(1))
                                      i2a = D2(1)
                                      a2a = D2(2)

                              else
                                      allocate(i2b(1),a2b(1))
                                      i2b = D2(1)
                                      a2b = D2(2)
                              end if

                      else if (exc2 == 2) then
                              if (spinB == 'A') then
                                      allocate(i2a(2),a2a(2))
                                      i2a = D2(1:2)
                                      a2a = D2(3:4)
                              else if (spinB == 'B') then
                                      allocate(i2a(1),i2b(1),a2a(1),a2b(1))
                                      i2a = D2(1)
                                      i2b = D2(2)
                                      a2a = D2(3)
                                      a2b = D2(4)
                              else
                                      allocate(i2b(2),a2b(2))
                                      i2b = D2(1:2)
                                      a2b = D2(3:4)
                              end if
                       else
                               print*,'EXCITATION RANK NOT SUPPORTED'
                       end if


              end subroutine get_orb_diff






end module cisd
