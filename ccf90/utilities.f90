module utilities


        implicit none

        contains

                function spatial_idx(p) result(pout)

                        integer, intent(in) :: p
                        integer :: pout

                        if (mod(p,2) == 1) then ! alpha
                                pout = (p + 1)/2
                        else
                                pout = p/2
                        end if

                end function spatial_idx

                subroutine convert_to_spinorb(t1a,t1b,t2a,t2b,t2c,t1,t2)

                        real, intent(in) :: t1a(:,:), t1b(:,:), &
                                            t2a(:,:,:,:), t2b(:,:,:,:), t2c(:,:,:,:)
                        real, allocatable, intent(out) :: t1(:,:), t2(:,:,:,:)
                        integer :: nua, noa, nub, nob, nu, no, i, j, a, b
                        
                        nua = size(t1a,1)
                        noa = size(t1a,2)
                        nub = size(t1b,1)
                        nob = size(t1b,2)

                        no = noa + nob
                        nu = nua + nub

                        allocate(t1(nu,no),t2(nu,nu,no,no))

                        do i = 1,no
                           do a = 1,nu
                              if ( (mod(i,2) == 1) .AND. (mod(a,2) == 1) ) then
                                      t1(a,i) = t1a(spatial_idx(a),spatial_idx(i))
                              else
                                      t1(a,i) = t1b(spatial_idx(a),spatial_idx(i))
                              end if
                           end do
                        end do

                        do i = 1,no
                           do j = i+1,no
                              do a = 1,nu
                                 do b = a+1,nu

                                    if ( (mod(i,2) == 1) .AND. (mod(a,2) == 1) .AND.&
                                            (mod(j,2) == 1) .AND. (mod(b,2) == 1) ) then ! aa
                                            t2(a,b,i,j) = t2a(spatial_idx(a),spatial_idx(b),spatial_idx(i),spatial_idx(j))
                                    else if ( (mod(i,2) == mod(a,2)) .AND. (mod(j,2) == mod(b,2)) ) then ! ab
                                            t2(a,b,i,j) = t2b(spatial_idx(a),spatial_idx(b),spatial_idx(i),spatial_idx(j))
                                    else if ( (mod(i,2) == 0) .AND. (mod(a,2) == 0) .AND.&
                                            (mod(j,2) == 0) .AND. (mod(b,2) == 0) ) then ! bb
                                            t2(a,b,i,j) = t2c(spatial_idx(a),spatial_idx(b),spatial_idx(i),spatial_idx(j))
                                    end if

                                    t2(b,a,i,j) = -t2(a,b,i,j)
                                    t2(a,b,j,i) = -t2(a,b,i,j)
                                    t2(b,a,j,i) = t2(a,b,i,j)
                                 end do
                              end do
                           end do
                        end do

                end subroutine convert_to_spinorb
                                   

end module utilities 


