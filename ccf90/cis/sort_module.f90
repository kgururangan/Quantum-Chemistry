module sort_module
 
implicit none
 
contains
 
        function argsort(a) result(idx)

                real, intent(in) :: a(:)
                integer :: i, j, idx(size(a)), temp

                do i = 1,size(a)
                    idx(i) = i
                end do

                do j = 0,size(a)-1
                    do i = 1,size(a)-j-1 
                            if (a(idx(i)) > a(idx(i+1))) then
                                    temp   = idx(i)
                                    idx(i) = idx(i+1)
                                    idx(i+1) = temp
                            end if
                    end do
                end do

        end function argsort

        function argsort_int(a) result(idx)

                integer, intent(in) :: a(:)
                integer :: i, j, idx(size(a)), temp

                do i = 1,size(a)
                    idx(i) = i
                end do

                do j = 0,size(a)-1
                    do i = 1,size(a)-j-1 
                            if (a(idx(i)) > a(idx(i+1))) then
                                    temp   = idx(i)
                                    idx(i) = idx(i+1)
                                    idx(i+1) = temp
                            end if
                    end do
                end do

        end function argsort_int

end module sort_module
