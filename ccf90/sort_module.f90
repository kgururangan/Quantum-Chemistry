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


       subroutine get_intersection_sorted(arr1,arr2,inter,idx1,idx2)

        ! Finds the intersection between two SORTED integer arrays
        ! Returns the intsersection as a sorted array as well as the
        ! indices of intersection in each sorted array

                integer, intent(in) :: arr1(:), arr2(:)
                integer, allocatable, intent(out) :: inter(:), idx1(:), idx2(:)
                integer :: m, n, k, i, j, ct
                integer, allocatable :: temp(:), tmp1(:), tmp2(:)

                m = size(arr1)
                n = size(arr2)
                k = min(m,n)

                allocate(temp(k),tmp1(k),tmp2(k))
        
                i = 1
                j = 1
                ct = 0
                do while ( (i <= m) .AND. (j <=n) )
                        if (arr1(i) < arr2(j)) then
                                i = i + 1
                        else if (arr2(j) < arr1(i)) then
                                j = j + 1
                        else
                                ct = ct + 1
                                tmp1(ct) = i
                                tmp2(ct) = j
                                temp(ct) = arr2(j)
                                j = j + 1
                                i = i + 1
                        end if
                 end do

                 allocate(inter(ct),idx1(ct),idx2(ct))
                 inter = temp(1:ct)
                 idx1 = tmp1(1:ct)
                 idx2 = tmp2(1:ct)
                 deallocate(temp,tmp1,tmp2)

        end subroutine get_intersection_sorted

        subroutine get_intersection(arr1,arr2,inter,idx1,idx2)

        ! Finds the intersection between two unsorted integer arrays
        ! Returns the intsersection as a sorted array as well as the
        ! indices of intersection in the original array

                integer, intent(in) :: arr1(:), arr2(:)
                integer, allocatable, intent(out) :: inter(:), idx1(:), idx2(:)
                integer :: m, n, k, i, j, ct
                integer, allocatable :: temp(:), tmp1(:), tmp2(:), arr1sort(:), arr2sort(:), &
                                        id1(:), id2(:)

                m = size(arr1)
                n = size(arr2)
                k = min(m,n)

                allocate(temp(k),tmp1(k),tmp2(k),arr1sort(m),arr2sort(n),id1(m),id2(n))

                id1 = argsort_int(arr1)
                id2 = argsort_int(arr2)
                arr1sort = arr1(id1)
                arr2sort = arr2(id2)
        
                i = 1
                j = 1
                ct = 0
                do while ( (i <= m) .AND. (j <=n) )
                        if (arr1sort(i) < arr2sort(j)) then
                                i = i + 1
                        else if (arr2sort(j) < arr1sort(i)) then
                                j = j + 1
                        else
                                ct = ct + 1
                                tmp1(ct) = id1(i)
                                tmp2(ct) = id2(j)
                                temp(ct) = arr2sort(j)
                                j = j + 1
                                i = i + 1
                        end if
                 end do

                 allocate(inter(ct),idx1(ct),idx2(ct))
                 inter = temp(1:ct)
                 idx1 = tmp1(1:ct)
                 idx2 = tmp2(1:ct)
                 deallocate(temp,tmp1,tmp2,id1,id2,arr1sort,arr2sort)

        end subroutine get_intersection


end module sort_module
