module permutils

    implicit none

    contains


    subroutine reorder_stripe(rank, mat_shape, mat_size, perm_str, A, B)

        integer, intent(in) :: rank
        integer, intent(in) :: mat_shape(:)
        integer, intent(in) :: mat_size
        character(len=*), intent(in) :: perm_str

        real, intent(in) :: A(mat_size)
        real, intent(inout) :: B(mat_size)

        integer :: perm(rank)

        integer :: ident(rank)

        integer :: stride_a(rank)
        integer :: stride_a_inner
        integer :: size_outer, size_inner
        integer :: offset_a, offset_b

        integer :: i, j, j_tmp
        integer :: current_index

        call gen_perm_array(perm_str, perm)

        do i=1, rank
            ident(i) = i
        enddo

        if (all(ident == perm)) then
            B = A
            return
        endif

        stride_a(1) = 1
        do i=2, rank
            stride_a(i) = stride_a(i-1) * mat_shape(i-1)
        enddo

        size_outer = 1
        do i=1, rank
            if (i /= perm(1)) size_outer = size_outer * mat_shape(i)
        enddo

        size_inner = mat_shape(perm(1))
        do j=0, size_outer-1
            offset_a = 0

            j_tmp = j

            do i=2, rank

                current_index = modulo(j_tmp, mat_shape(perm(i)))
                j_tmp = j_tmp / mat_shape(perm(i))
                offset_a = offset_a + (current_index * stride_a(perm(i)))

            enddo

            offset_b = j * size_inner
            stride_a_inner = stride_a(perm(1))

            do i=0, size_inner-1
                b(offset_b + i + 1) = a(offset_a + (i * stride_a_inner) + 1)
            enddo


        enddo

    end subroutine reorder_stripe

    subroutine gen_perm_array(perm_str, perm_array)

        character(len=*), intent(in) :: perm_str
        integer, intent(in out) :: perm_array(:)

        integer :: perm_rank
        integer :: i, idx

        perm_rank = len_trim(perm_str)

        do i=1, perm_rank
            idx = iachar(perm_str(i:i)) - 48
            perm_array(i) = idx
        enddo

    end subroutine gen_perm_array

    subroutine permsign(sgn,p)

        integer, intent(in) :: p(:)
        real, intent(out) :: sgn
        integer :: k, ct, L, n
        integer, allocatable :: visited(:)

        n = size(p)
        allocate(visited(1:n))

        do k = 1,n
           visited(k) = 0
        end do

        sgn = 1.0
        do k = 1,n
          if (visited(k) == 0) then
             ct = k
             L = 0
             do while (visited(ct) == 0) 
                L = L + 1
                visited(ct) = 1
                ct = p(ct)
             end do 
             if (mod(L,2) == 0) then
                sgn = -1.0 * sgn
             end if
          end if
        end do

        deallocate(visited)

    end subroutine permsign

    subroutine reorder_max_coincidence(D1,D2,sgn,k)

        use sort_module, only: get_intersection_sorted
        use sort_module, only: argsort_int

        integer, intent(inout) :: D1(:), D2(:)
        real, intent(out) :: sgn
        integer, intent(out) :: k
        integer, allocatable :: inter(:), idx1(:), idx2(:), idx3(:), idx4(:), &
                                perm1(:), perm2(:)
        integer :: m, n, i, cnt
        real :: sgn1, sgn2

        m = size(D1)
        n = size(D2)

        allocate(perm1(m),perm2(n))
        
        ! get the intersection of two sorted arrays D1 and D2
        call get_intersection_sorted(D1,D2,inter,idx1,idx2)

        ! length of intersection
        k = size(inter)

        allocate(idx3(m-k),idx4(n-k))

        ! get indices not intersected in D1
        cnt = 0
        do i = 1,m
           if (any(idx1==i)) then 
                   cycle
           else
                   cnt = cnt + 1
                   idx3(cnt) = i
           end if
        end do
        ! assemble permutation of D1
        perm1(1:k) = idx1
        perm1(k+1:m) = idx3
        ! get permuted D1 and sign of permutation
        D1 = D1(perm1)
        call permsign(sgn1,perm1)


        ! get indices not intersected in D2
        cnt = 0
        do i = 1,n
           if (any(idx2==i)) then
                   cycle
           else
                   cnt = cnt + 1
                   idx4(cnt) = i
           end if
        end do
        ! assemble permutation of D2
        perm2(1:k) = idx2
        perm2(k+1:n) = idx4
        ! get permuted D2 and sign of permutation
        D2 = D2(perm2)
        call permsign(sgn2,perm2)

        ! overall sign of permutations
        sgn = sgn1 * sgn2

        deallocate(inter,idx1,idx2,idx3,idx4,perm1,perm2)

   end subroutine reorder_max_coincidence
        


                


end module permutils
