program test

    implicit none

    real(kind=8), allocatable :: a(:, :, :)

    allocate(a(20,4,3))
    call fill_array(a)
    call print_array(a)
    deallocate(a)


contains
    subroutine fill_array(array)
        real(kind=8), allocatable, intent(in out) :: array(..)

        call fill_array_n(array, size(array))

    end subroutine fill_array

    subroutine fill_array_n(array, n)
        integer, intent(in) :: n
        real(kind=8), intent(in out) :: array(n)
        integer :: i

        do i=1, n
            array(i) = i
        end do
    end subroutine fill_array_n

    subroutine print_array(array)
        real(kind=8), allocatable, intent(in) :: array(..)
        !real(kind=8), allocatable :: array2

        !allocate(array2,mold=array)

        print*,'The shape of the assumed-rank array is ',shape(array)
        print*,'The number of elements in assumed-rank array is ',size(array)

        call print_array_n(array, size(array))

    end subroutine print_array

    subroutine print_array_n(array, n)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: array(n)
        integer :: i

        do i=1, n
            print '(i4, f10.4)', i, array(i)
        end do
    end subroutine print_array_n

end program test
