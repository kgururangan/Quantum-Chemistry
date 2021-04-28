module print_module
        implicit none
        contains
                 subroutine print_matrix(A)

                         real, intent(in) :: A(:,:)
                         integer :: N, i
                         N = size(A,1)
                         do i = 1,N
                            write(*,fmt=*) A(i,:)
                         end do
                 end subroutine print_matrix

end module print_module
