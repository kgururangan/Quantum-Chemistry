module utils
        
      implicit none

      contains
              function linspace(x_start, x_end, num_pts) result(xs)

                real, intent(in) :: x_start, x_end
                integer, intent(in) :: num_pts
                real, dimension(1:num_pts) :: xs
                integer :: i
                real :: dx

                dx = (x_end-x_start)/(real(num_pts-1))

                do i = 1,num_pts
                        xs(i) = x_start + real(i-1)*dx
                enddo

               end function linspace

end module utils

