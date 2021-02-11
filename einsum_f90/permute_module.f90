module permute_module

        implicit none


        contains

                ! ORDER_Y(k) = q where ORDER_X(q) = k

                subroutine permute2(x,x_out,order_x)

                        ! order is the array that tells you how x is permuted to give you x_out
                        ! we need to find order_y which is the array that tells you how x_out 
                        ! is permuted to give you x
                        ! then, reshape is called with x_out = reshape(x,shape=shape_xout,order=order_y)

                        real, intent(in) :: x(:,:,:,:)
                        integer, intent(in) :: order_x(1:4)
                        real, intent(out) :: x_out(:,:,:,:)
                        integer :: shape_x(1:4), shape_xout(1:4), order_y(1:4), i, id

                        do i = 1,4
                                id = findloc(order_x,i,1)
                                order_y(i) = id
                        end do


                        shape_x = shape(x)
                        shape_xout = shape_x(order_x)

                        x_out = reshape(x,shape_xout,order=order_y)

                end subroutine permute2

                subroutine permute3(x,x_out,order)
        
                        real, intent(in) :: x(:,:,:,:,:,:)
                        integer, intent(in) :: order(1:6)
                        real, intent(out) :: x_out(:,:,:,:,:,:)
                        integer :: shape_x(1:6), shape_xout(1:6)

                        shape_x = shape(x)
                        shape_xout = shape_x(order)

                        x_out = reshape(x,shape_xout,order=order)


                end subroutine permute3

                subroutine permute4(x,x_out,order)
        
                        real, intent(in) :: x(:,:,:,:,:,:,:,:)
                        integer, intent(in) :: order(1:8)
                        real, intent(out) :: x_out(:,:,:,:,:,:,:,:)
                        integer :: shape_x(1:8), shape_xout(1:8)
                        
                        shape_x = shape(x)
                        shape_xout = shape_x(order)

                        x_out = reshape(x,shape_xout,order=order)


                end subroutine permute4
 end module permute_module
