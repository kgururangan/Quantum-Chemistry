module utils_fcn

      implicit none
      
      real(16), parameter :: pi = 4 * atan (1.0_16)

      contains

              function boys(t) result(val)
                      ! Boys function
                      ! F0(t) = sqrt(pi/t)*erf(sqrt(t)) if t not = 0
                      ! F0(0) = 1

                      real, intent(in) :: t
                      real :: val

                      if (t .EQ. 0.0) then
                              val = 1.000
                      else
                              val = 0.5*sqrt(pi/t)*erf(sqrt(t))
                      endif
              end function boys

end module utils_fcn
