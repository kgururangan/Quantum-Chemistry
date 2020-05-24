module gaussian_integral_utils

      use constants

      implicit none
      
     ! real(16), parameter :: pi = 4 * atan (1.0_16)

      contains

              function fact(n) result(val)
                      integer, intent(in) :: n
                      real :: val
                      val = gamma(real(n+1))
              end function fact

              function nchoosek(n,k) result(val)
                      integer, intent(in) :: n, k
                      real :: val
                      val = fact(n)/(fact(k)*fact(n-k))
              end function nchoosek

              function fact2(n) result(val)
                      integer, intent(in) :: n
                      real :: val
                      integer :: i

                      val = 1.0
                      if (n .le. 0) then
                         val = 1.0
                      elseif (mod(n,2) .eq. 0) then
                         do i = 2,n,2
                                val = val*real(i)
                         enddo
                      else
                         do i = 1,n,2
                                val = val*real(i)
                         enddo
                      end if
              end function fact2

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

              function boys_v2(n,T) result(val)
                      use special_functions, only: hyp1F1

                      integer, intent(in) :: n
                      real, intent(in) :: T
                      real :: val

                      val = real(hyp1F1(real(n+0.5,8), real(n+1.5,8), -real(T,8))&
                                 /real(2.0*real(n,4)+1.0,8),4)

              end function boys_v2
               
              subroutine gaussian_product(c_exp,rC,C,alpha,rA,beta,rB)

                     real, intent(out) :: c_exp, C
                     real, dimension(1:3), intent(out) :: rC
                     real, intent(in) :: alpha, beta
                     real, dimension(1:3), intent(in) :: rA, rB


                     c_exp = alpha + beta
                     rC = (alpha*rA + beta*rB)/c_exp
                     C = exp(-alpha*beta/c*sum( (rA-rB)**2 ))

              end subroutine gaussian_product

              function gaussian_normfactor(shell, alpha) result(val)

                      integer, dimension(1:3), intent(in) :: shell
                      real, intent(in) :: alpha
                      real :: val
                      integer :: l, m, n

                      l = shell(1)
                      m = shell(2)
                      n = shell(3)

                      val = sqrt( (2.0*alpha/pi)**(3.0/2.0) * (4.0*alpha)**(real(l+m+n))/&
                                                           (fact2(2*l-1)*fact2(2*m-1)*fact2(2*n-1)) )

              end function gaussian_normfactor
              
              function gaussian_coeff(l,m,xa,xb,k) result(val)
                        real, intent(in) :: xa, xb
                        integer, intent(in) :: l, m, k
                        integer :: i, j
                        real :: val
                        
                        val = 0.0
                        do i = 0,l
                                do j = 0,m
                                        if (i+j == k) then
                                                val = val + nchoosek(l,i)*nchoosek(m,j)*&
                                                      xa**(real(l-i))*xB**(real(m-j))
                                        else
                                                cycle
                                        endif
                                enddo
                        enddo
              end function gaussian_coeff


end module gaussian_integral_utils
