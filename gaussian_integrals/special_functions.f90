module special_functions

      use constants

      implicit none

      contains

              function hyp1F1(a,b,x) result(hg)
                      
                       real (kind=8) :: a, b, x, a0, a1, hg, hg1, hg2,&
                                        r, r1, r2, rg, sum1, sum2, ta,&
                                        tb, tba, x0, xg, y0, y1
                       integer (kind=4) :: i, j, k, la, m, n, n1

                       a0 = a
                       a1 = a
                       x0 = x
                       hg = 0.0d+00

                       if (b == 0.0d+00 .OR. b == -abs(int(b))) then
                               hg = 1.0d+300
                       elseif (a==0.0d+00 .or. x==0.0d+00) then
                               hg = 1.0d+00
                       elseif (a==-1.0d+00) then
                               hg = 1.0d+00 - x/b
                       elseif (a==b) then
                               hg = exp(x)
                       elseif (a-b==1.0d+00) then
                               hg = (1.0d+00 + x/b)*exp(x)
                       elseif (a==1.0d+00 .AND. b==2.0d+00) then
                               hg = (exp(x)-1.0d+00)/x
                       elseif (a==int(a) .AND. a<0.0d+00) then
                               m = int(-a)
                               r = 1.0d+00
                               hg = 1.0d+00
                               do k = 1,m
                                  r = r*(a+k-1.0d+00)/k/(b+k-1.0d+00)*x
                                  hg = hg + r
                               enddo
                       endif

                       if (hg /= 0.0d+00) then
                              return
                       endif

                       if (x<0.0d+00) then
                                a = b - a
                                a0 = a
                                x = abs(x)
                       endif

                       if (a < 2.0d+00) then
                              n1 = 0
                       endif

                       if (2.0d+00 <= a) then
                              n1 = 1
                              la = int(a)
                              a = a - la - 1.0d+00
                       endif

                       do n = 0,n1
                           if (2.0d+00 <= a0) then
                                  a = a + 1.0d+00
                           endif

                          if (x <= 30.0e+00 + abs(b) .OR. a < 0.0d+00 ) then 
                                  hg = 1.0d+00
                                  rg = 1.0d+00
                                  do j = 1,500
                                        rg = rg*(a+j-1.0d+00)&
                                             /(j*(b+j-1.0d+00))*x
                                        hg = hg + rg
                                        if (abs(rg/hg) < 1.0d-15) then
                                                exit
                                        endif
                                  enddo
                          else
                                  ta = gamma(a)
                                  tb = gamma(b)
                                  xg = b-a
                                  tba = gamma(xg)
                                  sum1 = 1.0d+00
                                  sum2 = 1.0d+00
                                  r1 = 1.0d+00
                                  r2 = 1.0d+00
                                  do i = 1,8
                                        r1 =-r1*(a+i-1.0d+00)*(a-b+i)/(x*i)
                                        r2 =-r2*(b-a+i-1.0d+00)*(a-i)/(x*i)
                                        sum1 = sum1+r1
                                        sum2 = sum2+r2
                                  enddo
                                  hg1 = tb/tba*x**(-a)*cos(pi*a)*sum1
                                  hg2 = tb/ta*exp(x)*x**(a-b)*sum2
                                  hg=hg1+hg2

                           endif
                         enddo

                         if (2.0d+00 <= a0) then
                                 do i=1,la-1
                                        hg=((2.0d+00*a-b+x)*y1+(b-a)*y0)/a
                                        y0=y1
                                        y1=hg
                                        a=a+1.0d+00
                                 enddo
                         endif

                         if (x0<0.0d+00) then
                                 hg=hg*exp(x0)
                         endif

                         a=a1
                         x=x0

                         return

               end function hyp1F1


end module special_functions
