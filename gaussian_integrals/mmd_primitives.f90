module mmd_primitives

      use constants
      use gaussian_integral_utils, only : fact, fact2, boys_v2,&
                                           gaussian_product,&
                                           gaussian_normfactor 

      implicit none

      contains

              recursive function Efcn(la, lb, t, dist_AB, alpha, beta) result(val)

                      integer, intent(in) :: la, lb, t
                      real, intent(in) :: alpha, beta, dist_AB
                      real :: val, p, q

                      p = alpha + beta
                      q = alpha*beta/p

                      if (t < 0 .OR. t > la + lb) then
                              val = 0.0
                      elseif (la == lb .AND. la == 0 .AND. la == t) then
                              val = exp(-q*dist_AB*dist_AB)
                      elseif (lb == 0) then
                              val = 1.0/(2.0*p)*Efcn(la-1,lb,t-1,dist_AB,alpha,beta) -&
                                    (q*dist_AB/alpha)*Efcn(la-1,lb,t,dist_AB,alpha,beta) +&
                                    (t+1)*Efcn(la-1,lb,t+1,dist_AB,alpha,beta)
                      else
                             val = 1.0/(2.0*p)*Efcn(la,lb-1,t-1,dist_AB,alpha,beta) +&
                                   (q*dist_AB/beta)*Efcn(la,lb-1,t,dist_AB,alpha,beta) +&
                                   (t+1)*Efcn(la,lb-1,t+1,dist_AB,alpha,beta)
                      endif

              end function Efcn

              recursive function Rfcn(t,u,v,n,p,PCx,PCy,PCz,RPC) result(val)
                      use gaussian_integral_utils, only: boys_v2

                      integer, intent(in) :: t, u, v, n
                      real, intent(in) :: p, PCx, PCy, PCz, RPC
                      real :: val 
                      
                      val = 0.0
                      if (t == u .AND. u == v .AND. v == 0) then
                              val = val + (-2.0*p)**real(n)*boys_v2(n,p*RPC*RPC)
                      elseif (t == u .AND. u == 0) then
                              if (v > 1) then
                                      val = val + real(v-1)*Rfcn(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC)
                              endif
                              val = val + PCz*Rfcn(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)
                      elseif (t == 0) then
                              if (u > 1) then
                                      val = val + real(u-1)*Rfcn(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC)
                              endif
                              val = val + PCy*Rfcn(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)
                      else
                              if (t > 1) then
                                      val = val +real(t-1)*Rfcn(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC)
                              endif
                              val = val + PCx*Rfcn(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)
                      endif
               end function Rfcn


                function overlap(a,lmn1,rA,b,lmn2,rB) result(val)
                        real, intent(in) :: a, b
                        integer, dimension(1:3), intent(in) :: lmn1, lmn2
                        real, dimension(1:3), intent(in) :: rA, rB
                        real :: Sx, Sy, Sz, val

                        Sx = Efcn(lmn1(1),lmn2(1),0,rA(1)-rB(1),a,b)
                        Sy = Efcn(lmn1(2),lmn2(2),0,rA(2)-rB(2),a,b)
                        Sz = Efcn(lmn1(3),lmn2(3),0,rA(3)-rB(3),a,b)

                        val = Sx*Sy*Sz*(pi/(a+b))**(1.5)
                end function overlap

                function kinetic(a,lmn1,rA,b,lmn2,rB) result(val)
                        real, intent(in) :: a, b
                        integer, dimension(1:3), intent(in) :: lmn1, lmn2
                        real, dimension(1:3), intent(in) :: rA, rB
                        real :: val, term0, term1, term2, c1, c2, c3


                        term0 = b*(2.0*sum(real(lmn2))+3.0)*overlap(a,lmn1,rA,b,lmn2,rB)
                        term1 = -2.0*(b**2.0)*(overlap(a,lmn1,rA,b,[lmn2(1)+2,lmn2(2),lmn2(3)],rB)&
                                             +overlap(a,lmn1,rA,b,[lmn2(1),lmn2(2)+2,lmn2(3)],rB)&
                                             +overlap(a,lmn1,rA,b,[lmn2(1),lmn2(2),lmn2(3)+2],rB))
                        c1 = real(lmn2(1))*real(lmn2(1)-1)
                        c2 = real(lmn2(2))*real(lmn2(2)-1)
                        c3 = real(lmn2(3))*real(lmn2(3)-1)             
                        term2 = -0.5*(c1*overlap(a,lmn1,rA,b,[lmn2(1)-2,lmn2(2),lmn2(3)],rB)+&
                                      c2*overlap(a,lmn1,rA,b,[lmn2(1),lmn2(2)-2,lmn2(3)],rB)+&
                                      c3*overlap(a,lmn1,rA,b,[lmn2(1),lmn2(2),lmn2(3)-2],rB))
                        val = term0 + term1 + term2
                                           
                end function kinetic

                function nuclear_attraction(a,lmn1,rA,b,lmn2,rB,rC) result(val)

                        real, intent(in) :: a, b
                        integer, dimension(1:3), intent(in) :: lmn1, lmn2
                        real, dimension(1:3), intent(in) :: rA, rB, rC
                        real :: val, p, const, RPC
                        real, dimension(1:3) :: rP
                        integer :: t, u, v

                        call gaussian_product(p, rP, const, a, rA, b, rB)
                        RPC = sqrt(sum((rP-rC)**2.0))
                        val = 0.0

                        do t = 0,lmn1(1)+lmn2(1)
                                do u = 0,lmn1(2)+lmn2(2)
                                        do v = 0,lmn1(3)+lmn2(3)
                                                val = val + Efcn(lmn1(1),lmn2(1),t,rA(1)-rB(1),a,b)*&
                                                            Efcn(lmn1(2),lmn2(2),u,rA(2)-rB(2),a,b)*&
                                                            Efcn(lmn1(3),lmn2(3),v,rA(3)-rB(3),a,b)*&
                                                            Rfcn(t,u,v,0,p,rP(1)-rC(1),rP(2)-rC(2),rP(3)-rC(3),RPC)
                                        enddo
                                enddo
                        enddo

                        val = val*2.0*pi/p

                end function nuclear_attraction

                function electron_repulsion(a,lmn1,rA,b,lmn2,rB,c,lmn3,rC,d,lmn4,rD) result(val)

                        real, intent(in) :: a, b, c, d
                        real, dimension(1:3), intent(in) :: rA, rB, rC, rD
                        integer, dimension(1:3), intent(in) :: lmn1, lmn2, lmn3, lmn4
                        real :: val, p, q, RPQ, alpha, const_p, const_q
                        real, dimension(1:3) :: rP, rQ
                        integer :: t, u, v, tau, nu, phi,&
                                   l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4

                        l1 = lmn1(1)
                        m1 = lmn1(2)
                        n1 = lmn1(3)
                        l2 = lmn2(1)
                        m2 = lmn2(2)
                        n2 = lmn2(3)
                        l3 = lmn3(1)
                        m3 = lmn3(2)
                        n3 = lmn3(3)
                        l4 = lmn4(1)
                        m4 = lmn4(2)
                        n4 = lmn4(3)

                        call gaussian_product(p, rP, const_p, a, rA, b, rB)
                        call gaussian_product(q, rQ, const_q, c, rC, d, rD)
                        RPQ = sqrt(sum((rP-rQ)**2.0))

                        alpha = p*q/(p+q)

                        val = 0.0
                        do t = 0,l1+l2
                          do u = 0,m1+m2
                            do v = 0,n1+n2
                              do tau = 0,l3+l4
                                do nu = 0,m3+m4
                                  do phi = 0,n3+n4
                                    val = val + (-1.0)**(tau+nu+phi)*&
                                    Efcn(l1,l2,t,rA(1)-rB(1),a,b)*&
                                    Efcn(m1,m2,u,rA(2)-rB(2),a,b)*&
                                    Efcn(n1,n2,v,rA(3)-rB(3),a,b)*&
                                    Efcn(l3,l4,tau,rC(1)-rD(1),c,d)*&
                                    Efcn(m3,m4,nu,rC(2)-rD(2),c,d)*&
                                    Efcn(n3,n4,phi,rC(3)-rD(3),c,d)*&
                                    Rfcn(t+tau,u+nu,v+phi,0,alpha,rP(1)-rQ(1),rP(2)-rQ(2),rP(3)-rQ(3),RPQ)
                                  enddo
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo

                        val = val * 2.0*(pi**2.5)/(p*q*sqrt(p+q))

               end function electron_repulsion

end module mmd_primitives 

