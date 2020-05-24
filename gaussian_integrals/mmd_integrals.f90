module mmd_integrals

      use mmd_primitives, only: overlap, kinetic, nuclear_attraction, electron_repulsion
      use orbital_types, only: basisfcn_t
      use gaussian_integral_utils, only: gaussian_normfactor

      implicit none

      contains


              function sAB(orbA, orbB) result(val)

                      type(basisfcn_t), intent(in) :: orbA, orbB
                      integer :: nG_a, nG_b, i, j
                      real :: val

                      nG_a = size(orbA%exps)
                      nG_b = size(orbB%exps)

                      val = 0.0
                      do i = 1,nG_a
                           do j = 1,nG_b
                                val = val + gaussian_normfactor(orbA%shell,orbA%exps(i))*&
                                            gaussian_normfactor(orbB%shell,orbB%exps(j))*&
                                            orbA%coeff(i)*orbB%coeff(j)*&
                                            overlap(orbA%exps(i),orbA%shell,orbA%origin,&
                                                    orbB%exps(j),orbB%shell,orbB%origin)
                                
                           enddo
                      enddo

              end function sAB

              function tAB(orbA, orbB) result(val)

                      type(basisfcn_t), intent(in) :: orbA, orbB
                      integer :: nG_a, nG_b, i, j
                      real :: val

                      nG_a = size(orbA%exps)
                      nG_b = size(orbB%exps)

                      val = 0.0
                      do i = 1,nG_a
                           do j = 1,nG_b
                                val = val + gaussian_normfactor(orbA%shell,orbA%exps(i))*&
                                            gaussian_normfactor(orbB%shell,orbB%exps(j))*&
                                            orbA%coeff(i)*orbB%coeff(j)*&
                                            kinetic(orbA%exps(i),orbA%shell,orbA%origin,&
                                                    orbB%exps(j),orbB%shell,orbB%origin)
                                
                           enddo
                      enddo

              end function tAB

              
              function vAB(orbA, orbB, rC, Z) result(val)

                      type(basisfcn_t), intent(in) :: orbA, orbB
                      real, dimension(1:3), intent(in) :: rC
                      real, intent(in) :: Z
                      integer :: nG_a, nG_b, i, j
                      real :: val

                      nG_a = size(orbA%exps)
                      nG_b = size(orbB%exps)

                      val = 0.0
                      do i = 1,nG_a
                           do j = 1,nG_b
                                val = val + gaussian_normfactor(orbA%shell,orbA%exps(i))*&
                                            gaussian_normfactor(orbB%shell,orbB%exps(j))*&
                                            orbA%coeff(i)*orbB%coeff(j)*&
                                            nuclear_attraction(orbA%exps(i),orbA%shell,orbA%origin,&
                                                    orbB%exps(j),orbB%shell,orbB%origin,rC)
                                
                           enddo
                      enddo

                      val = val*(-Z)

              end function vAB


              function eri(orbA, orbB, orbC, orbD) result(val)

                      type(basisfcn_t), intent(in) :: orbA, orbB, orbC, orbD
                      integer :: nG_a, nG_b, nG_c, nG_d, i, j, k, l
                      real :: val

                      nG_a = size(orbA%exps)
                      nG_b = size(orbB%exps)
                      nG_c = size(orbC%exps)
                      nG_d = size(orbD%exps)

                      val = 0.0
                      do i = 1,nG_a
                           do j = 1,nG_b
                                do k = 1,nG_c
                                        do l = 1,nG_d
                                                val = val + gaussian_normfactor(orbA%shell,orbA%exps(i))*&
                                                            gaussian_normfactor(orbB%shell,orbB%exps(j))*&
                                                            gaussian_normfactor(orbC%shell,orbC%exps(k))*&
                                                            gaussian_normfactor(orbD%shell,orbD%exps(l))*&
                                                            orbA%coeff(i)*orbB%coeff(j)*orbC%coeff(k)*orbD%coeff(l)*&
                                                            electron_repulsion(orbA%exps(i),orbA%shell,orbA%origin,&
                                                                               orbB%exps(j),orbB%shell,orbB%origin,&
                                                                               orbC%exps(k),orbC%shell,orbC%origin,&
                                                                               orbD%exps(l),orbD%shell,orbD%origin)
                                        enddo
                                 enddo        
                           enddo
                      enddo

              end function eri

end module mmd_integrals


