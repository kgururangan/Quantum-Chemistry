module ao_integrals
        implicit none

        contains 

                function sAB(xiA_arr, coeffA_arr, rA, xiB_arr, coeffB_arr, rB) result(val)

                         use gaussian_integration, only: overlap_integral


                         integer :: i, j, num_gauss_A, num_gauss_B
                         real :: val
                         real, dimension(:) :: xiA_arr, coeffA_arr, xiB_arr, coeffB_arr
                         real, dimension(1:3) :: rA, rB

                         num_gauss_A = size(xiA_arr)
                         num_gauss_B = size(xiB_arr)
                         val = 0.00000
                         do i = 1,num_gauss_A
                                do j = 1,num_gauss_B
                                        val = val + coeffA_arr(i)*coeffB_arr(j)*overlap_integral(xiA_arr(i),rA,xiB_arr(j),rB)
                                enddo
                         enddo

                end function sAB
                
                function tAB(xiA_arr, coeffA_arr, rA, xiB_arr, coeffB_arr, rB) result(val)
                        
                         use gaussian_integration, only: KE_integral


                         integer :: i, j, num_gauss_A, num_gauss_B
                         real :: val
                         real, dimension(:) :: xiA_arr, coeffA_arr, xiB_arr, coeffB_arr
                         real, dimension(1:3) :: rA, rB

                         num_gauss_A = size(xiA_arr)
                         num_gauss_B = size(xiB_arr)
                         val = 0.00000
                         do i = 1,num_gauss_A
                                do j = 1,num_gauss_B
                                        val = val + coeffA_arr(i)*coeffB_arr(j)*KE_integral(xiA_arr(i),rA,xiB_arr(j),rB)
                                enddo
                         enddo

                 end function tAB

                function vAB(xiA_arr, coeffA_arr, rA, xiB_arr, coeffB_arr, rB, R_atom, Z) result(val)
                        
                         use gaussian_integration, only: PE_integral


                         integer :: i, j, k, num_gauss_A, num_gauss_B, num_atoms
                         real :: val
                         real, dimension(:) :: xiA_arr, coeffA_arr, xiB_arr, coeffB_arr, Z
                         real, dimension(1:3) :: rA, rB
                         real, dimension(:,:) :: R_atom

                         num_gauss_A = size(xiA_arr)
                         num_gauss_B = size(xiB_arr)
                         num_atoms = size(Z)

                         val = 0.00000
                         do i = 1,num_gauss_A
                                do j = 1,num_gauss_B
                                        do k = 1,num_atoms
                                                val = val + coeffA_arr(i)*coeffB_arr(j)*&
                                                        PE_integral(xiA_arr(i),rA,xiB_arr(j),rB,R_atom(k,:),Z(k))
                                        enddo
                                enddo
                         enddo

                 end function vAB

                function vABCD(xiA_arr, coeffA_arr, rA, &
                               xiB_arr, coeffB_arr, rB, &
                               xiC_arr, coeffC_arr, rC, &
                               xiD_arr, coeffD_arr, rD) result(val)
                        
                         use gaussian_integration, only: two_body_integral


                         integer :: i, j, k, l, num_gauss_A, num_gauss_B, num_gauss_C, num_gauss_D
                         real :: val
                         real, dimension(:) :: xiA_arr, coeffA_arr, xiB_arr, coeffB_arr, xiC_arr, coeffC_arr, xiD_arr, coeffD_arr
                         real, dimension(1:3) :: rA, rB, rC, rD

                         num_gauss_A = size(xiA_arr)
                         num_gauss_B = size(xiB_arr)
                         num_gauss_C = size(xiC_arr)
                         num_gauss_D = size(xiD_arr)

                         val = 0.00000
                         do i = 1,num_gauss_A
                                do j = 1,num_gauss_B
                                        do k = 1,num_gauss_C
                                                do l = 1,num_gauss_D
                                                        val = val + coeffA_arr(i)*coeffB_arr(j)*coeffC_arr(k)*coeffD_arr(l)*&
                                                     two_body_integral(xiA_arr(i),rA,xiB_arr(j),rB,xiC_arr(k),rC,xiD_arr(l),rD)
                                                enddo
                                        enddo
                                enddo
                         enddo

                 end function vABCD
end module ao_integrals
