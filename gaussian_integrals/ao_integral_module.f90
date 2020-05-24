module ao_integral_module

        contains


subroutine get_ao_integrals(Smat, Zmat, VVmat, num_orbs, num_atoms, orbs, atom_coords, Z)

     use mmd_integrals
     use orbital_types

     implicit none

     integer, intent(in) :: num_orbs, num_atoms
     real, dimension(1:num_atoms), intent(in) :: Z
     type(basisfcn_t), dimension(1:num_orbs), intent(in) :: orbs
     real, dimension(1:num_atoms,1:3), intent(in) :: atom_coords
     integer :: i, j, k, l, ij, kl

     ! Initialize output arrays
     real, dimension(1:num_orbs,1:num_orbs), intent(out) :: Smat, Zmat
     real, dimension(1:num_orbs,1:num_orbs,1:num_orbs,1:num_orbs),intent(out) :: VVmat


     do i = 1,num_orbs
        do j = i,num_orbs
           Smat(i,j) = sAB(orbs(i),orbs(j))
           Zmat(i,j) = tAB(orbs(i),orbs(j))
           do k = 1,num_atoms
              Zmat(i,j) = Zmat(i,j) + vAB(orbs(i),orbs(j),atom_coords(k,:),Z(k))
           enddo
           Smat(j,i) = Smat(i,j)
           Zmat(j,i) = Zmat(i,j)
        enddo
     enddo


     do i = 1,num_orbs
        do j = 1,i
           ij = i*(i-1)/2+j
              do k = 1,num_orbs
                 do l = 1,k
                    kl = k*(k-1)/2+l
                       if (ij>=kl) then
                               VVmat(i,j,k,l) = eri(orbs(i),orbs(j),orbs(k),orbs(l))
                               VVmat(j,i,k,l) = VVmat(i,j,k,l)
                               VVmat(i,j,l,k) = VVmat(i,j,k,l)
                               VVmat(j,i,l,k) = VVmat(i,j,k,l)
                               VVmat(k,l,i,j) = VVmat(i,j,k,l)
                               VVmat(l,k,i,j) = VVmat(i,j,k,l)
                               VVmat(k,l,j,i) = VVmat(i,j,k,l)
                               VVmat(l,k,j,i) = VVmat(i,j,k,l)
                        endif
                  enddo
              enddo
          enddo
      enddo


end subroutine get_ao_integrals


end module ao_integral_module
