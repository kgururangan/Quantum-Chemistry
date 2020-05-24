program main


      use orbital_types, only: basisfcn_t
      use ao_integral_module, only: get_ao_integrals

      implicit none

      type(basisfcn_t), dimension(1:4) :: orbs
      real, allocatable :: atom_coords(:,:), Z(:)
      integer :: i, j, k, l, ij, kl,  Norb, Natom
      real, allocatable :: Smat(:,:), Zmat(:,:)
      real, allocatable :: VVmat(:,:,:,:)
      real :: time_start, time_finish


      Natom = 2

      allocate(atom_coords(Natom,3))
      allocate(Z(Natom))

      Z = [1.0, 1.0]

      atom_coords(1,:) = [0.0, 0.0, -0.6991986157742016]
      atom_coords(2,:) = [0.0, 0.0, 0.6991986157742016]

      ! Derived types can be stored in a list
      ! 6-31G basis set for H2
      orbs(1)%origin = atom_coords(1,:)
      orbs(1)%exps = [18.73113696, 2.825394365, 0.6401216923]
      orbs(1)%coeff = [0.0334946043412761, 0.23472695350894282, 0.8137573261310032]
      orbs(1)%shell = [0,0,0]

      orbs(2)%origin = atom_coords(1,:)
      orbs(2)%exps = [0.1612777588]
      orbs(2)%coeff = [1.0]
      orbs(2)%shell = [0,0,0]

      orbs(3)%origin = atom_coords(2,:)
      orbs(3)%exps = [18.73113696, 2.825394365, 0.6401216923]
      orbs(3)%coeff = [0.0334946043412761, 0.23472695350894282, 0.8137573261310032]
      orbs(3)%shell = [0,0,0]

      orbs(4)%origin = atom_coords(2,:)
      orbs(4)%exps = [0.1612777588]
      orbs(4)%coeff = [1.0]
      orbs(4)%shell = [0,0,0]

      Norb = size(orbs)


      ! Starting time of integral routine
      call cpu_time(time_start)
      allocate(Smat(Norb,Norb))
      allocate(Zmat(Norb,Norb))
      allocate(VVmat(Norb,Norb,Norb,Norb))
      call get_ao_integrals(Smat,Zmat,VVmat,Norb,Natom,orbs,atom_coords,Z)
      call cpu_time(time_finish)
      print*,'AO integrals constructed in ',time_finish-time_start,'s'

      do i = 1,Norb
        do j = 1,Norb
            print*,'S(',i,',',j,') =',Smat(i,j)
        enddo
      enddo

      do i = 1,Norb
        do j = 1,Norb
            print*,'Z(',i,',',j,') =',Zmat(i,j)
        enddo
      enddo

      do i = 1,Norb
        do j = 1,Norb
          do k = 1,Norb
            do l = 1,Norb            
                print*,'VV(',i,',',j,',',k,',',l,') = ',VVmat(i,j,k,l)
            enddo
          enddo
        enddo
      enddo

end program main
