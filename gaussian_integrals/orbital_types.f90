module orbital_types

      implicit none

      type basisfcn_t

              integer, dimension(1:3) :: shell
              real, allocatable :: coeff(:)
              real, allocatable :: exps(:)
              real, dimension(1:3) :: origin
      end type basisfcn_t

end module orbital_types
