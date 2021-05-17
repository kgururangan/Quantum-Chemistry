module system_types

      implicit none

      type sys_t

              integer :: Norb
              integer :: Nelec
              integer :: Nfroz
              integer :: Nocc_a
              integer :: Nocc_b
              integer :: Nunocc_a
              integer :: Nunocc_b
              real :: Vnuc
              real :: Escf

      end type sys_t


end module system_types
