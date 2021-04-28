module integral_types

      implicit none

      type e1int_t

              real, allocatable :: oo(:,:)
              real, allocatable :: uu(:,:)
              real, allocatable :: ou(:,:)
              real, allocatable :: uo(:,:)

      end type e1int_t

      !type e2aint_t

      !        real, allocatable :: oooo(:,:,:,:)
      !        real, allocatable :: ooou(:,:,:,:)
      !        real, allocatable :: uooo(:,:,:,:)
      !        real, allocatable :: uoou(:,:,:,:)
      !        real, allocatable :: uuou(:,:,:,:)
      !        real, allocatable :: uuoo(:,:,:,:)
      !        real, allocatable :: oouu(:,:,:,:)
      !        real, allocatable :: uouu(:,:,:,:)
      !        real, allocatable :: uuuu(:,:,:,:)

      !end type e2aint_t

      type e2int_t
              
              real, allocatable :: oooo(:,:,:,:)
              real, allocatable :: ooou(:,:,:,:)
              real, allocatable :: oouo(:,:,:,:)
              real, allocatable :: uooo(:,:,:,:)
              real, allocatable :: ouoo(:,:,:,:)
              real, allocatable :: uoou(:,:,:,:)
              real, allocatable :: ouou(:,:,:,:)
              real, allocatable :: uouo(:,:,:,:)
              real, allocatable :: ouuo(:,:,:,:)
              real, allocatable :: oouu(:,:,:,:)
              real, allocatable :: uuoo(:,:,:,:)
              real, allocatable :: uuou(:,:,:,:)
              real, allocatable :: uuuo(:,:,:,:)
              real, allocatable :: uouu(:,:,:,:)
              real, allocatable :: ouuu(:,:,:,:)
              real, allocatable :: uuuu(:,:,:,:)

      end type e2int_t

end module integral_types
