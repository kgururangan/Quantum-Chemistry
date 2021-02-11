module tensor_type

      implicit none

      type tensor_t

              real, allocatable :: vec(:)
              integer, allocatable :: dimensions(:)

      end type tensor_t

end module tensor_type
