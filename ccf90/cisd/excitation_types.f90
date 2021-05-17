module excitation_types

        implicit none


        
        type excit_t

                integer :: nexcit
                integer :: nexcit_alpha
                integer :: nexcit_beta

                integer :: from_a(4), from_b(4), to_a(4), to_b(4)
                !integer :: from_orb(4), to_orb(4)
                real :: phase

        end type excit_t

end module excitation_types
