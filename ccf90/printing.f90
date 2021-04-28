module printing

    implicit none

    contains


            subroutine print_header()

                ! Print ccq information. This includes information about
                ! the compilation, the host, authors, etc.


                character(len=255) :: hostname
                character(len=255) :: cmd
                character(len=255) :: cwd
                character(len=255) :: user
                character(len=30) :: date


                write(*,'(a)') 'ccf90 Coupled-Cluster Program'
                write(*,'(a/)') '==========================='

                write(*, '(2x,a)') 'ccf90 repository online: '
                write(*, '(2x,a/)') "Program written in Piecuch's group at MSU: <https://www2.chemistry.msu.edu/faculty/piecuch/>"

                write(*, '(2x,a)') "Authors:"
                write(*, '(4x,a)') "Karthik Gururangan"


            end subroutine print_header

            subroutine print_calc_params(sys)

                ! Print calculation parameters. This includes
                ! molecular system data, CC settings, and other
                ! parameters and configurations related to the QM
                ! calculations.

                ! In:
                !   sys: molecular system data
                !   run: runtime configuration data
                !   cc: CC data, including vectors

                use system_types, only: sys_t

                type(sys_t), intent(in) :: sys

                write(*,'(a)') ''
                write(*,'(a)') 'System information'
                write(*,'(a)') '------------------'
                write(*,'(2x,a27,2x,i16)') 'No. correlated electrons', sys%Nelec
                write(*,'(2x,a27,2x,i16)') 'Frozen orbitals', sys%Nfroz
                write(*,'(2x,a27,2x,i16)') 'Occupied orbitals (alpha)',sys%Nocc_a
                write(*,'(2x,a27,2x,i16)') 'Unoccupied orbitals (alpha)', sys%Nunocc_a
                write(*,'(2x,a27,2x,i16)') 'Occupied orbitals (beta)', sys%Nocc_b
                write(*,'(2x,a27,2x,i16)') 'Unoccupied orbitals (beta)', sys%Nunocc_b
                write(*,'(2x,a27,2x,i16/)') 'Total orbitals', sys%Norb+sys%Nfroz

                !write(*,'(2x,a)') 'Active space partitioning:'
                !write(*,'(2x,a27,2x,i16)') 'Occupied (alpha)', sys%act_occ_a
                !write(*,'(2x,a27,2x,i16)') 'Occupied (beta)', sys%act_occ_b
                !write(*,'(2x,a27,2x,i16)') 'Unoccupied (alpha)', sys%act_unocc_a
                !write(*,'(2x,a27,2x,i16)') 'Unoccupied (beta)', sys%act_unocc_b


                !write(io,'(/a)') 'CC settings'
                !write(io,'(a)')  '-----------'
                !      write(io,'(2x,a27,2x,i16)') 'Number of excited states', nroot
                !write(io,'(2x,a27,2x,a16)') 'Calculation type', trim(run%calc_type)
                !write(io,'(2x,a27,2x,es16.2)') 'Convergence tolerance', run%tol
                !write(io,'(2x,a27,2x,i16)') 'Max. iterations', run%max_iter
                !write(io,'(2x,a27,2x,i16)') 'DIIS space', run%diis_space
                !write(io,'(2x,a27,2x,f16.4)') 'Shift energy', run%shift
                !write(io,'(2x,a27,2x,l16)') 'Restart', run%restart
                !write(io,'(2x,a27,2x,l16/)') 'Restricted CC', run%rhf

                !write(io,'(2x,a27,2x,i16)') 'Active triples indices', run%act_ind_t
                !write(io,'(2x,a27,2x,i16/)') 'Active quadruples indices', run%act_ind_q

                !if (run%ext_cor) then
                !    write(io,'(2x,a)') 'External correction parameters:'
                !    write(io,'(2x,a27,2x,l16/)') 'Use singles and doubles', run%ext_cor_sd
                !endif


                !write(io,'(2x,a)') 'ACC parameters:'
                !write(io,'(2x,a27,2x,5f6.2)') 'T2^2 -> T2 =', cc%acc%t2t2_t2
                !write(io,'(2x,a27,2x,2f6.2)') 'T3 -> T2 =', cc%acc%t3_t2
                !write(io,'(2x,a27,2x,4f6.2)') 'T1*T3 -> T2 =', cc%acc%t1t3_t2
                !write(io,'(2x,a27,2x,3f6.2)') 'T2^2 -> T3 =', cc%acc%t2t2_t3
                !write(io,'(2x,a27,2x,5f6.2)') 'T2*T3 -> T3 =', cc%acc%t2t3_t3

                write(*,'(/a)') 'Starting energies (Eh)'
                write(*,'(a)')  '----------------------'
                write(*,'(2x,a27,2x,f16.10)') 'Nuclear repulsion', sys%Vnuc
                write(*,'(2x,a27,2x,f16.10/)') 'Reference (HF)', sys%Escf

                !call flush(io)

            end subroutine print_calc_params

end module printing

