module printing

    implicit none

    contains


            subroutine print_header(io)

                ! Print ccq information. This includes information about
                ! the compilation, the host, authors, etc.


                integer, intent(in) :: io
                character(len=255) :: hostname
                character(len=255) :: cmd
                character(len=255) :: cwd
                character(len=255) :: user
                character(len=30) :: date


                write(io,'(a)') 'ccf90 Coupled-Cluster Program'
                write(io,'(a/)') '==========================='

                write(io, '(2x,a)') 'ccf90 repository online: '
                write(io, '(2x,a/)') "Program written in Piecuch's group at MSU: <https://www2.chemistry.msu.edu/faculty/piecuch/>"

                write(io, '(2x,a)') "Authors:"
                write(io, '(4x,a)') "Karthik Gururangan"


            end subroutine print_header

            subroutine print_calc_params(sys,io)

                ! Print calculation parameters. This includes
                ! molecular system data, CC settings, and other
                ! parameters and configurations related to the QM
                ! calculations.

                ! In:
                !   sys: molecular system data
                !   run: runtime configuration data
                !   cc: CC data, including vectors

                use system_types, only: sys_t

                integer, intent(in) :: io
                type(sys_t), intent(in) :: sys

                write(io,'(a)') ''
                write(io,'(a)') 'System information'
                write(io,'(a)') '------------------'
                write(io,'(2x,a27,2x,i16)') 'No. correlated electrons', sys%Nelec
                write(io,'(2x,a27,2x,i16)') 'Frozen orbitals', sys%Nfroz
                write(io,'(2x,a27,2x,i16)') 'Occupied orbitals (alpha)',sys%Nocc_a
                write(io,'(2x,a27,2x,i16)') 'Unoccupied orbitals (alpha)', sys%Nunocc_a
                write(io,'(2x,a27,2x,i16)') 'Occupied orbitals (beta)', sys%Nocc_b
                write(io,'(2x,a27,2x,i16)') 'Unoccupied orbitals (beta)', sys%Nunocc_b
                write(io,'(2x,a27,2x,i16/)') 'Total orbitals', sys%Norb+sys%Nfroz

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

                write(io,'(/a)') 'Starting energies (Eh)'
                write(io,'(a)')  '----------------------'
                write(io,'(2x,a27,2x,f16.10)') 'Nuclear repulsion', sys%Vnuc
                write(io,'(2x,a27,2x,f16.10/)') 'Reference (HF)', sys%Escf

                !call flush(io)

            end subroutine print_calc_params

            subroutine print_L_amps(sys,L,thresh)

                    use system_types, only: sys_t

                    type(sys_t), intent(in) :: sys
                    real, intent(in) :: L(:), thresh
                    integer :: a, b, i, j, n1a, n1b, n2a, n2b, n2c,&
                               noa, nob, nua, nub
                    real, allocatable :: l1a(:,:), l1b(:,:), l2a(:,:,:,:), l2b(:,:,:,:), l2c(:,:,:,:)

                    n1a = sys%Nocc_a * sys%Nunocc_a
                    n1b = sys%Nocc_b * sys%Nunocc_b
                    n2a = sys%Nocc_a**2 * sys%Nunocc_a**2
                    n2b = sys%Nocc_b * sys%Nocc_a * sys%Nunocc_b * sys%Nunocc_a
                    n2c = sys%Nocc_b**2 * sys%Nunocc_b**2

                    noa = sys%Nocc_a
                    nob = sys%Nocc_b
                    nua = sys%Nunocc_a
                    nub = sys%Nunocc_b

                    allocate(l1a(nua,noa),l1b(nub,nob),l2a(nua,nua,nob,nob),l2b(nua,nub,noa,nob),l2c(nub,nub,nob,nob))

                    l1a = reshape(L(1:n1a),(/nua,noa/))
                    l1b = reshape(L(n1a+1:n1a+n1b),(/nub,nob/))
                    l2a = reshape(L(n1a+n1b+1:n1a+n1b+n2a),(/nua,nua,noa,noa/))
                    l2b = reshape(L(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b),(/nua,nub,noa,nob/))
                    l2c = reshape(L(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c),(/nub,nub,nob,nob/))

                    do a = 1,sys%Nunocc_a
                       do i = 1,sys%Nocc_a
                          if (abs(l1a(a,i)) >= thresh) then
                                  print*,'l1a(',a,i,') = ',l1a(a,i)
                          end if
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_a
                       do i = 1,sys%Nocc_a
                          if (abs(l1b(a,i)) >= thresh) then
                                  print*,'l1b(',a,i,') = ',l1b(a,i)
                          end if
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_a
                       do b = 1,sys%Nunocc_a
                          do i = 1,sys%Nocc_a
                             do j = 1,sys%Nocc_a
                                if (abs(l2a(a,b,i,j)) >= thresh) then
                                     print*,'l2a(',a,b,i,j,') = ',l2a(a,b,i,j)
                                end if
                             end do
                          end do
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_a
                       do b = 1,sys%Nunocc_b
                          do i = 1,sys%Nocc_a
                             do j = 1,sys%Nocc_b
                                if (abs(l2b(a,b,i,j)) >= thresh) then
                                     print*,'l2b(',a,b,i,j,') = ',l2b(a,b,i,j)
                                end if
                             end do
                          end do
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_b
                       do b = 1,sys%Nunocc_b
                          do i = 1,sys%Nocc_b
                             do j = 1,sys%Nocc_b
                                if (abs(l2c(a,b,i,j)) >= thresh) then
                                     print*,'l2c(',a,b,i,j,') = ',l2c(a,b,i,j)
                                end if
                             end do
                          end do
                       end do
                    end do

            end subroutine print_L_amps

            subroutine print_R_amps(sys,R,thresh)

                    use system_types, only: sys_t

                    type(sys_t), intent(in) :: sys
                    real, intent(in) :: R(:), thresh
                    integer :: a, b, i, j, n1a, n1b, n2a, n2b, n2c,&
                               noa, nob, nua, nub
                    real, allocatable :: r1a(:,:), r1b(:,:), r2a(:,:,:,:), r2b(:,:,:,:), r2c(:,:,:,:)

                    n1a = sys%Nocc_a * sys%Nunocc_a
                    n1b = sys%Nocc_b * sys%Nunocc_b
                    n2a = sys%Nocc_a**2 * sys%Nunocc_a**2
                    n2b = sys%Nocc_b * sys%Nocc_a * sys%Nunocc_b * sys%Nunocc_a
                    n2c = sys%Nocc_b**2 * sys%Nunocc_b**2

                    noa = sys%Nocc_a
                    nob = sys%Nocc_b
                    nua = sys%Nunocc_a
                    nub = sys%Nunocc_b

                    allocate(r1a(nua,noa),r1b(nub,nob),r2a(nua,nua,nob,nob),r2b(nua,nub,noa,nob),r2c(nub,nub,nob,nob))

                    r1a = reshape(R(1:n1a),(/nua,noa/))
                    r1b = reshape(R(n1a+1:n1a+n1b),(/nub,nob/))
                    r2a = reshape(R(n1a+n1b+1:n1a+n1b+n2a),(/nua,nua,noa,noa/))
                    r2b = reshape(R(n1a+n1b+n2a+1:n1a+n1b+n2a+n2b),(/nua,nub,noa,nob/))
                    r2c = reshape(R(n1a+n1b+n2a+n2b+1:n1a+n1b+n2a+n2b+n2c),(/nub,nub,nob,nob/))

                    do a = 1,sys%Nunocc_a
                       do i = 1,sys%Nocc_a
                          if (abs(r1a(a,i)) >= thresh) then
                                  print*,'r1a(',a,i,') = ',r1a(a,i)
                          end if
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_a
                       do i = 1,sys%Nocc_a
                          if (abs(r1b(a,i)) >= thresh) then
                                  print*,'r1b(',a,i,') = ',r1b(a,i)
                          end if
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_a
                       do b = 1,sys%Nunocc_a
                          do i = 1,sys%Nocc_a
                             do j = 1,sys%Nocc_a
                                if (abs(r2a(a,b,i,j)) >= thresh) then
                                     print*,'r2a(',a,b,i,j,') = ',r2a(a,b,i,j)
                                end if
                             end do
                          end do
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_a
                       do b = 1,sys%Nunocc_b
                          do i = 1,sys%Nocc_a
                             do j = 1,sys%Nocc_b
                                if (abs(r2b(a,b,i,j)) >= thresh) then
                                     print*,'r2b(',a,b,i,j,') = ',r2b(a,b,i,j)
                                end if
                             end do
                          end do
                       end do
                    end do
                    print*,''
                    do a = 1,sys%Nunocc_b
                       do b = 1,sys%Nunocc_b
                          do i = 1,sys%Nocc_b
                             do j = 1,sys%Nocc_b
                                if (abs(r2c(a,b,i,j)) >= thresh) then
                                     print*,'r2c(',a,b,i,j,') = ',r2c(a,b,i,j)
                                end if
                             end do
                          end do
                       end do
                    end do

            end subroutine print_R_amps

end module printing

