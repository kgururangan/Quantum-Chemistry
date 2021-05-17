module slater

        use excitation_types, only: excit_t

        implicit none


        contains


                subroutine create_det(occ,N_int,I)

                        integer, intent(in) :: occ(:), N_int
                        integer(kind=8), intent(out) :: I(N_int,2)
                        integer :: j, tmp, ct, x
                        integer(kind=8), parameter :: TWO = 2

                        I = 0

                        do j = 1,size(occ)
                           if (mod(occ(j),2) == 0) then
                                   x = occ(j)/2
                                   ct = floor(real(x/64))+1
                                   tmp = x - (ct-1)*64
                                   I(ct,2) = I(ct,2) + TWO**(tmp-1)
                           else
                                   x = (occ(j)+1)/2
                                   ct = floor(real(x/64))+1
                                   tmp = x - (ct-1)*64
                                   I(ct,1) = I(ct,1) + TWO**(tmp-1)
                           end if
                        end do

                end subroutine create_det

                subroutine print_excitation(exc,io)

                        integer, intent(in) :: io
                        type(excit_t), intent(in) :: exc
                        integer :: i

                        write(io,fmt=*) 'EXCITATION DEGREE = ',exc%nexcit,'PHASE = ',exc%phase

                        do i = 1,exc%nexcit_alpha
                           write(io,fmt=*) exc%from_a(i),'A->',exc%to_a(i),'A'
                        end do
                        do i = 1,exc%nexcit_beta
                           write(io,fmt=*) exc%from_b(i),'B->',exc%to_b(i),'B'
                        end do

                end subroutine print_excitation
                        


                subroutine get_exc_degree(det1,det2,N_int,n_excitations)

                        integer, intent(in) :: N_int
                        integer(kind=8), intent(in) :: det1(N_int,2), det2(N_int,2)
                        integer, intent(out) :: n_excitations
                        integer :: l

                        n_excitations = &
                                popcnt(ieor(det1(1,1),det2(1,1))) + &
                                popcnt(ieor(det1(1,2),det2(1,2)))
                        do l = 2,N_int
                           n_excitations = n_excitations + &
                                   popcnt(ieor( det1(l,1), det2(l,1)) ) +&
                                   popcnt(ieor( det1(l,2), det2(l,2)) )
                        end do
                        n_excitations = ishft(n_excitations,-1)

               end subroutine get_exc_degree

               subroutine get_excitation(det1,det2,N_int,exc)

                       integer, intent(in) :: N_int
                       integer(kind=8), intent(in) :: det1(N_int,2), det2(N_int,2)
                       type(excit_t), intent(out) :: exc

                       integer :: excmat(0:2,2,2), degree
                       real :: phase

                       call get_exc_degree(det1,det2,N_int,degree)

                       exc%nexcit = degree

                       select case (degree)

                               case (3:)
                                       degree = -1
                                       return
                               case (2)
                                       call get_double_excitation(det1,det2,excmat,phase,N_int)

                                       exc%nexcit_alpha = excmat(0,1,1)
                                       exc%nexcit_beta = excmat(0,1,2)
                                       exc%from_a = (/excmat(1,1,1), excmat(2,1,1), 0, 0/)
                                       exc%from_b = (/excmat(1,1,2), excmat(2,1,2), 0, 0/)
                                       exc%to_a = (/excmat(1,2,1), excmat(2,2,1), 0, 0/)
                                       exc%to_b = (/excmat(1,2,2), excmat(2,2,2), 0, 0/)
                                       exc%phase = phase

                               case (1)
                                       call get_single_excitation(det1,det2,excmat,phase,N_int)

                                       exc%nexcit_alpha = excmat(0,1,1)
                                       exc%nexcit_beta = excmat(0,1,2)
                                       exc%from_a = (/excmat(1,1,1), excmat(2,1,1), 0, 0/)
                                       exc%from_b = (/excmat(1,1,2), excmat(2,1,2), 0, 0/)
                                       exc%to_a = (/excmat(1,2,1), excmat(2,2,1), 0, 0/)
                                       exc%to_b = (/excmat(1,2,2), excmat(2,2,2), 0, 0/)
                                       exc%phase = phase

                               case (0)
                                       return
                       end select

               end subroutine get_excitation

               subroutine get_single_excitation(det1,det2,exc,phase,N_int)

                       integer, intent(in) :: N_int
                       integer(kind=8), intent(in) :: det1(N_int,2), det2(N_int,2)
                       integer, intent(out) :: exc(0:2,2,2)
                       real, intent(out) :: phase
                       integer :: tz, l, ispin, ishift, nperm, i, j, k, m, n, high, low
                       integer(kind=8) :: hole, particle, tmp
                       real, parameter :: phase_dble(0:1) = (/ 1.0, -1.0 /)
                       integer(kind=8), parameter :: ZERO = 0, ONE = 1

                       exc = 0
                       exc(0,1,1) = 0
                       exc(0,2,1) = 0
                       exc(0,1,2) = 0
                       exc(0,2,2) = 0
                       do ispin = 1,2
                          ishift = -63
                          do l = 1,N_int
                             ishift = ishift + 64 
                             if (det1(l,ispin) == det2(l,ispin)) cycle
                             tmp = ieor( det1(l,ispin), det2(l,ispin) )
                             particle = iand(tmp, det2(l,ispin))
                             hole = iand(tmp, det1(l,ispin))
                             if (particle .ne. ZERO) then
                                     tz = trailz(particle)
                                     exc(0,2,ispin) = 1
                                     exc(1,2,ispin) = tz+ishift
                             end if 
                             if (hole .ne. ZERO) then
                                     tz = trailz(hole)
                                     exc(0,1,ispin) = 1
                                     exc(1,1,ispin) = tz+ishift
                             end if

                             if ( iand(exc(0,1,ispin),exc(0,2,ispin)) == 1 ) then
                                     low = min(exc(1,1,ispin),exc(1,2,ispin))
                                     high = max(exc(1,1,ispin),exc(1,2,ispin))
                                     j = ishft(low-1,-6)+1
                                     n = iand(low,63)
                                     k = ishft(high-1,-6)+1
                                     m = iand(high,63)
                                     if (j==k) then
                                             nperm = popcnt(iand(det1(j,ispin), &
                                                     iand( ibset(ZERO,m-1)-ONE, ibclr(-ONE,n)+ONE) ))
                                     else
                                             nperm = popcnt(iand(det1(k,ispin), ibset(ZERO,m-1)-ONE)) + &
                                                     popcnt(iand(det1(j,ispin), ibclr(-ONE,n)+ONE))
                                             do i = j+1,k-1
                                                nperm = nperm + popcnt(det1(i,ispin))
                                             end do
                                     end if
                                     phase = phase_dble(iand(nperm,1))
                                     return
                             end if
                          end do
                       end do

               end subroutine get_single_excitation

               subroutine get_double_excitation(det1,det2,exc,phase,N_int)

                       use sort_module, only: argsort_int

                       integer, intent(in)  :: N_int
                       integer(kind=8), intent(in)  :: det1(N_int,2), det2(N_int,2)
                       integer, intent(out) :: exc(0:2,2,2)
                       real, intent(out) :: phase
                       integer :: l, ispin, idx_hole, idx_particle, ishift
                       integer :: i,j,k,m,n,high, low,a,b,c,d,nperm,tz,nexc
                       integer(kind=8) :: hole, particle, tmp
                       real, parameter :: phase_dble(0:1) = (/ 1.0, -1.0 /)
                       integer(kind=8), parameter :: ZERO = 0, ONE = 1

                       integer :: idx(2)

                       !idx(1:exc(0,1,1)) = argsort_int(exc(1,1:exc(0,1,1),1))
                       !exc(1,1:exc(0,1,1),1) = exc(1,idx(1:exc(0,1,1)),1)


                       !exc = 0
                       exc(0,1,1) = 0
                       exc(0,2,1) = 0
                       exc(0,1,2) = 0
                       exc(0,2,2) = 0
                       nexc=0
                       nperm=0
                       do ispin = 1,2
                        idx_particle = 0
                        idx_hole = 0
                        ishift = -63
                        do l=1,N_int
                         ishift = ishift + 64
                         if (det1(l,ispin) == det2(l,ispin))  then
                           cycle
                         end if
                         tmp = ieor( det1(l,ispin), det2(l,ispin) )
                         particle = iand(tmp, det2(l,ispin))
                         hole     = iand(tmp, det1(l,ispin))
                         do while (particle .ne. 0_8)
                           tz = trailz(particle)
                           nexc = nexc+1
                           idx_particle = idx_particle + 1
                           exc(0,2,ispin) = exc(0,2,ispin) + 1
                           exc(idx_particle,2,ispin) = tz+ishift
                           particle = iand(particle,particle-1_8)
                         end do
                         do while (hole .ne. ZERO)
                           tz = trailz(hole)
                           nexc = nexc+1
                           idx_hole = idx_hole + 1
                           exc(0,1,ispin) = exc(0,1,ispin) + 1
                           exc(idx_hole,1,ispin) = tz+ishift
                           hole = iand(hole,hole-ONE)
                         end do
                         if (nexc == 4) exit
                        end do

                        do i = 1,exc(0,1,ispin)
                           low = min(exc(i,1,ispin),exc(i,2,ispin))
                           high = max(exc(i,1,ispin),exc(i,2,ispin))
                           j = ishft(low-1,-6)+1
                           n = iand(low,63)
                           k = ishft(high-1,-6)+1
                           m = iand(high,63)
                           if (j==k) then
                                   nperm = nperm + popcnt(iand(det1(j,ispin), &
                                           iand( ibset(ZERO,m-1)-ONE, ibclr(-ONE,n)+ONE) ))
                           else
                                   nperm = nperm + popcnt(iand(det1(k,ispin), &
                                                          ibset(ZERO,m-1)-ONE)) &
                                                 + popcnt(iand(det1(j,ispin), &
                                                          ibclr(-ONE,n)+ONE))
                                   do l = j+1,k-1
                                      nperm = nperm + popcnt(det1(l,ispin))
                                   end do
                           end if
                        end do
                        if (exc(0,1,ispin) == 2) then
                                a = min(exc(1,1,ispin), exc(1,2,ispin))
                                b = max(exc(1,1,ispin), exc(1,2,ispin))
                                c = min(exc(2,1,ispin), exc(2,2,ispin))
                                d = max(exc(2,1,ispin), exc(2,2,ispin))
                                if ( (c>a) .and. (c<b) .and. (d>b) ) then
                                        nperm = nperm + 1
                                end if
                                exit
                        end if
                     end do
                     phase = phase_dble(iand(nperm,1))


               end subroutine get_double_excitation



end module slater

