module integrals

        use integral_types, only: e1int_t, e2int_t
        use system_types, only: sys_t

        implicit none

        contains

                subroutine get_integrals(onebody_fid,twobody_fid,Norb,Nfroz,Nelec,sys,ERI_a,ERI_b,ERI_c,Fock_a,Fock_b)

                        character(*), intent(in) :: onebody_fid, twobody_fid
                        integer, intent(in) :: Norb, Nfroz, Nelec
                        type(e1int_t), intent(out) :: Fock_a, Fock_b
                        type(e2int_t), intent(out) :: ERI_a, ERI_b, ERI_c
                        type(sys_t), intent(out) :: sys
                        real :: e1int(Norb,Norb), e2int(Norb,Norb,Norb,Norb), Vnuc, fA(Norb,Norb), fB(Norb,Norb), &
                                vA(Norb,Norb,Norb,Norb), vB(Norb,Norb,Norb,Norb), vC(Norb,Norb,Norb,Norb), Escf
                        integer :: p, q, i, j, Nocc_a, Nocc_b, Nua, Nub, Noa, Nob
                        integer, allocatable :: ioa(:), iua(:), iob(:), iub(:)


                        write(*,'(/a)') 'Integral Files'
                        write(*,'(a)') '--------------'
                        write(*,'(2x,a)') onebody_fid 
                        call load_onebody(onebody_fid, Norb, e1int)
                        write(*,'(2x,a)') twobody_fid 
                        call load_twobody(twobody_fid, Norb, e2int, Vnuc)

                        if (mod(Nelec,2) == 0) then
                            Nocc_a = Nelec/2
                            Nocc_b = Nelec/2
                        else
                            Nocc_a = (Nelec+1)/2
                            Nocc_b = (Nelec-1)/2
                        end if

                        ! calculate HF energy
                        Escf = Vnuc
                        do i = 1,Nocc_a
                           Escf = Escf + e1int(i,i)
                        end do
                        do i = 1,Nocc_b
                           Escf = Escf + e1int(i,i)
                        end do
                        do i = 1,Nocc_a
                           do j = 1,Nocc_b
                              Escf = Escf + e2int(i,j,i,j)
                           end do
                        end do
                        do i = 1,Nocc_a
                           do j = 1,Nocc_a
                              Escf = Escf + 0.5*(e2int(i,j,i,j) - e2int(i,j,j,i))
                           end do
                        end do
                        do i = 1,Nocc_b
                           do j = 1,Nocc_b
                              Escf = Escf + 0.5*(e2int(i,j,i,j) - e2int(i,j,j,i))
                           end do
                        end do

                        ! get antisymmetrized V matrices
                        call get_V(Norb,e2int,vA,vB,vC)

                        ! calculate fock matrices
                        do p = 1,Norb
                          do q = 1,Norb
                             fA(p,q) = e1int(p,q)
                             fB(p,q) = e1int(p,q)
                             do i = 1,Nocc_a
                                fA(p,q) = fA(p,q) + vA(p,i,q,i)
                                fB(p,q) = fB(p,q) + vB(i,p,i,q)
                             end do
                             do i = 1,Nocc_b
                                fA(p,q) = fA(p,q) + vB(p,i,q,i)
                                fB(p,q) = fB(P,q) + vC(p,i,q,i)
                             end do
                           end do
                         end do

                         ! integral slicing arrays
                         Noa = Nocc_a - Nfroz
                         Nob = Nocc_b - Nfroz
                         Nua = Norb - Nocc_a
                         Nub = Norb - Nocc_b
                         allocate(ioa(Noa),iob(Nob),iua(Nua),iub(Nub))

                         ioa = (/(i, i=Nfroz+1,Nocc_a, 1)/)
                         iob = (/(i, i=Nfroz+1,Nocc_b, 1)/)
                         iua = (/(i, i=Nocc_a+1,Norb, 1)/)
                         iub = (/(i, i=Nocc_b+1,Norb, 1)/)

                         ! slice fA
                         Fock_a%oo = fA(ioa,ioa)  
                         Fock_a%ou = fA(ioa,iua)
                         Fock_a%uo = fA(iua,ioa)
                         Fock_a%uu = fA(iua,iua)
                        
                         ! slice fB
                         Fock_b%oo = fB(iob,iob)
                         Fock_b%ou = fB(iob,iub)
                         Fock_b%uo = fB(iub,iob)
                         Fock_b%uu = fB(iub,iub)   

                         ! slice vA
                         ERI_a%oooo = vA(ioa,ioa,ioa,ioa)
                         ERI_a%ooou = vA(ioa,ioa,ioa,iua)
                         ERI_a%oouo = vA(ioa,ioa,iua,ioa)
                         ERI_a%uooo = vA(iua,ioa,ioa,ioa)
                         ERI_a%ouoo = vA(ioa,iua,ioa,ioa)
                         ERI_a%uoou = vA(iua,ioa,ioa,iua)
                         ERI_a%ouou = vA(ioa,iua,ioa,iua)
                         ERI_a%uouo = vA(iua,ioa,iua,ioa)
                         ERI_a%ouuo = vA(ioa,iua,iua,ioa)
                         ERI_a%uuoo = vA(iua,iua,ioa,ioa)
                         ERI_a%oouu = vA(ioa,ioa,iua,iua)
                         ERI_a%uuou = vA(iua,iua,ioa,iua)
                         ERI_a%uuuo = vA(iua,iua,iua,ioa)
                         ERI_a%uouu = vA(iua,ioa,iua,iua)
                         ERI_a%ouuu = vA(ioa,iua,iua,iua)
                         ERI_a%uuuu = vA(iua,iua,iua,iua)
                              
                         ! slice vC              
                         ERI_c%oooo = vC(iob,iob,iob,iob)
                         ERI_c%ooou = vC(iob,iob,iob,iub)
                         ERI_c%oouo = vC(iob,iob,iub,iob)
                         ERI_c%uooo = vC(iub,iob,iob,iob)
                         ERI_c%ouoo = vC(iob,iub,iob,iob)
                         ERI_c%uoou = vC(iub,iob,iob,iub)
                         ERI_c%ouou = vC(iob,iub,iob,iub)
                         ERI_c%uouo = vC(iub,iob,iub,iob)
                         ERI_c%ouuo = vC(iob,iub,iub,iob)
                         ERI_c%uuoo = vC(iub,iub,iob,iob)
                         ERI_c%oouu = vC(iob,iob,iub,iub)
                         ERI_c%uuou = vC(iub,iub,iob,iub)
                         ERI_c%uuuo = vC(iub,iub,iub,iob)
                         ERI_c%uouu = vC(iub,iob,iub,iub)
                         ERI_c%ouuu = vC(iob,iub,iub,iub)
                         ERI_c%uuuu = vC(iub,iub,iub,iub)

                         ! slice vB
                         ERI_b%oooo = vB(ioa,iob,ioa,iob)
                         ERI_b%ooou = vB(ioa,iob,ioa,iub)
                         ERI_b%oouo = vB(ioa,iob,iua,iob)
                         ERI_b%uooo = vB(iua,iob,ioa,iob)
                         ERI_b%ouoo = vB(ioa,iub,ioa,iob)
                         ERI_b%uoou = vB(iua,iob,ioa,iub)
                         ERI_b%ouou = vB(ioa,iub,ioa,iub)
                         ERI_b%uouo = vB(iua,iob,iua,iob)
                         ERI_b%ouuo = vB(ioa,iub,iua,iob)
                         ERI_b%uuoo = vB(iua,iub,ioa,iob)
                         ERI_b%oouu = vB(ioa,iob,iua,iub)
                         ERI_b%uuou = vB(iua,iub,ioa,iub)
                         ERI_b%uuuo = vB(iua,iub,iua,iob)
                         ERI_b%uouu = vB(iua,iob,iua,iub)
                         ERI_b%ouuu = vB(ioa,iub,iua,iub)
                         ERI_b%uuuu = vB(iua,iub,iua,iub)

                         deallocate(ioa,iua,iob,iub)

                         ! populate correlated system information struct
                         sys%Norb = Norb - Nfroz
                         sys%Nelec = Nelec - 2*Nfroz
                         sys%Nfroz = Nfroz
                         sys%Nocc_a = Noa
                         sys%Nunocc_a = Nua
                         sys%Nocc_b = Nob
                         sys%Nunocc_b = Nub
                         sys%Vnuc = Vnuc
                         sys%Escf = Escf

                end subroutine get_integrals


                subroutine load_onebody(onebody_fid,Norb,e1int)

                        character(*), intent(in) :: onebody_fid
                        integer, intent(in) :: Norb
                        real, intent(out) :: e1int(:,:)
                        integer, parameter :: read_unit = 10
                        real :: val
                        integer :: ios, Norb2, cnt, p, q, err

                        open(unit=read_unit,file=onebody_fid,form='formatted',status='old',iostat=err)
                        if (err /= 0) then
                           print*,'Error: Onebody file not found!'
                        end if
                 
                        ios = 1
                        do while (ios /= -1) 
                                read(read_unit,fmt=*,iostat=ios) val, cnt 
                        end do
                        Norb2 = int(0.5*(-1.0 + sqrt(1.0+8.0*real(cnt))))

                        if (Norb2 /= Norb) then
                            print*, 'Error : Norb does not match number of integrals!'
                        end if

                        rewind(read_unit)

                        do p = 1,Norb
                          do q = 1,p
                            read(read_unit,fmt=*,iostat=ios) val, cnt
                            e1int(p,q) = val
                            e1int(q,p) = val
                          end do
                        end do

                        close(read_unit)
                     
                end subroutine load_onebody

                subroutine load_twobody(twobody_fid,Norb,e2int,Vnuc)

                        character(*), intent(in) :: twobody_fid
                        integer, intent(in) :: Norb
                        real, intent(out) :: e2int(Norb,Norb,Norb,Norb)
                        real :: val, Vnuc
                        integer :: p, q, r, s, ios, err
                        integer, parameter :: read_unit = 11

                        open(unit=read_unit,file=twobody_fid,form='formatted',status='old',iostat=err)
                        if (err /= 0) then
                            print*, 'Error: Twobody file not found!'
                        end if

                        ios = 1
                        do while (ios /= -1) 
                               read(read_unit,fmt=*,iostat=ios) p, r, q, s, val
                               if (p + q + r + s /= 0) then
                                   e2int(p,q,r,s) = val
                               else
                                   Vnuc = val
                               end if
                        end do

                        close(read_unit)

                 end subroutine load_twobody

                 subroutine get_V(Norb,e2int,vA,vB,vC)

                         integer, intent(in) :: Norb
                         real, intent(in) :: e2int(Norb,Norb,Norb,Norb)
                         real, intent(out) :: vA(Norb,Norb,Norb,Norb), vB(Norb,Norb,Norb,Norb), vC(Norb,Norb,Norb,Norb)
                         integer :: p, q, r, s
                         
                         do p = 1,Norb
                           do q = 1,Norb
                             do r = 1,Norb
                               do s = 1,Norb
                                 vA(p,q,r,s) = e2int(p,q,r,s) - e2int(p,q,s,r)
                                 vB(p,q,r,s) = e2int(p,q,r,s)
                                 vC(p,q,r,s) = e2int(p,q,r,s) - e2int(p,q,s,r)
                               end do
                             end do
                           end do 
                         end do

                  end subroutine get_V

                  

                        


 end module integrals
