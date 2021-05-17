module dgeev_module

        implicit none

        contains

                subroutine eig(A,VR,wR,wI)

                        real, intent(in) :: A(:,:)
                        real, allocatable, intent(out) :: VR(:,:), wR(:), wI(:)
                        real, allocatable :: work(:), VL(:,:)
                        integer :: N, lwork, info

                        N = ubound(A,1)
                        lwork = 4 * N
                        allocate(VR(N,N),VL(N,N),wR(N),wI(N),work(lwork))

                        call dgeev('N','V',N,A,N,wR,wI,VL,N,VR,N,work,lwork,info)
                        
                        deallocate(work,VL)

                        !if (info /= 0) then
                        !        write(*,'(a/)') 'Problem in DGEEV:'
                        !        write(*,fmt=*) 'info = ',info
                        !end if

                 end subroutine eig

                subroutine eigh(A,VR,wR)

                        real, intent(in) :: A(:,:)
                        real, allocatable, intent(out) :: VR(:,:), wR(:)
                        real, allocatable :: work(:)
                        integer :: N, lwork, info

                        N = ubound(A,1)
                        lwork = 4 * N
                        allocate(VR(N,N),wR(N),work(lwork))

                        call dsyev('V','L',N,A,N,wR,work,lwork,info)

                        VR = A
                        
                        deallocate(work)

                        !if (info /= 0) then
                        !        write(*,'(a/)') 'Problem in DGEEV:'
                        !        write(*,fmt=*) 'info = ',info
                        !end if

                 end subroutine eigh

 end module dgeev_module
