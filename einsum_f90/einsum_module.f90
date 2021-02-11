module einsum_module

    !use tensor_type, only: tensor_t
    use blas_module, only: kgemm 
    use sort_module, only: argsort_int, argsort
    use permute_module, only: permute2

    implicit none

    contains

        subroutine einsum444(str,A,B,C)

            character, intent(in) :: str(1:15)
            real, intent(in) :: A(:,:,:,:), B(:,:,:,:)
            real, intent(out) :: C(:,:,:,:)
            real, allocatable :: A2(:,:), B2(:,:), C2(:,:),&
                                 Ap(:,:,:,:), Bp(:,:,:,:), Cp(:,:,:,:)
            character :: s1(1:4), s2(1:4), s3(1:4)
            integer :: i, j, k, l, idx, ct1, ct2, ct3, &
                       idxA(1:4), idxB(1:4), idxC(1:4), idxC2(1:4)
            integer :: shapeA(1:4), shapeB(1:4), shapeC(1:4), n1, n2, n3
            integer :: temp1(1:4), temp2(1:4)
            real :: xsum

            s1 = str(1:4)
            s2 = str(6:9)
            s3 = str(12:)

            shapeA = shape(A)
            shapeB = shape(B)
            shapeC = shape(C)

            ct1 = 1
            ct2 = 1
            ct3 = 1

            n1 = 1
            n2 = 1
            n3 = 1

            do i = 1,4
                idx = findloc(s2,s1(i),1)
                if (idx == 0) then ! i is an output index
                    idx = findloc(s3,s1(i),1)
                    idxC(ct3) = idx
                    ct3 = ct3 + 1
                    idxA(ct1) = i 
                    n1 = n1 * shapeA(i)
                    ct1 = ct1 + 1 
                else ! i is contracted
                    idxB(ct2) = idx 
                    ct2 = ct2 + 1
                end if
            end do

            do i = 1,4
                idx = findloc(s3,s2(i),1)
                if (idx /= 0) then ! idx is an output index
                    idxC(ct3) = idx
                    ct3 = ct3 + 1
                end if
            end do

            temp1 = (/1,2,3,4/)
            ct3 = 1
            do i = 1,4
                if (any(temp1(i) == idxA(1:ct1-1))) then
                    cycle
                else 
                    idxA(ct1+ct3-1) = i 
                    ct3 = ct3 + 1
                    n2 = n2 * shapeA(i)
                end if
            end do

            temp2 = (/1,2,3,4/)
            ct3 = 1
            do i = 1,4 
                if (any(temp2(i) == idxB(1:ct2-1))) then 
                    cycle 
                else 
                    idxB(ct2+ct3-1) = i 
                    ct3 = ct3 + 1 
                    n3 = n3 * shapeB(i)
                end if 
            end do

            shapeA = shapeA(idxA)
            shapeB = shapeB(idxB)
            shapeC = shapeC(idxC)

            ! print*,'n1 = ',n1,'n2 = ',n2,'n3 = ',n3
            ! print*,'idxA = ',idxA
            ! print*,'shapeA = ',shapeA
            ! print*,'idxB = ',idxB
            ! print*,'shapeB = ',shapeB
            ! print*,'idxC = ',idxC
            ! print*,'shapeC = ',shapeC

            allocate(Ap(shapeA(1),shapeA(2),shapeA(3),shapeA(4)))
            allocate(Bp(shapeB(1),shapeB(2),shapeB(3),shapeB(4)))
            allocate(Cp(shapeC(1),shapeC(2),shapeC(3),shapeC(4)))
            allocate(A2(n1,n2),B2(n2,n3),C2(n1,n3))

            !Ap = reshape(A,shape=shapeA,order=idxA)
            call permute2(A,Ap,idxA)

            ! There's a problem here with how B is being permuted or something...
            !Bp = reshape(B,shape=shapeB,order=idxB) 

            call permute2(B,Bp,idxB)

            ! temp1 = shape(Bp)
            ! do i = 1,temp1(1)
            !     do j = 1,temp1(2)
            !         do k = 1,temp1(3)
            !             do l = 1,temp1(4)
            !                 print*,'B(',i,j,k,l,') = ',Bp(i,j,k,l)
            !             end do 
            !         end do
            !     end do 
            ! end do  

            A2 = reshape(Ap,shape=(/n1,n2/))
            B2 = reshape(Bp,shape=(/n2,n3/))
            C2 = kgemm(A2,B2)
            Cp = reshape(C2,shape=shapeC)
            idxC2 = argsort_int(idxC)
            C = reshape(Cp,shape=shapeC(idxC2),order=idxC2)


        deallocate(Ap,Bp,A2,B2,C2)

        end subroutine einsum444



end module einsum_module
