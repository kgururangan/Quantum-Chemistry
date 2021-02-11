program main

    use einsum_module, only: einsum444
    use blas_module, only: kgemm
    use tensor_type, only: tensor_t
    use permute_module, only: permute2

    implicit none

    integer, parameter :: nu = 20, no = 5
    integer :: a, b, c, d, i, j, k, l, m, n, e, f
    real :: Voovv(no,no,nu,nu), T(nu,nu,no,no), Vvoov(nu,no,no,nu), Vvvvv(nu,nu,nu,nu)
    real, allocatable :: Z(:,:,:,:)
    real, allocatable :: Vp(:,:,:,:), Tp(:,:,:,:), V2(:,:), T2(:,:), Z2(:,:), Zp(:,:,:,:)
    real :: xsum

    call get_matrices(no,nu,Voovv,Vvoov,Vvvvv,T)

    print*,'++++++++++++++++TEST 1: Z(abef) = 0.5*V(mnef)T(abmn)++++++++++++++++'

    allocate(Z(nu,nu,nu,nu))
    xsum = 0.0
    do a = 1,nu
        do b = 1,nu
            do e = 1,nu
                do f = 1,nu
                    Z(a,b,e,f) = 0.0
                    do m = 1,no
                        do n = 1,no
                            Z(a,b,e,f) = Z(a,b,e,f) + &
                            0.5*Voovv(m,n,e,f)*T(a,b,m,n)
                        end do
                    end do 
                    xsum = xsum + Z(a,b,e,f)
                end do 
            end do 
        end do 
    end do 
    print*,'LOOP contraction = ',xsum
    deallocate(Z)

    allocate(Z(nu,nu,nu,nu))
    call einsum444('mnef,abmn->abef',0.5*Voovv,T,Z)
    xsum = 0.0
    do a = 1,nu
        do b = 1,nu
            do e = 1,nu
                do f = 1,nu
                    xsum = xsum + Z(a,b,e,f)
                end do 
            end do 
        end do 
    end do 
    print*,'EINSUM contraction = ',xsum
    deallocate(Z)

    print*,'++++++++++++++++TEST 2: Z(abij) = V(amie)T(bejm)++++++++++++++++'

    allocate(Z(nu,nu,no,no))
    xsum = 0.0
    do a = 1,nu
        do b = 1,nu
            do i = 1,no
                do j = 1,no
                    Z(a,b,i,j) = 0.0
                    do m = 1,no
                        do e = 1,nu
                            Z(a,b,i,j) = Z(a,b,i,j) + &
                            Vvoov(a,m,i,e)*T(b,e,j,m)
                        end do
                    end do 
                    xsum = xsum + Z(a,b,i,j)
                end do 
            end do 
        end do 
    end do 
    print*,'LOOP contraction = ',xsum
    deallocate(Z)

    !!! order in reshape(SOURCE,SHAPE,ORDER) works very strangely...
    ! ORDER is an array defined such that if the RESHAPED array is taken with 
    ! indices given by ORDER, it will return SOURCE

    ! i know what the problem is
    ! idx = [icontr, iuncontr] is identified using the SOURCE ordering \
    ! but reshape requires ORDER to be in terms of the RESHAPED ordering
    ! e.g. consider this
    ! bejm -> mebj
    ! b goes to position 3 so order(1) = 3
    ! e goes to position 2 so order(2) = 2
    ! j goes to position 4 so order(3) = 4
    ! m goes to position 1 so order(4) = 1
    ! whereas looking at the contraction string, we identify 
    ! icontr = [4,2] and iuncotr = [1,3] hence we had 
    ! ORDER = [4,2,1,3]

    ! if we have ORDER_X = [4,2,1,3], we must
    ! this defines mebj
    ! ORDER_Y(k) = q where SOURCE(ORDER_X(q)) = SOURCE(k)
    ! e.g. k = 1 -> SOURCE(ORDER_X(q)) = SOURCE(1) -> ORDER_X(q) = 1
    !               ORDER_X equals 1 at position 3 so q = 3

    allocate(Z(nu,nu,no,no),Vp(nu,no,no,nu),Tp(no,nu,nu,no))
    Vp = reshape(Vvoov,shape=(/nu,no,no,nu/),order=(/1,3,2,4/)) ! amie -> aime
    Tp = reshape(T,shape=(/no,nu,nu,no/),order=(/3,2,4,1/)) ! bejm -> mebj
    ! 4,2,1,3
    ! do i = 1,2 
    !     do j = 1,3   
    !         do k = 1,3
    !             do l = 1,2
    !                 print*,'B(',i,j,k,l,') = ',Tp(i,j,k,l)
    !             end do 
    !         end do 
    !     end do 
    ! end do
    allocate(V2(nu*no,nu*no),T2(nu*no,nu*no),Z2(nu*no,nu*no))
    V2 = reshape(Vp,(/nu*no,nu*no/))
    T2 = reshape(Tp,(/nu*no,nu*no/))
    Z2 = kgemm(V2,T2)
    ! xsum = 0.0
    ! do i = 1,nu*no
    !     do j = 1,nu*no 
    !         xsum = xsum + Z2(i,j) 
    !     end do 
    ! end do 
    ! print*,xsum 
    Z = reshape(Z2,(/nu,no,nu,no/)) ! aibj 
    Z = reshape(Z,(/nu,nu,no,no/),order=(/1,3,2,4/))
    xsum = 0.0
    do a = 1,nu
        do b = 1,nu
            do i = 1,no
                do j = 1,no
                    xsum = xsum + Z(a,b,i,j)
                end do 
            end do 
        end do 
    end do 
    print*,'BLAS contraction = ',xsum
    deallocate(Z)

    allocate(Z(nu,nu,no,no))
    call einsum444('amie,bejm->abij',Vvoov,T,Z)
    xsum = 0.0
    do a = 1,nu
        do b = 1,nu
            do i = 1,no
                do j = 1,no
                    xsum = xsum + Z(a,b,i,j)
                end do 
            end do 
        end do 
    end do 
    print*,'EINSUM contraction = ',xsum
    deallocate(Z)

    print*,'++++++++++++++++TEST 3: Z(abij) = 0.5*V(abef)T(efij)++++++++++++++++'

    allocate(Z(nu,nu,no,no))
    xsum = 0.0
    do a = 1,nu
        do b = 1,nu
            do i = 1,no
                do j = 1,no
                    Z(a,b,i,j) = 0.0
                    do e = 1,nu
                        do f = 1,nu
                            Z(a,b,i,j) = Z(a,b,i,j) + &
                            0.5*Vvvvv(a,b,e,f)*T(e,f,i,j)
                        end do
                    end do 
                    xsum = xsum + Z(a,b,i,j)
                end do 
            end do 
        end do 
    end do 
    print*,'LOOP contraction = ',xsum
    deallocate(Z)

    allocate(Z(nu,nu,no,no))
    call einsum444('abfe,feij->abij',0.5*Vvvvv,T,Z)
    xsum = 0.0
    do a = 1,nu
        do b = 1,nu
            do i = 1,no
                do j = 1,no
                    xsum = xsum + Z(a,b,i,j)
                end do 
            end do 
        end do 
    end do 
    print*,'EINSUM contraction = ',xsum
    deallocate(Z)

    contains 

        subroutine get_matrices(no,nu,Voovv,Vvoov,Vvvvv,T)

            integer, intent(in) :: no, nu
            real, intent(out) :: Voovv(no,no,nu,nu), Vvoov(nu,no,no,nu), Vvvvv(nu,nu,nu,nu), T(nu,nu,no,no)
            integer :: a, b, i, j, m, e, f
            real :: r, xsum, ct

            xsum = 0.0
            ct = 1.0
            do a = 1,nu
                do b = 1,nu
                    do i = 1,no
                        do j = 1,no

                            Voovv(j,i,b,a) = ct 
                            ct = ct + 1.0

                            xsum = xsum + Voovv(j,i,b,a)

                        end do
                    end do
                end do
            end do
            print*,xsum

            xsum = 0.0
            ct = 1.0;
            do a = 1,nu
                do m = 1,no
                    do i = 1,no
                        do e = 1,nu
                            Vvoov(e,i,m,a) = ct 
                            ct = ct + 1.0

                            xsum = xsum + Vvoov(e,i,m,a)

                        end do
                    end do
                end do
            end do
            print*,xsum

            xsum = 0.0
            ct = 1.0
            do a = 1,nu
                do b = 1,nu
                    do e = 1,nu 
                        do f = 1,nu 
                            Vvvvv(f,e,b,a) = ct 
                            ct = ct + 1.0

                            xsum = xsum + Vvvvv(f,e,b,a)
  
                        end do 
                    end do 
                end do 
            end do
            print*,xsum

            xsum = 0.0
            ct = 1.0
            do i = 1,no
                do j = 1,no
                    do a = 1,nu
                        do b = 1,nu

                            T(b,a,j,i) = ct
                            ct = ct + 1.0

                            xsum = xsum + T(b,a,j,i)

                        end do
                    end do
                end do
            end do
            print*,xsum

        end subroutine get_matrices

end program main