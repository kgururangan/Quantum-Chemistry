program main

    use ao_integrals, only : sAB, tAB, vAB, vABCD
    use gaussian_integration, only: overlap_integral

    real, dimension(2,2) :: Smat, Tmat, Vmat
    real, dimension(2,2,2,2) :: VVmat
    real, dimension(2,3) :: atom_coords
    real, dimension(1:2) :: Z
    real, dimension(1:3) :: xi1, xi2, cm1, cm2
    real :: zeta_H

    Z = [1.0, 1.0]

    zeta_H = 1.24

    atom_coords(1,:) = [-0.5, 0.0, 0.0]
    atom_coords(2,:) = [0.5, 0.0, 0.0]

    xi1 = [0.109818*zeta_H**2.0, 0.405771*zeta_H**2.0, 2.22766*zeta_H**2.0]
    cm1 = [0.444635, 0.535328, 0.154329]
    xi2 = [0.109818*zeta_H**2.0, 0.405771*zeta_H**2.0, 2.22766*zeta_H**2.0]
    cm2 = [0.444635, 0.535328, 0.154329]

    Smat(1,1) = sAB(xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:))
    Smat(1,2) = sAB(xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:))
    Smat(2,1) = sAB(xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:))
    Smat(2,2) = sAB(xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:))
    
    Tmat(1,1) = tAB(xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:))
    Tmat(1,2) = tAB(xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:))
    Tmat(2,1) = tAB(xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:))
    Tmat(2,2) = tAB(xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:))

    Vmat(1,1) = vAB(xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),atom_coords,Z)
    Vmat(1,2) = vAB(xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),atom_coords,Z)
    Vmat(2,1) = vAB(xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),atom_coords,Z)
    Vmat(2,2) = vAB(xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),atom_coords,Z)

    VVmat(1,1,1,1) = vABCD(xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:))
    VVmat(1,1,1,2) = vABCD(xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:))
    VVmat(1,1,2,1) = vABCD(xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:))
    VVmat(1,2,1,1) = vABCD(xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:))
    VVmat(2,1,1,1) = vABCD(xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:))
    VVmat(2,2,1,1) = vABCD(xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:))
    VVmat(2,1,2,1) = vABCD(xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:))
    VVmat(2,1,1,2) = vABCD(xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:))
    VVmat(1,2,2,1) = vABCD(xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:))
    VVmat(1,2,1,2) = vABCD(xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:))
    VVmat(1,1,2,2) = vABCD(xi1,cm1,atom_coords(1,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:))
    VVmat(2,2,2,1) = vABCD(xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:))
    VVmat(2,2,1,2) = vABCD(xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:))
    VVmat(2,1,2,2) = vABCD(xi2,cm2,atom_coords(2,:),xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:))
    VVmat(1,2,2,2) = vABCD(xi1,cm1,atom_coords(1,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:))
    VVmat(2,2,2,2) = vABCD(xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:),xi2,cm2,atom_coords(2,:))


    print*, '==============OVERLAP MATRIX==============='
    do i = 1,2
        do j = 1,2
            print*, 'S(',i,',',j,') = ',Smat(i,j)
        enddo
    enddo
    
    print*, '==============KINETIC ENERGY MATRIX==============='
    do i = 1,2
        do j = 1,2
            print*, 'T(',i,',',j,') = ',Tmat(i,j)
        enddo
    enddo

    print*, '==============POTENTIAL ENERGY MATRIX==============='
    do i = 1,2
        do j = 1,2
            print*, 'V(',i,',',j,') = ',Vmat(i,j)
        enddo
    enddo


    print*, '==============ERI MATRIX==============='
    do i = 1,2
        do j = 1,2
           do k = 1,2
              do l = 1,2
                print*, 'VV(',i,',',j,',',k,',',l,') = ',VVmat(i,j,k,l)
              enddo
           enddo
        enddo
    enddo

end program main
