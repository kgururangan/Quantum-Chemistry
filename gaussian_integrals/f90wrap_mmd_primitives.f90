! Module mmd_primitives defined in file mmd_primitives.f90

subroutine f90wrap_efcn(la, lb, t, dist_ab, alpha, ret_val, beta)
    use mmd_primitives, only: efcn
    implicit none
    
    integer, intent(in) :: la
    integer, intent(in) :: lb
    integer, intent(in) :: t
    real, intent(in) :: dist_ab
    real, intent(in) :: alpha
    real, intent(out) :: ret_val
    real, intent(in) :: beta
    ret_val = efcn(la=la, lb=lb, t=t, dist_AB=dist_ab, alpha=alpha, beta=beta)
end subroutine f90wrap_efcn

subroutine f90wrap_rfcn(t, u, v, n, p, pcx, pcy, pcz, ret_val, rpc)
    use mmd_primitives, only: rfcn
    implicit none
    
    integer, intent(in) :: t
    integer, intent(in) :: u
    integer, intent(in) :: v
    integer, intent(in) :: n
    real, intent(in) :: p
    real, intent(in) :: pcx
    real, intent(in) :: pcy
    real, intent(in) :: pcz
    real, intent(out) :: ret_val
    real, intent(in) :: rpc
    ret_val = rfcn(t=t, u=u, v=v, n=n, p=p, PCx=pcx, PCy=pcy, PCz=pcz, RPC=rpc)
end subroutine f90wrap_rfcn

subroutine f90wrap_overlap(a, lmn1, ra, b, lmn2, ret_val, rb, n0, n1, n2, n3)
    use mmd_primitives, only: overlap
    implicit none
    
    real, intent(in) :: a
    integer, intent(in), dimension(n0) :: lmn1
    real, intent(in), dimension(n1) :: ra
    real, intent(in) :: b
    integer, intent(in), dimension(n2) :: lmn2
    real, intent(out) :: ret_val
    real, intent(in), dimension(n3) :: rb
    integer :: n0
    !f2py intent(hide), depend(lmn1) :: n0 = shape(lmn1,0)
    integer :: n1
    !f2py intent(hide), depend(ra) :: n1 = shape(ra,0)
    integer :: n2
    !f2py intent(hide), depend(lmn2) :: n2 = shape(lmn2,0)
    integer :: n3
    !f2py intent(hide), depend(rb) :: n3 = shape(rb,0)
    ret_val = overlap(a=a, lmn1=lmn1, rA=ra, b=b, lmn2=lmn2, rB=rb)
end subroutine f90wrap_overlap

subroutine f90wrap_kinetic(a, lmn1, ra, b, lmn2, ret_val, rb, n0, n1, n2, n3)
    use mmd_primitives, only: kinetic
    implicit none
    
    real, intent(in) :: a
    integer, intent(in), dimension(n0) :: lmn1
    real, intent(in), dimension(n1) :: ra
    real, intent(in) :: b
    integer, intent(in), dimension(n2) :: lmn2
    real, intent(out) :: ret_val
    real, intent(in), dimension(n3) :: rb
    integer :: n0
    !f2py intent(hide), depend(lmn1) :: n0 = shape(lmn1,0)
    integer :: n1
    !f2py intent(hide), depend(ra) :: n1 = shape(ra,0)
    integer :: n2
    !f2py intent(hide), depend(lmn2) :: n2 = shape(lmn2,0)
    integer :: n3
    !f2py intent(hide), depend(rb) :: n3 = shape(rb,0)
    ret_val = kinetic(a=a, lmn1=lmn1, rA=ra, b=b, lmn2=lmn2, rB=rb)
end subroutine f90wrap_kinetic

subroutine f90wrap_nuclear_attraction(a, lmn1, ra, b, lmn2, rb, ret_val, rc, n0, n1, n2, n3, n4)
    use mmd_primitives, only: nuclear_attraction
    implicit none
    
    real, intent(in) :: a
    integer, intent(in), dimension(n0) :: lmn1
    real, intent(in), dimension(n1) :: ra
    real, intent(in) :: b
    integer, intent(in), dimension(n2) :: lmn2
    real, intent(in), dimension(n3) :: rb
    real, intent(out) :: ret_val
    real, intent(in), dimension(n4) :: rc
    integer :: n0
    !f2py intent(hide), depend(lmn1) :: n0 = shape(lmn1,0)
    integer :: n1
    !f2py intent(hide), depend(ra) :: n1 = shape(ra,0)
    integer :: n2
    !f2py intent(hide), depend(lmn2) :: n2 = shape(lmn2,0)
    integer :: n3
    !f2py intent(hide), depend(rb) :: n3 = shape(rb,0)
    integer :: n4
    !f2py intent(hide), depend(rc) :: n4 = shape(rc,0)
    ret_val = nuclear_attraction(a=a, lmn1=lmn1, rA=ra, b=b, lmn2=lmn2, rB=rb, rC=rc)
end subroutine f90wrap_nuclear_attraction

subroutine f90wrap_electron_repulsion(a, lmn1, ra, b, lmn2, rb, c, lmn3, rc, d, lmn4, ret_val, rd, n0, n1, n2, n3, n4, &
    n5, n6, n7)
    use mmd_primitives, only: electron_repulsion
    implicit none
    
    real, intent(in) :: a
    integer, intent(in), dimension(n0) :: lmn1
    real, intent(in), dimension(n1) :: ra
    real, intent(in) :: b
    integer, intent(in), dimension(n2) :: lmn2
    real, intent(in), dimension(n3) :: rb
    real, intent(in) :: c
    integer, intent(in), dimension(n4) :: lmn3
    real, intent(in), dimension(n5) :: rc
    real, intent(in) :: d
    integer, intent(in), dimension(n6) :: lmn4
    real, intent(out) :: ret_val
    real, intent(in), dimension(n7) :: rd
    integer :: n0
    !f2py intent(hide), depend(lmn1) :: n0 = shape(lmn1,0)
    integer :: n1
    !f2py intent(hide), depend(ra) :: n1 = shape(ra,0)
    integer :: n2
    !f2py intent(hide), depend(lmn2) :: n2 = shape(lmn2,0)
    integer :: n3
    !f2py intent(hide), depend(rb) :: n3 = shape(rb,0)
    integer :: n4
    !f2py intent(hide), depend(lmn3) :: n4 = shape(lmn3,0)
    integer :: n5
    !f2py intent(hide), depend(rc) :: n5 = shape(rc,0)
    integer :: n6
    !f2py intent(hide), depend(lmn4) :: n6 = shape(lmn4,0)
    integer :: n7
    !f2py intent(hide), depend(rd) :: n7 = shape(rd,0)
    ret_val = electron_repulsion(a=a, lmn1=lmn1, rA=ra, b=b, lmn2=lmn2, rB=rb, c=c, lmn3=lmn3, rC=rc, d=d, lmn4=lmn4, rD=rd)
end subroutine f90wrap_electron_repulsion

! End of module mmd_primitives defined in file mmd_primitives.f90

