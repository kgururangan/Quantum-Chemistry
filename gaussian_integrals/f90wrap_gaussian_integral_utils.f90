! Module gaussian_integral_utils defined in file gaussian_integral_utils.f90

subroutine f90wrap_fact(ret_val, n)
    use gaussian_integral_utils, only: fact
    implicit none
    
    real, intent(out) :: ret_val
    integer, intent(in) :: n
    ret_val = fact(n=n)
end subroutine f90wrap_fact

subroutine f90wrap_nchoosek(n, ret_val, k)
    use gaussian_integral_utils, only: nchoosek
    implicit none
    
    integer, intent(in) :: n
    real, intent(out) :: ret_val
    integer, intent(in) :: k
    ret_val = nchoosek(n=n, k=k)
end subroutine f90wrap_nchoosek

subroutine f90wrap_fact2(ret_val, n)
    use gaussian_integral_utils, only: fact2
    implicit none
    
    real, intent(out) :: ret_val
    integer, intent(in) :: n
    ret_val = fact2(n=n)
end subroutine f90wrap_fact2

subroutine f90wrap_boys(ret_val, t)
    use gaussian_integral_utils, only: boys
    implicit none
    
    real, intent(out) :: ret_val
    real, intent(in) :: t
    ret_val = boys(t=t)
end subroutine f90wrap_boys

subroutine f90wrap_boys_v2(n, ret_val, t)
    use gaussian_integral_utils, only: boys_v2
    implicit none
    
    integer, intent(in) :: n
    real, intent(out) :: ret_val
    real, intent(in) :: t
    ret_val = boys_v2(n=n, T=t)
end subroutine f90wrap_boys_v2

subroutine f90wrap_gaussian_product(c_exp, rc, c, alpha, ra, beta, rb, n0, n1, n2)
    use gaussian_integral_utils, only: gaussian_product
    implicit none
    
    real, intent(out) :: c_exp
    real, intent(inout), dimension(n0) :: rc
    real, intent(out) :: c
    real, intent(in) :: alpha
    real, intent(in), dimension(n1) :: ra
    real, intent(in) :: beta
    real, intent(in), dimension(n2) :: rb
    integer :: n0
    !f2py intent(hide), depend(rc) :: n0 = shape(rc,0)
    integer :: n1
    !f2py intent(hide), depend(ra) :: n1 = shape(ra,0)
    integer :: n2
    !f2py intent(hide), depend(rb) :: n2 = shape(rb,0)
    call gaussian_product(c_exp=c_exp, rC=rc, C=c, alpha=alpha, rA=ra, beta=beta, rB=rb)
end subroutine f90wrap_gaussian_product

subroutine f90wrap_gaussian_normfactor(shell, ret_val, alpha, n0)
    use gaussian_integral_utils, only: gaussian_normfactor
    implicit none
    
    integer, intent(in), dimension(n0) :: shell
    real, intent(out) :: ret_val
    real, intent(in) :: alpha
    integer :: n0
    !f2py intent(hide), depend(shell) :: n0 = shape(shell,0)
    ret_val = gaussian_normfactor(shell=shell, alpha=alpha)
end subroutine f90wrap_gaussian_normfactor

subroutine f90wrap_gaussian_coeff(l, m, xa, xb, ret_val, k)
    use gaussian_integral_utils, only: gaussian_coeff
    implicit none
    
    integer, intent(in) :: l
    integer, intent(in) :: m
    real, intent(in) :: xa
    real, intent(in) :: xb
    real, intent(out) :: ret_val
    integer, intent(in) :: k
    ret_val = gaussian_coeff(l=l, m=m, xa=xa, xb=xb, k=k)
end subroutine f90wrap_gaussian_coeff

! End of module gaussian_integral_utils defined in file gaussian_integral_utils.f90

