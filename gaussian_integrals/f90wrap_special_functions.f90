! Module special_functions defined in file special_functions.f90

subroutine f90wrap_hyp1f1(a, b, ret_hg, x)
    use special_functions, only: hyp1f1
    implicit none
    
    real(4) :: a
    real(4) :: b
    real(4), intent(out) :: ret_hg
    real(4) :: x
    ret_hg = hyp1f1(a=a, b=b, x=x)
end subroutine f90wrap_hyp1f1

! End of module special_functions defined in file special_functions.f90

