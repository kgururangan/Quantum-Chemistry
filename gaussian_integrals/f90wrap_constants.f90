! Module constants defined in file constants.f90

subroutine f90wrap_constants__get__pi(f90wrap_pi)
    use constants, only: constants_pi => pi
    implicit none
    real, intent(out) :: f90wrap_pi
    
    f90wrap_pi = constants_pi
end subroutine f90wrap_constants__get__pi

subroutine f90wrap_constants__get__hbar(f90wrap_hbar)
    use constants, only: constants_hbar => hbar
    implicit none
    real(4), intent(out) :: f90wrap_hbar
    
    f90wrap_hbar = constants_hbar
end subroutine f90wrap_constants__get__hbar

subroutine f90wrap_constants__get__e(f90wrap_e)
    use constants, only: constants_e => e
    implicit none
    real, intent(out) :: f90wrap_e
    
    f90wrap_e = constants_e
end subroutine f90wrap_constants__get__e

! End of module constants defined in file constants.f90

