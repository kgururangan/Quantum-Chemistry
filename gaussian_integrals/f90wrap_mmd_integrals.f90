! Module mmd_integrals defined in file mmd_integrals.f90

subroutine f90wrap_sab(orba, ret_val, orbb)
    use mmd_integrals, only: sab
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    type(basisfcn_t_ptr_type) :: orba_ptr
    integer, intent(in), dimension(2) :: orba
    real, intent(out) :: ret_val
    type(basisfcn_t_ptr_type) :: orbb_ptr
    integer, intent(in), dimension(2) :: orbb
    orba_ptr = transfer(orba, orba_ptr)
    orbb_ptr = transfer(orbb, orbb_ptr)
    ret_val = sab(orbA=orba_ptr%p, orbB=orbb_ptr%p)
end subroutine f90wrap_sab

subroutine f90wrap_tab(orba, ret_val, orbb)
    use orbital_types, only: basisfcn_t
    use mmd_integrals, only: tab
    implicit none
    
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    type(basisfcn_t_ptr_type) :: orba_ptr
    integer, intent(in), dimension(2) :: orba
    real, intent(out) :: ret_val
    type(basisfcn_t_ptr_type) :: orbb_ptr
    integer, intent(in), dimension(2) :: orbb
    orba_ptr = transfer(orba, orba_ptr)
    orbb_ptr = transfer(orbb, orbb_ptr)
    ret_val = tab(orbA=orba_ptr%p, orbB=orbb_ptr%p)
end subroutine f90wrap_tab

subroutine f90wrap_vab(orba, orbb, rc, ret_val, z, n0)
    use mmd_integrals, only: vab
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    type(basisfcn_t_ptr_type) :: orba_ptr
    integer, intent(in), dimension(2) :: orba
    type(basisfcn_t_ptr_type) :: orbb_ptr
    integer, intent(in), dimension(2) :: orbb
    real, intent(in), dimension(n0) :: rc
    real, intent(out) :: ret_val
    real, intent(in) :: z
    integer :: n0
    !f2py intent(hide), depend(rc) :: n0 = shape(rc,0)
    orba_ptr = transfer(orba, orba_ptr)
    orbb_ptr = transfer(orbb, orbb_ptr)
    ret_val = vab(orbA=orba_ptr%p, orbB=orbb_ptr%p, rC=rc, Z=z)
end subroutine f90wrap_vab

subroutine f90wrap_eri(orba, orbb, orbc, ret_val, orbd)
    use mmd_integrals, only: eri
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    type(basisfcn_t_ptr_type) :: orba_ptr
    integer, intent(in), dimension(2) :: orba
    type(basisfcn_t_ptr_type) :: orbb_ptr
    integer, intent(in), dimension(2) :: orbb
    type(basisfcn_t_ptr_type) :: orbc_ptr
    integer, intent(in), dimension(2) :: orbc
    real, intent(out) :: ret_val
    type(basisfcn_t_ptr_type) :: orbd_ptr
    integer, intent(in), dimension(2) :: orbd
    orba_ptr = transfer(orba, orba_ptr)
    orbb_ptr = transfer(orbb, orbb_ptr)
    orbc_ptr = transfer(orbc, orbc_ptr)
    orbd_ptr = transfer(orbd, orbd_ptr)
    ret_val = eri(orbA=orba_ptr%p, orbB=orbb_ptr%p, orbC=orbc_ptr%p, orbD=orbd_ptr%p)
end subroutine f90wrap_eri

! End of module mmd_integrals defined in file mmd_integrals.f90

