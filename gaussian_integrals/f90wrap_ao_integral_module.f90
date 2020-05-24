! Module ao_integral_module defined in file ao_integral_module.f90

subroutine f90wrap_get_ao_integrals(smat, zmat, vvmat, num_orbs, num_atoms, orbs, atom_coords, z, n0, n1, n2, n3, n4, &
    n5, n6, n7, n8, n9, n10)
    use orbital_types, only: basisfcn_t
    use ao_integral_module, only: get_ao_integrals
    implicit none
    
    type basisfcn_t_x1:num_orbs_array
        type(basisfcn_t), dimension(1:num_orbs) :: items
    end type basisfcn_t_x1:num_orbs_array
    
    type basisfcn_t_x1:num_orbs_array_ptr_type
        type(basisfcn_t_x1:num_orbs_array), pointer :: p => NULL()
    end type basisfcn_t_x1:num_orbs_array_ptr_type
    real, intent(inout), dimension(n0,n1) :: smat
    real, intent(inout), dimension(n2,n3) :: zmat
    real, intent(inout), dimension(n4,n5,n6,n7) :: vvmat
    integer, intent(in) :: num_orbs
    integer, intent(in) :: num_atoms
    type(basisfcn_t_x1:num_orbs_array_ptr_type) :: orbs_ptr
    integer, intent(in), dimension(2) :: orbs
    real, intent(in), dimension(n8,n9) :: atom_coords
    real, intent(in), dimension(n10) :: z
    integer :: n0
    !f2py intent(hide), depend(smat) :: n0 = shape(smat,0)
    integer :: n1
    !f2py intent(hide), depend(smat) :: n1 = shape(smat,1)
    integer :: n2
    !f2py intent(hide), depend(zmat) :: n2 = shape(zmat,0)
    integer :: n3
    !f2py intent(hide), depend(zmat) :: n3 = shape(zmat,1)
    integer :: n4
    !f2py intent(hide), depend(vvmat) :: n4 = shape(vvmat,0)
    integer :: n5
    !f2py intent(hide), depend(vvmat) :: n5 = shape(vvmat,1)
    integer :: n6
    !f2py intent(hide), depend(vvmat) :: n6 = shape(vvmat,2)
    integer :: n7
    !f2py intent(hide), depend(vvmat) :: n7 = shape(vvmat,3)
    integer :: n8
    !f2py intent(hide), depend(atom_coords) :: n8 = shape(atom_coords,0)
    integer :: n9
    !f2py intent(hide), depend(atom_coords) :: n9 = shape(atom_coords,1)
    integer :: n10
    !f2py intent(hide), depend(z) :: n10 = shape(z,0)
    orbs_ptr = transfer(orbs, orbs_ptr)
    call get_ao_integrals(Smat=smat, Zmat=zmat, VVmat=vvmat, num_orbs=num_orbs, num_atoms=num_atoms, orbs=orbs_ptr%p%items, &
        atom_coords=atom_coords, Z=z)
end subroutine f90wrap_get_ao_integrals

! End of module ao_integral_module defined in file ao_integral_module.f90

