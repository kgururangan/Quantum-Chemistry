! Module orbital_types defined in file orbital_types.f90

subroutine f90wrap_basisfcn_t__array__shell(this, nd, dtype, dshape, dloc)
    use orbital_types, only: basisfcn_t
    implicit none
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    integer, intent(in) :: this(2)
    type(basisfcn_t_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%shell)
    dloc = loc(this_ptr%p%shell)
end subroutine f90wrap_basisfcn_t__array__shell

subroutine f90wrap_basisfcn_t__array__coeff(this, nd, dtype, dshape, dloc)
    use orbital_types, only: basisfcn_t
    implicit none
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    integer, intent(in) :: this(2)
    type(basisfcn_t_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%coeff)) then
        dshape(1:1) = shape(this_ptr%p%coeff)
        dloc = loc(this_ptr%p%coeff)
    else
        dloc = 0
    end if
end subroutine f90wrap_basisfcn_t__array__coeff

subroutine f90wrap_basisfcn_t__array__exps(this, nd, dtype, dshape, dloc)
    use orbital_types, only: basisfcn_t
    implicit none
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    integer, intent(in) :: this(2)
    type(basisfcn_t_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%exps)) then
        dshape(1:1) = shape(this_ptr%p%exps)
        dloc = loc(this_ptr%p%exps)
    else
        dloc = 0
    end if
end subroutine f90wrap_basisfcn_t__array__exps

subroutine f90wrap_basisfcn_t__array__origin(this, nd, dtype, dshape, dloc)
    use orbital_types, only: basisfcn_t
    implicit none
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    integer, intent(in) :: this(2)
    type(basisfcn_t_ptr_type) :: this_ptr
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%origin)
    dloc = loc(this_ptr%p%origin)
end subroutine f90wrap_basisfcn_t__array__origin

subroutine f90wrap_basisfcn_t_initialise(this)
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    type(basisfcn_t_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_basisfcn_t_initialise

subroutine f90wrap_basisfcn_t_finalise(this)
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    type(basisfcn_t_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_basisfcn_t_finalise

subroutine f90wrap_basisfcn_t_x1:num_orbs_array__array_getitem__items(f90wrap_this, f90wrap_i, itemsitem)
    
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_x1:num_orbs_array
        type(basisfcn_t), dimension(1:num_orbs) :: items
    end type basisfcn_t_x1:num_orbs_array
    
    type basisfcn_t_x1:num_orbs_array_ptr_type
        type(basisfcn_t_x1:num_orbs_array), pointer :: p => NULL()
    end type basisfcn_t_x1:num_orbs_array_ptr_type
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(basisfcn_t_x1:num_orbs_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: itemsitem(2)
    type(basisfcn_t_ptr_type) :: items_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr%p => this_ptr%p%items(f90wrap_i)
        itemsitem = transfer(items_ptr,itemsitem)
    endif
end subroutine f90wrap_basisfcn_t_x1:num_orbs_array__array_getitem__items

subroutine f90wrap_basisfcn_t_x1:num_orbs_array__array_setitem__items(f90wrap_this, f90wrap_i, itemsitem)
    
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_x1:num_orbs_array
        type(basisfcn_t), dimension(1:num_orbs) :: items
    end type basisfcn_t_x1:num_orbs_array
    
    type basisfcn_t_x1:num_orbs_array_ptr_type
        type(basisfcn_t_x1:num_orbs_array), pointer :: p => NULL()
    end type basisfcn_t_x1:num_orbs_array_ptr_type
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(basisfcn_t_x1:num_orbs_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: itemsitem(2)
    type(basisfcn_t_ptr_type) :: items_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr = transfer(itemsitem,items_ptr)
        this_ptr%p%items(f90wrap_i) = items_ptr%p
    endif
end subroutine f90wrap_basisfcn_t_x1:num_orbs_array__array_setitem__items

subroutine f90wrap_basisfcn_t_x1:num_orbs_array__array_len__items(f90wrap_this, f90wrap_n)
    
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_x1:num_orbs_array
        type(basisfcn_t), dimension(1:num_orbs) :: items
    end type basisfcn_t_x1:num_orbs_array
    
    type basisfcn_t_x1:num_orbs_array_ptr_type
        type(basisfcn_t_x1:num_orbs_array), pointer :: p => NULL()
    end type basisfcn_t_x1:num_orbs_array_ptr_type
    type basisfcn_t_ptr_type
        type(basisfcn_t), pointer :: p => NULL()
    end type basisfcn_t_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(basisfcn_t_x1:num_orbs_array_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    f90wrap_n = size(this_ptr%p%items)
end subroutine f90wrap_basisfcn_t_x1:num_orbs_array__array_len__items

subroutine f90wrap_basisfcn_t_x1:num_orbs_array_initialise(this)
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_x1:num_orbs_array
        type(basisfcn_t), dimension(1:num_orbs) :: items
    end type basisfcn_t_x1:num_orbs_array
    
    type basisfcn_t_x1:num_orbs_array_ptr_type
        type(basisfcn_t_x1:num_orbs_array), pointer :: p => NULL()
    end type basisfcn_t_x1:num_orbs_array_ptr_type
    type(basisfcn_t_x1:num_orbs_array_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_basisfcn_t_x1:num_orbs_array_initialise

subroutine f90wrap_basisfcn_t_x1:num_orbs_array_finalise(this)
    use orbital_types, only: basisfcn_t
    implicit none
    
    type basisfcn_t_x1:num_orbs_array
        type(basisfcn_t), dimension(1:num_orbs) :: items
    end type basisfcn_t_x1:num_orbs_array
    
    type basisfcn_t_x1:num_orbs_array_ptr_type
        type(basisfcn_t_x1:num_orbs_array), pointer :: p => NULL()
    end type basisfcn_t_x1:num_orbs_array_ptr_type
    type(basisfcn_t_x1:num_orbs_array_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_basisfcn_t_x1:num_orbs_array_finalise

! End of module orbital_types defined in file orbital_types.f90

