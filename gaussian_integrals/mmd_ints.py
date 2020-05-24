from __future__ import print_function, absolute_import, division
import _mmd_ints
import f90wrap.runtime
import logging

class Ao_Integral_Module(f90wrap.runtime.FortranModule):
    """
    Module ao_integral_module
    
    
    Defined at ao_integral_module.f90 lines 1-45
    
    """
    @staticmethod
    def get_ao_integrals(smat, zmat, vvmat, num_orbs, num_atoms, orbs, atom_coords, \
        z):
        """
        get_ao_integrals(smat, zmat, vvmat, num_orbs, num_atoms, orbs, atom_coords, z)
        
        
        Defined at ao_integral_module.f90 lines 3-45
        
        Parameters
        ----------
        smat : float array
        zmat : float array
        vvmat : float array
        num_orbs : int
        num_atoms : int
        orbs : Basisfcn_T_X1:Num_Orbs_Array
        	super-type
        
        atom_coords : float array
        z : float array
        
        """
        _mmd_ints.f90wrap_get_ao_integrals(smat=smat, zmat=zmat, vvmat=vvmat, \
            num_orbs=num_orbs, num_atoms=num_atoms, orbs=orbs._handle, \
            atom_coords=atom_coords, z=z)
    
    _dt_array_initialisers = []
    

ao_integral_module = Ao_Integral_Module()

class Mmd_Integrals(f90wrap.runtime.FortranModule):
    """
    Module mmd_integrals
    
    
    Defined at mmd_integrals.f90 lines 1-87
    
    """
    @staticmethod
    def sab(self, orbb):
        """
        val = sab(self, orbb)
        
        
        Defined at mmd_integrals.f90 lines 7-22
        
        Parameters
        ----------
        orba : Basisfcn_T
        orbb : Basisfcn_T
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_sab(orba=self._handle, orbb=orbb._handle)
        return val
    
    @staticmethod
    def tab(self, orbb):
        """
        val = tab(self, orbb)
        
        
        Defined at mmd_integrals.f90 lines 24-39
        
        Parameters
        ----------
        orba : Basisfcn_T
        orbb : Basisfcn_T
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_tab(orba=self._handle, orbb=orbb._handle)
        return val
    
    @staticmethod
    def vab(self, orbb, rc, z):
        """
        val = vab(self, orbb, rc, z)
        
        
        Defined at mmd_integrals.f90 lines 41-59
        
        Parameters
        ----------
        orba : Basisfcn_T
        orbb : Basisfcn_T
        rc : float array
        z : float
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_vab(orba=self._handle, orbb=orbb._handle, rc=rc, z=z)
        return val
    
    @staticmethod
    def eri(self, orbb, orbc, orbd):
        """
        val = eri(self, orbb, orbc, orbd)
        
        
        Defined at mmd_integrals.f90 lines 61-86
        
        Parameters
        ----------
        orba : Basisfcn_T
        orbb : Basisfcn_T
        orbc : Basisfcn_T
        orbd : Basisfcn_T
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_eri(orba=self._handle, orbb=orbb._handle, \
            orbc=orbc._handle, orbd=orbd._handle)
        return val
    
    _dt_array_initialisers = []
    

mmd_integrals = Mmd_Integrals()

class Mmd_Primitives(f90wrap.runtime.FortranModule):
    """
    Module mmd_primitives
    
    
    Defined at mmd_primitives.f90 lines 1-148
    
    """
    @staticmethod
    def efcn(la, lb, t, dist_ab, alpha, beta):
        """
        val = efcn(la, lb, t, dist_ab, alpha, beta)
        
        
        Defined at mmd_primitives.f90 lines 8-26
        
        Parameters
        ----------
        la : int
        lb : int
        t : int
        dist_ab : float
        alpha : float
        beta : float
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_efcn(la=la, lb=lb, t=t, dist_ab=dist_ab, alpha=alpha, \
            beta=beta)
        return val
    
    @staticmethod
    def rfcn(t, u, v, n, p, pcx, pcy, pcz, rpc):
        """
        val = rfcn(t, u, v, n, p, pcx, pcy, pcz, rpc)
        
        
        Defined at mmd_primitives.f90 lines 28-51
        
        Parameters
        ----------
        t : int
        u : int
        v : int
        n : int
        p : float
        pcx : float
        pcy : float
        pcz : float
        rpc : float
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_rfcn(t=t, u=u, v=v, n=n, p=p, pcx=pcx, pcy=pcy, pcz=pcz, \
            rpc=rpc)
        return val
    
    @staticmethod
    def overlap(a, lmn1, ra, b, lmn2, rb):
        """
        val = overlap(a, lmn1, ra, b, lmn2, rb)
        
        
        Defined at mmd_primitives.f90 lines 53-61
        
        Parameters
        ----------
        a : float
        lmn1 : int array
        ra : float array
        b : float
        lmn2 : int array
        rb : float array
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_overlap(a=a, lmn1=lmn1, ra=ra, b=b, lmn2=lmn2, rb=rb)
        return val
    
    @staticmethod
    def kinetic(a, lmn1, ra, b, lmn2, rb):
        """
        val = kinetic(a, lmn1, ra, b, lmn2, rb)
        
        
        Defined at mmd_primitives.f90 lines 63-78
        
        Parameters
        ----------
        a : float
        lmn1 : int array
        ra : float array
        b : float
        lmn2 : int array
        rb : float array
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_kinetic(a=a, lmn1=lmn1, ra=ra, b=b, lmn2=lmn2, rb=rb)
        return val
    
    @staticmethod
    def nuclear_attraction(a, lmn1, ra, b, lmn2, rb, rc):
        """
        val = nuclear_attraction(a, lmn1, ra, b, lmn2, rb, rc)
        
        
        Defined at mmd_primitives.f90 lines 80-100
        
        Parameters
        ----------
        a : float
        lmn1 : int array
        ra : float array
        b : float
        lmn2 : int array
        rb : float array
        rc : float array
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_nuclear_attraction(a=a, lmn1=lmn1, ra=ra, b=b, \
            lmn2=lmn2, rb=rb, rc=rc)
        return val
    
    @staticmethod
    def electron_repulsion(a, lmn1, ra, b, lmn2, rb, c, lmn3, rc, d, lmn4, rd):
        """
        val = electron_repulsion(a, lmn1, ra, b, lmn2, rb, c, lmn3, rc, d, lmn4, rd)
        
        
        Defined at mmd_primitives.f90 lines 102-147
        
        Parameters
        ----------
        a : float
        lmn1 : int array
        ra : float array
        b : float
        lmn2 : int array
        rb : float array
        c : float
        lmn3 : int array
        rc : float array
        d : float
        lmn4 : int array
        rd : float array
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_electron_repulsion(a=a, lmn1=lmn1, ra=ra, b=b, \
            lmn2=lmn2, rb=rb, c=c, lmn3=lmn3, rc=rc, d=d, lmn4=lmn4, rd=rd)
        return val
    
    _dt_array_initialisers = []
    

mmd_primitives = Mmd_Primitives()

class Gaussian_Integral_Utils(f90wrap.runtime.FortranModule):
    """
    Module gaussian_integral_utils
    
    
    Defined at gaussian_integral_utils.f90 lines 1-88
    
    """
    @staticmethod
    def fact(n):
        """
        val = fact(n)
        
        
        Defined at gaussian_integral_utils.f90 lines 6-9
        
        Parameters
        ----------
        n : int
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_fact(n=n)
        return val
    
    @staticmethod
    def nchoosek(n, k):
        """
        val = nchoosek(n, k)
        
        
        Defined at gaussian_integral_utils.f90 lines 11-14
        
        Parameters
        ----------
        n : int
        k : int
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_nchoosek(n=n, k=k)
        return val
    
    @staticmethod
    def fact2(n):
        """
        val = fact2(n)
        
        
        Defined at gaussian_integral_utils.f90 lines 16-31
        
        Parameters
        ----------
        n : int
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_fact2(n=n)
        return val
    
    @staticmethod
    def boys(t):
        """
        val = boys(t)
        
        
        Defined at gaussian_integral_utils.f90 lines 33-43
        
        Parameters
        ----------
        t : float
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_boys(t=t)
        return val
    
    @staticmethod
    def boys_v2(n, t):
        """
        val = boys_v2(n, t)
        
        
        Defined at gaussian_integral_utils.f90 lines 45-51
        
        Parameters
        ----------
        n : int
        t : float
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_boys_v2(n=n, t=t)
        return val
    
    @staticmethod
    def gaussian_product(rc, alpha, ra, beta, rb):
        """
        c_exp, c = gaussian_product(rc, alpha, ra, beta, rb)
        
        
        Defined at gaussian_integral_utils.f90 lines 53-60
        
        Parameters
        ----------
        rc : float array
        alpha : float
        ra : float array
        beta : float
        rb : float array
        
        Returns
        -------
        c_exp : float
        c : float
        
        """
        c_exp, c = _mmd_ints.f90wrap_gaussian_product(rc=rc, alpha=alpha, ra=ra, \
            beta=beta, rb=rb)
        return c_exp, c
    
    @staticmethod
    def gaussian_normfactor(shell, alpha):
        """
        val = gaussian_normfactor(shell, alpha)
        
        
        Defined at gaussian_integral_utils.f90 lines 62-71
        
        Parameters
        ----------
        shell : int array
        alpha : float
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_gaussian_normfactor(shell=shell, alpha=alpha)
        return val
    
    @staticmethod
    def gaussian_coeff(l, m, xa, xb, k):
        """
        val = gaussian_coeff(l, m, xa, xb, k)
        
        
        Defined at gaussian_integral_utils.f90 lines 73-88
        
        Parameters
        ----------
        l : int
        m : int
        xa : float
        xb : float
        k : int
        
        Returns
        -------
        val : float
        
        """
        val = _mmd_ints.f90wrap_gaussian_coeff(l=l, m=m, xa=xa, xb=xb, k=k)
        return val
    
    _dt_array_initialisers = []
    

gaussian_integral_utils = Gaussian_Integral_Utils()

class Orbital_Types(f90wrap.runtime.FortranModule):
    """
    Module orbital_types
    
    
    Defined at orbital_types.f90 lines 1-7
    
    """
    @f90wrap.runtime.register_class("mmd_ints.basisfcn_t")
    class basisfcn_t(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=basisfcn_t)
        
        
        Defined at orbital_types.f90 lines 3-7
        
        """
        def __init__(self, handle=None):
            """
            self = Basisfcn_T()
            
            
            Defined at orbital_types.f90 lines 3-7
            
            
            Returns
            -------
            this : Basisfcn_T
            	Object to be constructed
            
            
            Automatically generated constructor for basisfcn_t
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _mmd_ints.f90wrap_basisfcn_t_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Basisfcn_T
            
            
            Defined at orbital_types.f90 lines 3-7
            
            Parameters
            ----------
            this : Basisfcn_T
            	Object to be destructed
            
            
            Automatically generated destructor for basisfcn_t
            """
            if self._alloc:
                _mmd_ints.f90wrap_basisfcn_t_finalise(this=self._handle)
        
        @property
        def shell(self):
            """
            Element shell ftype=integer pytype=int
            
            
            Defined at orbital_types.f90 line 4
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _mmd_ints.f90wrap_basisfcn_t__array__shell(self._handle)
            if array_handle in self._arrays:
                shell = self._arrays[array_handle]
            else:
                shell = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _mmd_ints.f90wrap_basisfcn_t__array__shell)
                self._arrays[array_handle] = shell
            return shell
        
        @shell.setter
        def shell(self, shell):
            self.shell[...] = shell
        
        @property
        def coeff(self):
            """
            Element coeff ftype=real pytype=float
            
            
            Defined at orbital_types.f90 line 5
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _mmd_ints.f90wrap_basisfcn_t__array__coeff(self._handle)
            if array_handle in self._arrays:
                coeff = self._arrays[array_handle]
            else:
                coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _mmd_ints.f90wrap_basisfcn_t__array__coeff)
                self._arrays[array_handle] = coeff
            return coeff
        
        @coeff.setter
        def coeff(self, coeff):
            self.coeff[...] = coeff
        
        @property
        def exps(self):
            """
            Element exps ftype=real pytype=float
            
            
            Defined at orbital_types.f90 line 6
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _mmd_ints.f90wrap_basisfcn_t__array__exps(self._handle)
            if array_handle in self._arrays:
                exps = self._arrays[array_handle]
            else:
                exps = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _mmd_ints.f90wrap_basisfcn_t__array__exps)
                self._arrays[array_handle] = exps
            return exps
        
        @exps.setter
        def exps(self, exps):
            self.exps[...] = exps
        
        @property
        def origin(self):
            """
            Element origin ftype=real pytype=float
            
            
            Defined at orbital_types.f90 line 7
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _mmd_ints.f90wrap_basisfcn_t__array__origin(self._handle)
            if array_handle in self._arrays:
                origin = self._arrays[array_handle]
            else:
                origin = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _mmd_ints.f90wrap_basisfcn_t__array__origin)
                self._arrays[array_handle] = origin
            return origin
        
        @origin.setter
        def origin(self, origin):
            self.origin[...] = origin
        
        def __str__(self):
            ret = ['<basisfcn_t>{\n']
            ret.append('    shell : ')
            ret.append(repr(self.shell))
            ret.append(',\n    coeff : ')
            ret.append(repr(self.coeff))
            ret.append(',\n    exps : ')
            ret.append(repr(self.exps))
            ret.append(',\n    origin : ')
            ret.append(repr(self.origin))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("mmd_ints.Basisfcn_T_X1:Num_Orbs_Array")
    class Basisfcn_T_X1:Num_Orbs_Array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=basisfcn_t_x1:num_orbs_array)
        
        
        Defined at orbital_types.f90 lines 3-7
        
        super-type
        Automatically generated to handle derived type arrays as a new derived type
        """
        def __init__(self, handle=None):
            """
            self = Basisfcn_T_X1:Num_Orbs_Array()
            
            
            Defined at orbital_types.f90 lines 3-7
            
            
            Returns
            -------
            this : Basisfcn_T_X1:Num_Orbs_Array
            	Object to be constructed
            
            
            Automatically generated constructor for basisfcn_t_x1:num_orbs_array
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _mmd_ints.f90wrap_basisfcn_t_x1:num_orbs_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Basisfcn_T_X1:Num_Orbs_Array
            
            
            Defined at orbital_types.f90 lines 3-7
            
            Parameters
            ----------
            this : Basisfcn_T_X1:Num_Orbs_Array
            	Object to be destructed
            
            
            Automatically generated destructor for basisfcn_t_x1:num_orbs_array
            """
            if self._alloc:
                _mmd_ints.f90wrap_basisfcn_t_x1:num_orbs_array_finalise(this=self._handle)
        
        def init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _mmd_ints.f90wrap_basisfcn_t_x1:num_orbs_array__array_getitem__items,
                                            _mmd_ints.f90wrap_basisfcn_t_x1:num_orbs_array__array_setitem__items,
                                            _mmd_ints.f90wrap_basisfcn_t_x1:num_orbs_array__array_len__items,
                                            """
            Element items ftype=type(basisfcn_t) pytype=Basisfcn_T
            
            
            Defined at  line 0
            
            """, orbital_types.basisfcn_t)
            return self.items
        
        _dt_array_initialisers = [init_array_items]
        
    
    _dt_array_initialisers = []
    

orbital_types = Orbital_Types()

class Special_Functions(f90wrap.runtime.FortranModule):
    """
    Module special_functions
    
    
    Defined at special_functions.f90 lines 1-99
    
    """
    @staticmethod
    def hyp1f1(a, b, x):
        """
        hg = hyp1f1(a, b, x)
        
        
        Defined at special_functions.f90 lines 5-99
        
        Parameters
        ----------
        a : float
        b : float
        x : float
        
        Returns
        -------
        hg : float
        
        """
        hg = _mmd_ints.f90wrap_hyp1f1(a=a, b=b, x=x)
        return hg
    
    _dt_array_initialisers = []
    

special_functions = Special_Functions()

class Constants(f90wrap.runtime.FortranModule):
    """
    Module constants
    
    
    Defined at constants.f90 lines 1-4
    
    """
    @property
    def pi(self):
        """
        Element pi ftype=real pytype=float
        
        
        Defined at constants.f90 line 3
        
        """
        return _mmd_ints.f90wrap_constants__get__pi()
    
    @property
    def hbar(self):
        """
        Element hbar ftype=real(kind=8) pytype=float
        
        
        Defined at constants.f90 line 4
        
        """
        return _mmd_ints.f90wrap_constants__get__hbar()
    
    @property
    def e(self):
        """
        Element e ftype=real pytype=float
        
        
        Defined at constants.f90 line 5
        
        """
        return _mmd_ints.f90wrap_constants__get__e()
    
    def __str__(self):
        ret = ['<constants>{\n']
        ret.append('    pi : ')
        ret.append(repr(self.pi))
        ret.append(',\n    hbar : ')
        ret.append(repr(self.hbar))
        ret.append(',\n    e : ')
        ret.append(repr(self.e))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

constants = Constants()

