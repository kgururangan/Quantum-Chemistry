import numpy as np
import pprint

from scipy.special import factorial2 as fact2

import basis_set_exchange as bse

ATOMIC_NUMBER = {
        'H': 1,
        'He': 2,
        'C': 6,
        'N': 7,
        'O': 8,
        }

class Molecule:
    """
    A class that keeps information for a molecular system

    """

    def __init__(self,
            geometry=None,
            basis=None,
            angs=False
            ):

        if angs:
            self.convert = 1.8897259886
        else:
            self.convert = 1.0

        # get geometry
        self.geometry = None
        self.nel = None
        self.parse_molecule(geometry)

        # get basis set from BSE
        self.basis = basis
        self.basis_set = None
        self.get_basis()

        # build orbitals
        self.orbitals = None
        self.build_orbitals()

        self.norb = len(self.orbitals)

        # compute nuclear repulsion energy
        self.e_nn = None
        self.nuclear_repulsion()

    def parse_molecule(self, geometry):
        """
        Parse molecular information from string
        """

        geom = []
        nel = 0
        if isinstance(geometry, str):
            for line in geometry.split("\n"):
                fields = line.split()
                x = np.float64(fields[1]) * self.convert
                y = np.float64(fields[2]) * self.convert
                z = np.float64(fields[3]) * self.convert

                Z = ATOMIC_NUMBER[fields[0]]
                nel += Z

                geom.append({
                    'label': fields[0],
                    'Z': Z,
                    'x': x,
                    'y': y,
                    'z': z,
                    'r': (x,y,z)
                    })

        self.geometry = geom
        self.nel = nel

    def build_orbitals(self):

        orbitals = []
        cnt = 0
        for atom in self.geometry:
            for shell in self.basis_set[atom['label']]:
                for momentum in shell['shells']:
                    cnt += 1
                    orbitals.append(
                            BasisFunction(
                                origin = [atom['x'], atom['y'],
                                    atom['z']],
                                shell = momentum,
                                exps = shell['exps'],
                                coefs = shell['coefs']
                                )
                            )

        self.orbitals = orbitals


    def get_basis(self):

        atom_types = set([i['label'] for i in self.geometry])
        zs = [ATOMIC_NUMBER[i] for i in atom_types]

        bs = bse.get_basis(self.basis, elements=zs)
        self.bs = bs

        basis_set = dict()
        orbitals = []

        for element in atom_types:

            basis_set[element] = []

            bs_data = bs['elements'][str(ATOMIC_NUMBER[element])]['electron_shells']

            for shell in bs_data:


                exps = list(map(np.float64,
                    shell['exponents']))

                for i, momentum in \
                        enumerate(shell['angular_momentum']):
                    shells = []
                    if momentum == 0:
                        shells.append((0,0,0))
                    elif momentum == 1:
                        shells.append((1,0,0))
                        shells.append((0,1,0))
                        shells.append((0,0,1))
                    elif momentum == 2:
                        shells.append((2,0,0))
                        shells.append((0,2,0))
                        shells.append((0,0,2))
                        shells.append((1,1,0))
                        shells.append((1,0,1))
                        shells.append((0,1,1))


                    coefs = list(map(
                        np.float64,shell['coefficients'][i]))

                    basis_set[element].append({
                            'shells': shells,
                            'exps': exps,
                            'coefs': coefs,
                            })

        self.basis_set = basis_set


    def nuclear_repulsion(self):

        e_nn = 0.0
        geom = self.geometry
        for i in range(len(geom)):
            Zi = geom[i]['Z']
            Ri = np.array(geom[i]['r'])
            for j in range(i+1,len(geom)):
                Zj = geom[j]['Z']
                Rj = np.array(geom[j]['r'])

                r = np.linalg.norm(Ri - Rj)

                e_nn += Zi * Zj / r

        self.e_nn = e_nn


class BasisFunction:

    ''' A class that contains all our basis function data
        Attributes:
        origin: array/list containing the coordinates of the Gaussian origin
        shell:  tuple of angular momentum
        exps:   list of primitive Gaussian exponents
        coefs:  list of primitive Gaussian coefficients
        norm:   list of normalization factors for Gaussian primitives
    '''

    def __init__(self,
            origin=[0.0,0.0,0.0],
            shell=(0,0,0),
            exps=[],
            coefs=[]
            ):

        self.origin = np.asarray(origin)
        self.shell = shell
        self.exps  = exps
        self.coefs = coefs
        self.norm = None
        self.normalize()

    def normalize(self):

        ''' Routine to normalize the basis functions, in case they
            do not integrate to unity.
        '''

        l, m, n = self.shell
        L = l+m+n
        # self.norm is a list of length equal to number primitives
        # normalize primitives first (PGBFs)
        self.norm = np.sqrt(np.power(2,2*(l+m+n)+1.5)*
                        np.power(self.exps,l+m+n+1.5)/
                        fact2(2*l-1)/fact2(2*m-1)/
                        fact2(2*n-1)/np.power(np.pi,1.5))

        # now normalize the contracted basis functions (CGBFs)
        # Eq. 1.44 of Valeev integral whitepaper
        prefactor = np.power(np.pi,1.5)*\
            fact2(2*l - 1)*fact2(2*m - 1)*fact2(2*n - 1)/np.power(2.0,L)

        N = 0.0
        num_exps = len(self.exps)
        for ia in range(num_exps):
            for ib in range(num_exps):
                N += self.norm[ia]*self.norm[ib]*self.coefs[ia]*self.coefs[ib]/\
                         np.power(self.exps[ia] + self.exps[ib],L+1.5)

        N *= prefactor
        N = np.power(N,-0.5)
        for ia in range(num_exps):
            self.coefs[ia] *= N

    def show(self):
        print("Origin: ", self.origin)
        print("Shell: ", self.shell)
        print("Exponents: ", self.exps)
        print("Coefficients: ", self.coefs)
        print("Norm: ", self.norm)


