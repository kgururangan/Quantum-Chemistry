import numpy as np

import numba as nb
from numba import jit

from scipy.special import hyp1f1

# helper functions
def boys(n, T):
    return hyp1f1(n+0.5, n+1.5, -T) / (2.0 * n + 1.0)

def gaussian_product_center(a, A, b, B):
    return (a*A + b*B) / (a + b)


@nb.njit(
        nb.float64(
            nb.int64,
            nb.int64,
            nb.int64,
            nb.float64,
            nb.float64,
            nb.float64,
            )
        )
def E(i,j,t,Qx,a,b):

    ''' Recursive definition of Hermite Gaussian coefficients.
        Returns a float.
        a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
        i,j: orbital angular momentum number on Gaussian 'a' and 'b'
        t: number nodes in Hermite (depends on type of integral,
           e.g. always zero for overlap integrals)
        Qx: distance between origins of Gaussian 'a' and 'b'
    '''

    p = a + b
    q = a*b/p
    if (t < 0) or (t > (i + j)):
        # out of bounds for t
        return 0.0
    elif i == j == t == 0:
        # base case
        return np.exp(-q*Qx*Qx) # K_AB
    elif j == 0:
        # decrement index i
        return (1/(2*p))*E(i-1,j,t-1,Qx,a,b) - \
                (q*Qx/a)*E(i-1,j,t,Qx,a,b)    + \
                (t+1)*E(i-1,j,t+1,Qx,a,b)
    else:
        # decrement index j
        return (1/(2*p))*E(i,j-1,t-1,Qx,a,b) + \
                (q*Qx/b)*E(i,j-1,t,Qx,a,b)    + \
                (t+1)*E(i,j-1,t+1,Qx,a,b)

@nb.njit(
        nb.float64(
            nb.float64,
            nb.types.UniTuple(nb.int64, 3),
            nb.float64[::1],
            nb.float64,
            nb.types.UniTuple(nb.int64, 3),
            nb.float64[::1]
            )
        )
def overlap(a, lmn1, A, b, lmn2, B):

    ''' Evaluates overlap integral between two Gaussians
        Returns a float.
        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
              for Gaussian 'a'
        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
        B:    list containing origin of Gaussian 'b'
    '''

    l1,m1,n1 = lmn1 # shell angular momentum on Gaussian 'a'
    l2,m2,n2 = lmn2 # shell angular momentum on Gaussian 'b'

    #import pdb
    #pdb.set_trace()

    S1 = E(l1,l2,0,A[0]-B[0],a,b) # X
    S2 = E(m1,m2,0,A[1]-B[1],a,b) # Y
    S3 = E(n1,n2,0,A[2]-B[2],a,b) # Z

    return S1*S2*S3*np.power(np.pi/(a+b),1.5)



def S(a,b):

    '''Evaluates overlap between two contracted Gaussians
       Returns float.
       Arguments:
       a: contracted Gaussian 'a', BasisFunction object
       b: contracted Gaussian 'b', BasisFunction object
    '''

    s = 0.0
    for ia, ca in enumerate(a.coefs):
        for ib, cb in enumerate(b.coefs):
            s += (
                    a.norm[ia]
                    * b.norm[ib]
                    * ca
                    * cb
                    * overlap(a.exps[ia], a.shell, a.origin,
                     b.exps[ib], b.shell, b.origin)
                    )

    return s


@nb.njit
def kinetic(a, lmn1, A, b, lmn2, B):

    ''' Evaluates kinetic energy integral between two Gaussians
        Returns a float.
        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
              for Gaussian 'a'
        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
        B:    list containing origin of Gaussian 'b'
    '''

    l1,m1,n1 = lmn1
    l2,m2,n2 = lmn2
    term0 = b*(2*(l2+m2+n2)+3)*\
                            overlap(a,(l1,m1,n1),A,b,(l2,m2,n2),B)
    term1 = -2*np.power(b,2)*\
                           (overlap(a,(l1,m1,n1),A,b,(l2+2,m2,n2),B) +
                            overlap(a,(l1,m1,n1),A,b,(l2,m2+2,n2),B) +
                            overlap(a,(l1,m1,n1),A,b,(l2,m2,n2+2),B))
    term2 = -0.5*(l2*(l2-1)*overlap(a,(l1,m1,n1),A,b,(l2-2,m2,n2),B) +
                  m2*(m2-1)*overlap(a,(l1,m1,n1),A,b,(l2,m2-2,n2),B) +
                  n2*(n2-1)*overlap(a,(l1,m1,n1),A,b,(l2,m2,n2-2),B))

    return term0+term1+term2


def T(a,b):

    '''Evaluates kinetic energy between two contracted Gaussians
       Returns float.
       Arguments:
       a: contracted Gaussian 'a', BasisFunction object
       b: contracted Gaussian 'b', BasisFunction object
    '''

    t = 0.0
    for ia, ca in enumerate(a.coefs):
        for ib, cb in enumerate(b.coefs):
            t += a.norm[ia]*b.norm[ib]*ca*cb*\
                     kinetic(a.exps[ia],a.shell,a.origin,\
                     b.exps[ib],b.shell,b.origin)
    return t



def R(t,u,v,n,p,PCx,PCy,PCz,RPC):

    ''' Returns the Coulomb auxiliary Hermite integrals
        Returns a float.
        Arguments:
        t,u,v:   order of Coulomb Hermite derivative in x,y,z
                 (see defs in Helgaker and Taylor)
        n:       order of Boys function
        PCx,y,z: Cartesian vector distance between Gaussian
                 composite center P and nuclear center C
        RPC:     Distance between P and C
    '''

    T = p*RPC*RPC
    val = 0.0
    if t == u == v == 0:
        val += np.power(-2*p,n)*boys(n,T)
    elif t == u == 0:
        if v > 1:
            val += (v-1)*R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC)
        val += PCz*R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)
    elif t == 0:
        if u > 1:
            val += (u-1)*R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC)
        val += PCy*R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)
    else:
        if t > 1:
            val += (t-1)*R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC)
        val += PCx*R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)
    return val



def nuclear_attraction(a,lmn1,A,b,lmn2,B,C):

    ''' Evaluates kinetic energy integral between two Gaussians
         Returns a float.
         a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
         b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
         lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
               for Gaussian 'a'
         lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
         A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
         B:    list containing origin of Gaussian 'b'
         C:    list containing origin of nuclear center 'C'
    '''

    l1,m1,n1 = lmn1
    l2,m2,n2 = lmn2
    p = a + b
    P = gaussian_product_center(a,A,b,B) # Gaussian composite center
    RPC = np.linalg.norm(P-C)

    val = 0.0
    for t in range(l1+l2+1):
        for u in range(m1+m2+1):
            for v in range(n1+n2+1):
                val += E(l1,l2,t,A[0]-B[0],a,b) * \
                       E(m1,m2,u,A[1]-B[1],a,b) * \
                       E(n1,n2,v,A[2]-B[2],a,b) * \
                       R(t,u,v,0,p,P[0]-C[0],P[1]-C[1],P[2]-C[2],RPC)
    val *= 2*np.pi/p

    return val


def V(a,b,C):

    '''Evaluates overlap between two contracted Gaussians
       Returns float.
       Arguments:
       a: contracted Gaussian 'a', BasisFunction object
       b: contracted Gaussian 'b', BasisFunction object
       C: center of nucleus
    '''

    v = 0.0
    for ia, ca in enumerate(a.coefs):
        for ib, cb in enumerate(b.coefs):
            v += a.norm[ia]*b.norm[ib]*ca*cb*\
                     nuclear_attraction(a.exps[ia],a.shell,a.origin,
                     b.exps[ib],b.shell,b.origin,C)
    return v


def electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,lmn4,D):

    ''' Evaluates electronic repulsioni energy integral between two Gaussians
        Returns a float.
        a,b,c,d:   orbital exponent on Gaussian 'a','b','c','d'
        lmn1,lmn2
        lmn3,lmn4: int tuple containing orbital angular momentum
                   for Gaussian 'a','b','c','d', respectively
        A,B,C,D:   list containing origin of Gaussian 'a','b','c','d'
    '''

    l1,m1,n1 = lmn1
    l2,m2,n2 = lmn2
    l3,m3,n3 = lmn3
    l4,m4,n4 = lmn4
    p = a+b # composite exponent for P (from Gaussians 'a' and 'b')
    q = c+d # composite exponent for Q (from Gaussians 'c' and 'd')
    alpha = p*q/(p+q)
    P = gaussian_product_center(a,A,b,B) # A and B composite center
    Q = gaussian_product_center(c,C,d,D) # C and D composite center
    RPQ = np.linalg.norm(P-Q)

    val = 0.0
    for t in range(l1+l2+1):
        for u in range(m1+m2+1):
            for v in range(n1+n2+1):
                for tau in range(l3+l4+1):
                    for nu in range(m3+m4+1):
                        for phi in range(n3+n4+1):
                            val += E(l1,l2,t,A[0]-B[0],a,b) * \
                                   E(m1,m2,u,A[1]-B[1],a,b) * \
                                   E(n1,n2,v,A[2]-B[2],a,b) * \
                                   E(l3,l4,tau,C[0]-D[0],c,d) * \
                                   E(m3,m4,nu ,C[1]-D[1],c,d) * \
                                   E(n3,n4,phi,C[2]-D[2],c,d) * \
                                   np.power(-1,tau+nu+phi) * \
                                   R(t+tau,u+nu,v+phi,0,\
                                       alpha,P[0]-Q[0],P[1]-Q[1],P[2]-Q[2],RPQ)

    val *= 2*np.power(np.pi,2.5)/(p*q*np.sqrt(p+q))

    return val

def ERI(a,b,c,d):

    '''Evaluates overlap between two contracted Gaussians
        Returns float.
        Arguments:
        a: contracted Gaussian 'a', BasisFunction object
        b: contracted Gaussian 'b', BasisFunction object
        c: contracted Gaussian 'b', BasisFunction object
        d: contracted Gaussian 'b', BasisFunction object
    '''

    eri = 0.0
    for ja, ca in enumerate(a.coefs):
        for jb, cb in enumerate(b.coefs):
            for jc, cc in enumerate(c.coefs):
                for jd, cd in enumerate(d.coefs):
                    eri += a.norm[ja]*b.norm[jb]*c.norm[jc]*d.norm[jd]*\
                             ca*cb*cc*cd*\
                             electron_repulsion(a.exps[ja],a.shell,a.origin,\
                                                b.exps[jb],b.shell,b.origin,\
                                                c.exps[jc],c.shell,c.origin,\
                                                d.exps[jd],d.shell,d.origin)
    return eri

def build_onebody(mol):

    norb = mol.norb

    Smat = np.zeros((norb, norb))
    Tmat = np.zeros((norb, norb))
    Amat = np.zeros((norb, norb))

    for i in range(norb):
        iorb = mol.orbitals[i]
        for j in range(i,norb):
            jorb = mol.orbitals[j]
            val = S(iorb,jorb)
            Smat[i,j] = val
            Smat[j,i] = val

            val = T(iorb,jorb)
            Tmat[i,j] = val
            Tmat[j,i] = val

            val = 0.0
            for k in mol.geometry:
                val -= V(iorb, jorb, k['r']) * k['Z']
            Amat[i,j] = val
            Amat[j,i] = val

    return Smat, Tmat + Amat

def build_twobody(mol):

    norb = mol.norb

    Vmat = np.zeros((norb, norb, norb, norb))
    # Calculate twobody integrals
    for i in range(norb):
        iorb = mol.orbitals[i]
        for j in range(i+1):
            jorb = mol.orbitals[j]
            ij = i*(i+1)/2 + j
            for k in range(norb):
                korb = mol.orbitals[k]
                for l in range(k+1):
                    lorb = mol.orbitals[l]
                    kl = k*(k+1)//2 + l
                    if ij >= kl:
                        val = ERI(iorb,jorb,korb,lorb)
                        Vmat[i,k,j,l] = val
                        Vmat[j,k,i,l] = val
                        Vmat[i,l,j,k] = val
                        Vmat[k,i,l,j] = val
                        Vmat[l,i,k,j] = val
                        Vmat[l,j,k,i] = val
                        Vmat[j,l,i,k] = val
                        Vmat[k,j,l,i] = val
    return Vmat
