import click

import numpy as np
import time

from basis import BasisFunction, Molecule
from integrals import build_onebody, build_twobody
from utils import orthomat, diis_pulay_solver

import sys

WATER = """H -0.0000000000       -1.5152630000       -1.0588980000
H 0.0000000000        1.5152630000       -1.0588980000
O 0.0000000000        0.0000000000       -0.0090000000"""


@click.command()
def main():

    mol = Molecule(geometry=WATER, basis='STO-3G', angs=False)
    nocc = int(mol.nel / 2)
    norb = mol.norb

    import pdb
    pdb.set_trace()

    # build integrals
    print("Building integrals")
    start = time.time()
    Smat, Hcore = build_onebody(mol)
    #Vmat = build_twobody(mol)
    #np.save("Vmat.npy", Vmat)
    #sys.exit()
    Vmat = np.load("Vmat.npy")
    end = time.time()
    print("Took {:.3f}s".format(end-start))

    X = orthomat(Smat, 0.0, 'symmetric')

    P = np.zeros((norb, norb))

    G = np.einsum('prqs,rs->pq', Vmat, P, optimize=True) \
            - 0.5*np.einsum('prsq,rs->pq', Vmat, P, optimize=True)

    F = Hcore + G

    eps_mo, Cp = np.linalg.eigh(np.dot(X.T, np.dot(F, X)))
    idx = eps_mo.argsort()
    eps_mo = eps_mo[idx]
    C = np.dot(X, Cp[:,idx])
    P = 2*np.einsum('pi,qi->pq', C[:,:nocc], C[:,:nocc], optimize=True)

    # DIIS containers
    F_list = []
    resid_F = []


    for it in range(30):

        G = np.einsum('prqs,rs->pq', Vmat, P, optimize=True) \
                - 0.5*np.einsum('prsq,rs->pq', Vmat, P, optimize=True)
        F = Hcore + G

        # Build DIIS Residual
        diis_r = X.dot(F.dot(P).dot(Smat) - Smat.dot(P).dot(F)).dot(X)
        scf_resid = np.mean(diis_r**2)**0.5

        #if scf_resid <= tol_scf:
        #    flag_conv_scf = True
        #    it_scf = it
        #    break

        # Append trial & residual vectors to lists
        # (Note: do NOT use extrapolated matrices in F_list!)
        if len(F_list) == 3:
            del F_list[0]
            del resid_F[0]
        F_list.append(F)
        resid_F.append(diis_r)

        # DIIS extrapolated Fock matrix
        if it >= 2:
            F = diis_pulay_solver(F_list, resid_F)

        eps_mo, Cp = np.linalg.eigh(np.dot(X.T, np.dot(F, X)))
        idx = eps_mo.argsort()
        eps_mo = eps_mo[idx]
        C = np.dot(X, Cp[:,idx])
        P = 2*np.einsum('pi,qi->pq', C[:,:nocc], C[:,:nocc], optimize=True)

        H00 = 0.5*np.einsum('qp,pq->', P, Hcore + F,
                optimize=True) + mol.e_nn
        print(it, H00)







    #eri = ERI()


if __name__ == "__main__":
    main()
