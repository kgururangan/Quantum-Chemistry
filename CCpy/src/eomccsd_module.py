import numpy as np
from HBar_module import HBar_CCSD

#from fortran_cis import cis_hamiltonian
#print(cis_hamiltonian.__doc__)


def eomccsd(cc_t,ints,sys):

    # Number of roots
    nroot = 1
    # Number of guess vectors
    guess_size = 4
    # Convergence tolerance of the residual norm
    ctol = 1.e-7
    # Convergence tolerance of the energy
    etol = 1.e-9

    # Build the CCSD similarity-transformed Hamiltonian
    H1A, H1B, H2A, H2B, H2C = HBar_CCSD(cc_t,ints,sys)
    
    # Obtain the vectors of a CIS calculation
    gevecs, gevals = cis(nroot,ints,sys)

    # Maximum number of vectors
    max_guess = 200

    # Current number of vectors
    num_vecs = guess_size

    # Copy gvec data into in ci_gvecs
    arr_cvecs = np.asarray(cvecs)
    for x in range(guess_size):
        arr_cvecs[:] = gvecs[x]
        cvecs.write(x, 0)
        cvecs.symnormalize(1 / np.linalg.norm(gvecs[x]), x)

    delta_c = np.zeros(nroot)

    G = np.zeros((max_guess, max_guess))

    # Begin Davidson iterations
    for CI_ITER in range(max_guess - 1):

        # Subspace Matrix, Gij = < bi | H | bj >
        for i in range(0, num_vecs):
            # Build sigma for each b
            cvecs.read(i, 0)
            svecs.read(i, 0)
            ciwfn.sigma(cvecs, svecs, i, i)
            for j in range(i, num_vecs):
                # G_ij = (b_i, sigma_j)
                cvecs.read(i, 0)
                svecs.read(j, 0)
                G[j, i] = G[i, j] = svecs.vdot(cvecs, i, j)

        evals, evecs = np.linalg.eigh(G[:num_vecs, :num_vecs])
        CI_E = evals

        # Use average over roots as convergence criteria
        avg_energy = 0.0
        avg_dc = 0.0
        for n in range(nroot):
            avg_energy += evals[n]
            avg_dc += delta_c[n]
        avg_energy /= nroot
        avg_dc /= nroot
        avg_energy += mol.nuclear_repulsion_energy()

        print('CI Iteration %3d: Energy = %4.16f   dE = % 1.5E   dC = %1.5E' % (CI_ITER, avg_energy, (avg_energy - Eold),
                                                                                avg_dc))
        if (abs(avg_energy - Eold) < etol) and (avg_dc < ctol) and (CI_ITER > 3):
            print('CI has converged!\n')
            break
        Eold = avg_energy

        # Build new vectors as linear combinations of the subspace matrix, H
        for n in range(nroot):

            # Build as linear combinations of previous vectors
            dvecs.zero()
            dvecs.write(dwork_vec, 0)
            for c in range(len(evecs[:, n])):
                dvecs.axpy(evecs[c, n], cvecs, dwork_vec, c)

            # Build new vector new_vec = ((H * cvec) - evals[n] * cvec) / (evals[n] - Hd)
            ciwfn.sigma(dvecs, svecs, dwork_vec, swork_vec)
            svecs.axpy(-1 * evals[n], dvecs, swork_vec, dwork_vec)
            norm = svecs.dcalc(evals[n], Hd, swork_vec)

            if (norm < 1e-9):
                continue

            svecs.symnormalize(1 / norm, swork_vec)
            delta_c[n] = norm

            # Build a new vector that is orthornormal to all previous vectors
            dvecs.copy(svecs, n, swork_vec)
            norm = dvecs.norm(n)
            dvecs.symnormalize(1 / norm, n)

            total_proj = 0
            for i in range(num_vecs):
                proj = svecs.vdot(cvecs, swork_vec, i)
                total_proj += proj
                dvecs.axpy(-proj, cvecs, n, i)

            norm = dvecs.norm(n)
            dvecs.symnormalize(1 / norm, n)

            # This *should* screen out contributions that are projected out by above
            if True:
                cvecs.write(num_vecs, 0)
                cvecs.copy(dvecs, num_vecs, n)
                num_vecs += 1

    return 


def flatten_R(r1a,r1b,r2a,r2b,r2c):
    return np.concatenate((r1a.flatten(),r1b.flatten(),r2a.flatten(),r2b.flatten(),r2c.flatten()),axis=0)

def unflatten_R(R,sys):

    n1a = sys['Nocc_a'] * sys['Nunocc_a']
    n1b = sys['Nocc_b'] * sys['Nunocc_b']
    n2a = sys['Nocc_a'] ** 2 * sys['Nunocc_a'] ** 2
    n2b = sys['Nocc_a'] * sys['Nocc_b'] * sys['Nunocc_a'] * sys['Nunocc_b']
    n2c = sys['Nocc_b'] ** 2 * sys['Nunocc_b'] ** 2
    idx_1a = slice(0,n1a)
    idx_1b = slice(n1a,n1a+n1b)
    idx_2a = slice(n1a+n1b,n1a+n1b+n2a)
    idx_2b = slice(n1a+n1b+n2a,n1a+n1b+n2a+n2b)
    idx_2c = slice(n1a+n1b+n2a+n2b,n1a+n1b+n2a+n2b+n2c)

    r1a  = np.reshape(R[idx_1a],(sys['Nunocc_a'],sys['Nocc_a']))
    r1b  = np.reshape(R[idx_1b],(sys['Nunocc_b'],sys['Nocc_b']))
    r2a  = np.reshape(R[idx_2a],(sys['Nunocc_a'],sys['Nunocc_a'],sys['Nocc_a'],sys['Nocc_a']))
    r2b  = np.reshape(R[idx_2b],(sys['Nunocc_a'],sys['Nunocc_b'],sys['Nocc_a'],sys['Nocc_b']))
    r2c  = np.reshape(R[idx_2c],(sys['Nunocc_b'],sys['Nunocc_b'],sys['Nocc_b'],sys['Nocc_b']))

    return r1a, r1b, r2a, r2b, r2c


def HR(R,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys):

    r1a, r1b, r2a, r2b, r2c = unflatten_R(R,sys)

    X1A = build_HR_1A(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    X1B = build_HR_1B(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    X2A = build_HR_2A(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    X2B = build_HR_2B(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    X2C = build_HR_2C(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)

    return flatten_R(X1A, X1B, X2A, X2B, X2C)

def build_HR_1A(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys):

    X1A = 0.0
    X1A -= np.einsum('mi,am->ai',H1A['oo'],r1a,optimize=True)
    X1A += np.einsum('ae,ei->ai',H1A['vv'],r1a,optimize=True)
    X1A += np.einsum('amie,em->ai',H2A['voov'],r1a,optimize=True)
    X1A += np.einsum('amie,em->ai',H2B['voov'],r1b,optimize=True)
    X1A  -= 0.5*np.einsum('mnif,afmn->ai',H2A['ooov'],r2a,optimize=True)
    X1A -= np.einsum('mnif,afmn->ai',H2B['ooov'],r2b,optimize=True)
    X1A += 0.5*np.einsum('anef,efin->ai',H2A['vovv'],r2a,optimize=True)
    X1A += np.einsum('anef,efin->ai',H2B['vovv'],r2b,optimize=True)
    X1A += np.einsum('me,aeim->ai',H1A['ov'],r2a,optimize=True)
    X1A += np.einsum('me,aeim->ai',H1B['ov'],r2b,optimize=True)

    return X1A

def build_HR_1B(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys):

    X1B = 0.0
    X1B -= np.einsum('mi,am->ai',H1B['oo'],r1b,optimize=True)
    X1B += np.einsum('ae,ei->ai',H1B['vv'],r1b,optimize=True)
    X1B += np.einsum('maei,em->ai',H2B['ovvo'],r1a,optimize=True)
    X1B += np.einsum('amie,em->ai',H2C['voov'],r1b,optimize=True)
    X1B -= np.einsum('nmfi,fanm->ai',H2B['oovo'],r2b,optimize=True)
    X1B -= 0.5*np.einsum('mnif,afmn->ai',H2C['ooov'],r2c,optimize=True)
    X1B += np.einsum('nafe,feni->ai',H2B['ovvv'],r2b,optimize=True)
    X1B += 0.5*np.einsum('anef,efin->ai',H2C['vovv'],r2c,optimize=True)
    X1B += np.einsum('me,eami->ai',H1A['ov'],r2b,optimize=True)
    X1B += np.einsum('me,aeim->ai',H1B['ov'],r2c,optimize=True)

    return X1B

def build_HR_2A(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys):

    t2a = cc_t['t2a']
    vA = ints['vA']
    vB = ints['vB']

    X2A = 0.0
    D1 = -np.einsum('mi,abmj->abij',H1A['oo'],r2a,optimize=True) # A(ij) 
    D2 = np.einsum('ae,ebij->abij',H1A['vv'],r2a,optimize=True) # A(ab)
    X2A += 0.5*np.einsum('mnij,abmn->abij',H2A['oooo'],r2a,optimize=True)
    X2A += 0.5*np.einsum('abef,efij->abij',H2A['vvvv'],r2a,optimize=True)
    D3 = np.einsum('amie,ebmj->abij',H2A['voov'],r2a,optimize=True) # A(ij)A(ab)
    D4 = np.einsum('amie,bejm->abij',H2B['voov'],r2b,optimize=True) # A(ij)A(ab)
    D5 = -np.einsum('bmji,am->abij',H2A['vooo'],r1a,optimize=True) # A(ab)
    D6 = np.einsum('baje,ei->abij',H2A['vvov'],r1a,optimize=True) # A(ij)

    Q1 = -0.5*np.einsum('mnef,bfmn->eb',vA['oovv'],r2a,optimize=True)
    D7 = np.einsum('eb,aeij->abij',Q1,t2a,optimize=True) # A(ab)
    Q2 = -np.einsum('mnef,bfmn->eb',vB['oovv'],r2b,optimize=True)
    D8 = np.einsum('eb,aeij->abij',Q2,t2a,optimize=True) # A(ab)

    Q1 = 0.5*np.einsum('mnef,efjn->mj',vA['oovv'],r2a,optimize=True)
    D9 = -np.einsum('mj,abim->abij',Q1,t2a,optimize=True) # A(ij)
    Q2 = np.einsum('mnef,efjn->mj',vB['oovv'],r2b,optimize=True)
    D10 = -np.einsum('mj,abim->abij',Q2,t2a,optimize=True) # A(ij)

    Q1 = np.einsum('amfe,em->af',H2A['vovv'],r1a,optimize=True)
    D11 = np.einsum('af,fbij->abij',Q1,t2a,optimize=True) # A(ab)
    Q2 = np.einsum('nmie,em->ni',H2A['ooov'],r1a,optimize=True)
    D12 = -np.einsum('ni,abnj->abij',Q2,t2a,optimize=True) # A(ij)

    Q1 = np.einsum('amfe,em->af',H2B['vovv'],r1b,optimize=True)
    D13 = np.einsum('af,fbij->abij',Q1,t2a,optimize=True) # A(ab)
    Q2 = np.einsum('nmie,em->ni',H2B['ooov'],r1b,optimize=True)
    D14 = -np.einsum('ni,abnj->abij',Q2,t2a,optimize=True) # A(ij)

    D_ij = D1 + D6 + D9 + D10 + D12 + D14
    D_ab = D2 + D5 + D7 + D8  + D11 + D13
    D_abij = D3 + D4

    D_ij -= np.einsum('abij->abji',D_ij,optimize=True)
    D_ab -= np.einsum('abij->baij',D_ab,optimize=True)
    D_abij += -np.einsum('abij->baij',D_abij,optimize=True)\
    -np.einsum('abij->abji',D_abij,optimize=True)\
    +np.einsum('abij->baji',D_abij,optimize=True)    

    X2A += D_ij + D_ab + D_abij

    return X2A

def build_HR_2B(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys):

    t2b = cc_t['t2b']
    vA = ints['vA']
    vB = ints['vB']
    vC = ints['vC']

    X2B = 0.0
    X2B += np.einsum('ae,ebij->abij',H1A['vv'],r2b,optimize=True)
    X2B += np.einsum('be,aeij->abij',H1B['vv'],r2b,optimize=True)
    X2B -= np.einsum('mi,abmj->abij',H1A['oo'],r2b,optimize=True)
    X2B -= np.einsum('mj,abim->abij',H1B['oo'],r2b,optimize=True)
    X2B += np.einsum('mnij,abmn->abij',H2B['oooo'],r2b,optimize=True)
    X2B += np.einsum('abef,efij->abij',H2B['vvvv'],r2b,optimize=True)
    X2B += np.einsum('amie,ebmj->abij',H2A['voov'],r2b,optimize=True)
    X2B += np.einsum('amie,ebmj->abij',H2B['voov'],r2c,optimize=True)
    X2B += np.einsum('mbej,aeim->abij',H2B['ovvo'],r2a,optimize=True)
    X2B += np.einsum('bmje,aeim->abij',H2C['voov'],r2b,optimize=True)
    X2B -= np.einsum('mbie,aemj->abij',H2B['ovov'],r2b,optimize=True)
    X2B -= np.einsum('amej,ebim->abij',H2B['vovo'],r2b,optimize=True)
    X2B += np.einsum('abej,ei->abij',H2B['vvvo'],r1a,optimize=True)
    X2B += np.einsum('abie,ej->abij',H2B['vvov'],r1b,optimize=True)
    X2B -= np.einsum('mbij,am->abij',H2B['ovoo'],r1a,optimize=True)
    X2B -= np.einsum('amij,bm->abij',H2B['vooo'],r1b,optimize=True)

    Q1 = -0.5*np.einsum('mnef,afmn->ae',vA['oovv'],r2a,optimize=True)
    X2B += np.einsum('ae,ebij->abij',Q1,t2b,optimize=True)
    Q2 = 0.5*np.einsum('mnef,efin->mi',vA['oovv'],r2a,optimize=True)
    X2B -= np.einsum('mi,abmj->abij',Q2,t2b,optimize=True)

    Q1 = -np.einsum('nmfe,fbnm->be',vB['oovv'],r2b,optimize=True)
    X2B += np.einsum('be,aeij->abij',Q1,t2b,optimize=True)
    Q2 = -np.einsum('mnef,afmn->ae',vB['oovv'],r2b,optimize=True)
    X2B += np.einsum('ae,ebij->abij',Q2,t2b,optimize=True)
    Q3 = np.einsum('nmfe,fenj->mj',vB['oovv'],r2b,optimize=True)
    X2B -= np.einsum('mj,abim->abij',Q3,t2b,optimize=True)
    Q4 = np.einsum('mnef,efin->mi',vB['oovv'],r2b,optimize=True)
    X2B -= np.einsum('mi,abmj->abij',Q4,t2b,optimize=True)

    Q1 = -0.5*np.einsum('mnef,bfmn->be',vC['oovv'],r2c,optimize=True)
    X2B += np.einsum('be,aeij->abij',Q1,t2b,optimize=True)
    Q2 = 0.5*np.einsum('mnef,efjn->mj',vC['oovv'],r2c,optimize=True)
    X2B -= np.einsum('mj,abim->abij',Q2,t2b,optimize=True)

    Q1 = np.einsum('mbef,em->bf',H2B['ovvv'],r1a,optimize=True)
    X2B += np.einsum('bf,afij->abij',Q1,t2b,optimize=True)
    Q2 = np.einsum('mnej,em->nj',H2B['oovo'],r1a,optimize=True)
    X2B -= np.einsum('nj,abin->abij',Q2,t2b,optimize=True)
    Q3 = np.einsum('amfe,em->af',H2A['vovv'],r1a,optimize=True)
    X2B += np.einsum('af,fbij->abij',Q3,t2b,optimize=True)
    Q4 = np.einsum('nmie,em->ni',H2A['ooov'],r1a,optimize=True)
    X2B -= np.einsum('ni,abnj->abij',Q4,t2b,optimize=True)

    Q1 = np.einsum('amfe,em->af',H2B['vovv'],r1b,optimize=True)
    X2B += np.einsum('af,fbij->abij',Q1,t2b,optimize=True)
    Q2 = np.einsum('nmie,em->ni',H2B['ooov'],r1b,optimize=True)
    X2B -= np.einsum('ni,abnj->abij',Q2,t2b,optimize=True)
    Q3 = np.einsum('bmfe,em->bf',H2C['vovv'],r1b,optimize=True)
    X2B += np.einsum('bf,afij->abij',Q3,t2b,optimize=True)
    Q4 = np.einsum('nmje,em->nj',H2C['ooov'],r1b,optimize=True)
    X2B -= np.einsum('nj,abin->abij',Q4,t2b,optimize=True)


    return X2B

def build_HR_2C(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys):

    t2c = cc_t['t2c']
    vC = ints['vC']
    vB = ints['vB']

    X2C = 0.0
    D1 = -np.einsum('mi,abmj->abij',H1B['oo'],r2c,optimize=True) # A(ij) 
    D2 = np.einsum('ae,ebij->abij',H1B['vv'],r2c,optimize=True) # A(ab)
    X2C += 0.5*np.einsum('mnij,abmn->abij',H2C['oooo'],r2c,optimize=True)
    X2C += 0.5*np.einsum('abef,efij->abij',H2C['vvvv'],r2c,optimize=True)
    D3 = np.einsum('amie,ebmj->abij',H2C['voov'],r2c,optimize=True) # A(ij)A(ab)
    D4 = np.einsum('maei,ebmj->abij',H2B['ovvo'],r2b,optimize=True) # A(ij)A(ab)
    D5 = -np.einsum('bmji,am->abij',H2C['vooo'],r1b,optimize=True) # A(ab)
    D6 = np.einsum('baje,ei->abij',H2C['vvov'],r1b,optimize=True) # A(ij)

    Q1 = -0.5*np.einsum('mnef,bfmn->eb',vC['oovv'],r2c,optimize=True)
    D7 = np.einsum('eb,aeij->abij',Q1,t2c,optimize=True) # A(ab)
    Q2 = -np.einsum('nmfe,fbnm->eb',vB['oovv'],r2b,optimize=True)
    D8 = np.einsum('eb,aeij->abij',Q2,t2c,optimize=True) # A(ab)

    Q1 = 0.5*np.einsum('mnef,efjn->mj',vC['oovv'],r2c,optimize=True)
    D9 = -np.einsum('mj,abim->abij',Q1,t2c,optimize=True) # A(ij)
    Q2 = np.einsum('nmfe,fenj->mj',vB['oovv'],r2b,optimize=True)
    D10 = -np.einsum('mj,abim->abij',Q2,t2c,optimize=True) # A(ij)

    Q1 = np.einsum('amfe,em->af',H2C['vovv'],r1b,optimize=True)
    D11 = np.einsum('af,fbij->abij',Q1,t2c,optimize=True) # A(ab)
    Q2 = np.einsum('nmie,em->ni',H2C['ooov'],r1b,optimize=True)
    D12 = -np.einsum('ni,abnj->abij',Q2,t2c,optimize=True) # A(ij)

    Q1 = np.einsum('maef,em->af',H2B['ovvv'],r1a,optimize=True)
    D13 = np.einsum('af,fbij->abij',Q1,t2c,optimize=True) # A(ab)
    Q2 = np.einsum('mnei,em->ni',H2B['oovo'],r1a,optimize=True)
    D14 = -np.einsum('ni,abnj->abij',Q2,t2c,optimize=True) # A(ij)

    D_ij = D1 + D6 + D9 + D10 + D12 + D14
    D_ab = D2 + D5 + D7 + D8  + D11 + D13
    D_abij = D3 + D4

    D_ij -= np.einsum('abij->abji',D_ij,optimize=True)
    D_ab -= np.einsum('abij->baij',D_ab,optimize=True)
    D_abij += -np.einsum('abij->baij',D_abij,optimize=True)\
    -np.einsum('abij->abji',D_abij,optimize=True)\
    +np.einsum('abij->baji',D_abij,optimize=True)    

    X2C += D_ij + D_ab + D_abij

    return X2C

def cis(nroot,ints,sys):

    fA = ints['fA']
    fB = ints['fB']
    vA = ints['vA']
    vB = ints['vB']
    vC = ints['vC']

    n1a = sys['Nocc_a'] * sys['Nunocc_a']
    n1b = sys['Nocc_b'] * sys['Nunocc_b']

    HAA = np.zeros((n1a,n1a))
    HAB = np.zeros((n1a,n1b))
    HBA = np.zeros((n1b,n1a))
    HBB = np.zeros((n1b,n1b))

    ct1 = 0 
    for a in range(sys['Nunocc_a']):
        for i in range(sys['Nocc_a']):
            ct2 = 0 
            for b in range(sys['Nunocc_a']):
                for j in range(sys['Nocc_a']):
                    HAA[ct1,ct2] += vA['voov'][a,j,i,b]
                    HAA[ct1,ct2] += (i == j) * fA['vv'][a,b]
                    HAA[ct1,ct2] -= (a == b) * fA['oo'][j,i]
                    ct2 += 1
            ct1 += 1
    ct1 = 0
    for a in range(sys['Nunocc_a']):
        for i in range(sys['Nocc_a']):
            ct2 = 0 
            for b in range(sys['Nunocc_b']):
                for j in range(sys['Nocc_b']):
                    HAB[ct1,ct2] += vB['voov'][a,j,i,b]
                    ct2 += 1
            ct1 += 1
    ct1 = 0
    for a in range(sys['Nunocc_b']):
        for i in range(sys['Nocc_b']):
            ct2 = 0 
            for b in range(sys['Nunocc_a']):
                for j in range(sys['Nocc_a']):
                    HBA[ct1,ct2] += vB['ovvo'][j,a,b,i]
                    ct2 += 1
            ct1 += 1
    ct1 = 0 
    for a in range(sys['Nunocc_b']):
        for i in range(sys['Nocc_b']):
            ct2 = 0 
            for b in range(sys['Nunocc_b']):
                for j in range(sys['Nocc_b']):
                    HBB[ct1,ct2] += vC['voov'][a,j,i,b]
                    HBB[ct1,ct2] += (i == j) * fB['vv'][a,b]
                    HBB[ct1,ct2] -= (a == b) * fB['oo'][j,i]
                    ct2 += 1
            ct1 += 1

    H = np.hstack( (np.vstack((HAA,HBA)), np.vstack((HAB,HBB))) )

    E_cis, C = np.linalg.eig(H) 
    idx = np.argsort(E_cis)
    E_cis = E_cis[idx[:nroot]]
    C = C[:,idx[:nroot]]

    print(E_cis)

    return C, E_cis

def calc_cc_energy(cc_t,ints):

    vA = ints['vA']
    vB = ints['vB']
    vC = ints['vC']
    fA = ints['fA']
    fB = ints['fB']
    t1a = cc_t['t1a']
    t1b = cc_t['t1b']
    t2a = cc_t['t2a']
    t2b = cc_t['t2b']
    t2c = cc_t['t2c']

    Ecorr = 0.0
    Ecorr += np.einsum('me,em->',fA['ov'],t1a,optimize=True)
    Ecorr += np.einsum('me,em->',fB['ov'],t1b,optimize=True)
    Ecorr += 0.25*np.einsum('mnef,efmn->',vA['oovv'],t2a,optimize=True)
    Ecorr += np.einsum('mnef,efmn->',vB['oovv'],t2b,optimize=True)
    Ecorr += 0.25*np.einsum('mnef,efmn->',vC['oovv'],t2c,optimize=True)
    Ecorr += 0.5*np.einsum('mnef,fn,em->',vA['oovv'],t1a,t1a,optimize=True)
    Ecorr += 0.5*np.einsum('mnef,fn,em->',vC['oovv'],t1b,t1b,optimize=True)
    Ecorr += np.einsum('mnef,em,fn->',vB['oovv'],t1a,t1b,optimize=True)

    return Ecorr

def test_updates(matfile,ints,sys):

    from scipy.io import loadmat
    from HBar_module import HBar_CCSD
    from fortran_cis import cis_hamiltonian

    print('')
    print('TEST SUBROUTINE:')
    print('Loading Matlab .mat file from {}'.format(matfile))
    print('')

    data_dict = loadmat(matfile)
    cc_t = data_dict['cc_t']

    r1a = data_dict['r1a']
    r1b = data_dict['r1b']
    r2a = data_dict['r2a']
    r2b = data_dict['r2b']
    r2c = data_dict['r2c']

    t1a = cc_t['t1a'][0,0]
    t1b = cc_t['t1b'][0,0]
    t2a = cc_t['t2a'][0,0]
    t2b = cc_t['t2b'][0,0]
    t2c = cc_t['t2c'][0,0]

    fA = ints['fA']
    fB = ints['fB']
    vA = ints['vA']
    vB = ints['vB']
    vC = ints['vC']

    cc_t = {'t1a' : t1a, 't1b' : t1b, 't2a' : t2a, 't2b' : t2b, 't2c' : t2c}
    Ecorr = calc_cc_energy(cc_t,ints)
    print('Correlation energy = {}'.format(Ecorr))

    H1A,H1B,H2A,H2B,H2C = HBar_CCSD(cc_t,ints,sys)

    # test r1a update
    X1A = build_HR_1A(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    print('|X1A| = {}'.format(np.linalg.norm(X1A)))

    # test r1b update
    X1B = build_HR_1B(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    print('|X1B| = {}'.format(np.linalg.norm(X1B)))

    # test r2a update
    X2A = build_HR_2A(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    print('|X2A| = {}'.format(np.linalg.norm(X2A)))

    # test r2b update
    X2B = build_HR_2B(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    print('|X2B| = {}'.format(np.linalg.norm(X2B)))

    # test t2c update
    X2C = build_HR_2C(r1a,r1b,r2a,r2b,r2c,cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    print('|X2C| = {}'.format(np.linalg.norm(X2C)))

    # test CIS 
    nroot = 10
    # H_cis = cis_hamiltonian.build_cis(fA['oo'],fA['vv'],fB['oo'],fB['vv'],\
    #     vA['voov'],vB['voov'],vB['ovvo'],vC['voov'])
    C, E_cis = cis(nroot,ints,sys)
    print('CIS Guess')
    for i in range(nroot):
        print('E%d = %f'% (i+1, E_cis[i]))


    return