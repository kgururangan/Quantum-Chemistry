import numpy as np 
import time 

def crcc24(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys):

    print('\n==================================++Entering CR-CC(2,4) Routine++=============================')
    
    t_start = time.time()
    
    # get fock matrices
    fA = ints['fA']; fB = ints['fB']
    
    # get the 3-body HBar triples diagonal
    D3A, D3B, D3C, D3D = triples_3body_diagonal(cc_t,ints,sys)

    # MM correction containers
    deltaA = 0.0; # using MP denominator -(f_aa - f_ii + f_bb - f_jj + f_cc - f_kk)
    deltaB = 0.0; # using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) | phi_{ijk}^{abc}>
    deltaC = 0.0; # using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) | phi_{ijk}^{abc}>
    deltaD = 0.0; # using EN denominator -<phi_{ijk}^{abc} | H_1(CCSD) + H_2(CCSD) + H_3(CCSD) | phi_{ijk}^{abc}>
        
    # MM24A correction
    MM24A = build_MM24A(cc_t,H2A,sys)
    L4A = build_L4A(cc_t,ints,sys)
    for a in range(sys['Nunocc_a']):
        for b in range(a+1,sys['Nunocc_a']):
            for c in range(b+1,sys['Nunocc_a']):
                for d in range(c+1,sys['Nunocc_a']):
                    for i in range(sys['Nocc_a']):
                        for j in range(i+1,sys['Nocc_a']):
                            for k in range(j+1,sys['Nocc_a']):
                                for l in range(k+1,sys['Nocc_a']):

                                    LM = L4A[a,b,c,d,i,j,k,l] * MM24A[a,b,c,d,i,j,k,l]
                                    
                                    DMP = fA['vv'][a,a] + fA['vv'][b,b] + fA['vv'][c,c] + fA['vv'][d,d]\
                                    -fA['oo'][i,i] - fA['oo'][j,j] - fA['oo'][k,k] - fA['oo'][l,l]
                    
                                    D_A = -1.0 * DMP

                                    deltaA += LM/D_A

    
    # MM24B correction
    MM24B = build_MM24B(cc_t,H2B,H2A,sys)
    L4B = build_L4B(cc_t,ints,sys)
    for a in range(sys['Nunocc_a']):
        for b in range(a+1,sys['Nunocc_a']):
            for c in range(b+1,sys['Nunocc_a']):
                for d in range(sys['Nunocc_b']):
                    for i in range(sys['Nocc_a']):
                        for j in range(i+1,sys['Nocc_a']):
                            for k in range(j+1,sys['Nocc_a']):
                                for l in range(sys['Nocc_b']):

                                    LM = L4B[a,b,c,d,i,j,k,l] * MM24B[a,b,c,d,i,j,k,l]
                                    
                                    DMP = fA['vv'][a,a] + fA['vv'][b,b] + fA['vv'][c,c] + fB['vv'][d,d]\
                                    -fA['oo'][i,i] - fA['oo'][j,j] - fA['oo'][k,k] - fB['oo'][l,l]
                    
                                    D_A = -1.0 * DMP

                                    deltaA += LM/D_A
    
    # MM24C correction
    MM24C = build_MM24C(cc_t,H2B,H2A,H2C,sys)
    L4C = build_L4C(cc_t,ints,sys)
    for a in range(sys['Nunocc_a']):
        for b in range(a+1,sys['Nunocc_a']):
            for c in range(sys['Nunocc_b']):
                for d in range(c+1,sys['Nunocc_b']):
                    for i in range(sys['Nocc_a']):
                        for j in range(i+1,sys['Nocc_a']):
                            for k in range(sys['Nocc_b']):
                                for l in range(k+1,sys['Nocc_b']):

                                    LM = L4C[a,b,c,d,i,j,k,l] * MM24C[a,b,c,d,i,j,k,l]
                                    
                                    DMP = fA['vv'][a,a] + fA['vv'][b,b] + fB['vv'][c,c] + fB['vv'][d,d]\
                                    -fA['oo'][i,i] - fA['oo'][j,j] - fB['oo'][k,k] - fB['oo'][l,l]
                    
                                    D_A = -1.0 * DMP

                                    deltaA += LM/D_A
    
    # MM24D correction
    MM24D = build_MM24D(cc_t,H2B,H2C,sys)
    L4D = build_L4D(cc_t,ints,sys)
    for a in range(sys['Nunocc_a']):
        for b in range(sys['Nunocc_b']):
            for c in range(b+1,sys['Nunocc_b']):
                for d in range(c+1,sys['Nunocc_b']):
                    for i in range(sys['Nocc_a']):
                        for j in range(sys['Nocc_b']):
                            for k in range(j+1,sys['Nocc_b']):
                                for l in range(k+1,sys['Nocc_b']):

                                    LM = L4D[a,b,c,d,i,j,k,l] * MM24D[a,b,c,d,i,j,k,l]
                                    
                                    DMP = fA['vv'][a,a] + fB['vv'][b,b] + fB['vv'][c,c] + fB['vv'][d,d]\
                                    -fA['oo'][i,i] - fB['oo'][j,j] - fB['oo'][k,k] - fB['oo'][l,l]
                    
                                    D_A = -1.0 * DMP

                                    deltaA += LM/D_A

    # MM24E correction
    MM24E = build_MM24E(cc_t,H2C,sys)
    L4E = build_L4E(cc_t,ints,sys)
    for a in range(sys['Nunocc_b']):
        for b in range(a+1,sys['Nunocc_b']):
            for c in range(b+1,sys['Nunocc_b']):
                for d in range(c+1,sys['Nunocc_b']):
                    for i in range(sys['Nocc_b']):
                        for j in range(i+1,sys['Nocc_b']):
                            for k in range(j+1,sys['Nocc_b']):
                                for l in range(k+1,sys['Nocc_b']):

                                    LM = L4E[a,b,c,d,i,j,k,l] * MM24E[a,b,c,d,i,j,k,l]
                                    
                                    DMP = fB['vv'][a,a] + fB['vv'][b,b] + fB['vv'][c,c] + fB['vv'][d,d]\
                                    -fB['oo'][i,i] - fB['oo'][j,j] - fB['oo'][k,k] - fB['oo'][l,l]
                    
                                    D_A = -1.0 * DMP

                                    deltaA += LM/D_A
    
    Ecorr = calc_cc_energy(cc_t,ints)

    EcorrA = Ecorr + deltaA; EcorrB = Ecorr + deltaB; EcorrC = Ecorr + deltaC; EcorrD = Ecorr + deltaD

    E23A = ints['Escf'] + EcorrA
 
    print('CR-CC(2,3)_A = {} Eh     Ecorr_A = {} Eh     Delta_A = {} Eh'.format(E23A,EcorrA,deltaA))

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return 

def build_MM24A(cc_t,H2A,sys):

    print('MMCC(2,4)A construction... ')

    MM24A = np.zeros((sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_a'],\
    sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_a']))
    

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return MM24A

def build_MM24B(cc_t,H2A,H2B,sys):

    print('MMCC(2,4)B construction... ')

    MM24B = np.zeros((sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_b'],\
    sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_b']))
    

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return MM24B

def build_MM24C(cc_t,H2A,H2B,H2C,sys):

    print('MMCC(2,4)C construction... ')

    MM24C = np.zeros((sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_b'],sys['Nunocc_b'],\
    sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_b'],sys['Nocc_b']))
    

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return MM24C

def build_MM24D(cc_t,H2B,H2C,sys):

    print('MMCC(2,4)D construction... ')

    MM24D = np.zeros((sys['Nunocc_a'],sys['Nunocc_b'],sys['Nunocc_b'],sys['Nunocc_b'],\
    sys['Nocc_a'],sys['Nocc_b'],sys['Nocc_b'],sys['Nocc_b']))
    

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return MM24D


def build_MM24E(cc_t,H2C,sys):

    print('MMCC(2,4)E construction... ')

    MM24E = np.zeros((sys['Nunocc_b'],sys['Nunocc_b'],sys['Nunocc_b'],sys['Nunocc_b'],\
    sys['Nocc_b'],sys['Nocc_b'],sys['Nocc_b'],sys['Nocc_b']))
    

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return MM24E

def build_L4A(cc_t,ints,sys):

    print('Approximate L4A construction... ')
    
    t_start = time.time()

    vA = ints['vA']
    l2a = cc_t['l2a']

    L4A = np.zeros((sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_a'],\
    sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_a']))

    L4A += np.einsum('ijab,cdkl->abcdijkl',vA['oovv'],l2a,optimize=True)
    L4A = L4A - permute(L4A,[1,2,3,4,7,6,5,8]) - permute(L4A,[1,2,3,4,8,6,7,5])\
        -permute(L4A,[1,2,3,4,5,7,6,8]) - permute(L4A,[1,2,3,4,5,8,7,6])\
        +permute(L4A,[1,2,3,4,7,8,5,6])
    L4A = L4A - permute(L4A,[3,2,1,4,5,6,7,8]) - permute(L4A,[4,2,3,1,5,6,7,8])\
        -permute(L4A,[1,3,2,4,5,6,7,8]) - permute(L4A,[1,4,3,2,5,6,7,8])\
        +permute(L4A,[3,4,1,2,5,6,7,8])

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return L4A

def build_L4B(cc_t,ints,sys):

    print('Approximate L4B construction... ')
    
    t_start = time.time()

    vA = ints['vA']
    vB = ints['vB']
    l2a = cc_t['l2a']
    l2b = cc_t['l2b']

    L4B = np.zeros((sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_b'],\
    sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_b']))

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return L4B

def build_L4C(cc_t,ints,sys):

    print('Approximate L4C construction... ')
    
    t_start = time.time()

    vA = ints['vA']
    vB = ints['vB']
    vC = ints['vC']
    l2a = cc_t['l2a']
    l2b = cc_t['l2b']
    l2c = cc_t['l2c']

    L4C = np.zeros((sys['Nunocc_a'],sys['Nunocc_a'],sys['Nunocc_b'],sys['Nunocc_b'],\
    sys['Nocc_a'],sys['Nocc_a'],sys['Nocc_b'],sys['Nocc_b']))

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return L4C


def build_L4D(cc_t,ints,sys):

    print('Approximate L4D construction... ')
    
    t_start = time.time()

    vB = ints['vB']
    vC = ints['vC']
    l2b = cc_t['l2b']
    l2c = cc_t['l2c']

    L4D = np.zeros((sys['Nunocc_a'],sys['Nunocc_b'],sys['Nunocc_b'],sys['Nunocc_b'],\
    sys['Nocc_a'],sys['Nocc_b'],sys['Nocc_b'],sys['Nocc_b']))

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return L4D

def build_L4E(cc_t,ints,sys):

    print('Approximate L4E construction... ')
    
    t_start = time.time()

    vC = ints['vC']
    l2c = cc_t['l2c']

    L4E = np.zeros((sys['Nunocc_b'],sys['Nunocc_b'],sys['Nunocc_b'],sys['Nunocc_b'],\
    sys['Nocc_b'],sys['Nocc_b'],sys['Nocc_b'],sys['Nocc_b']))

    L4E += np.einsum('ijab,cdkl->abcdijkl',vC['oovv'],l2c,optimize=True)
    L4E = L4E - permute(L4E,[1,2,3,4,7,6,5,8]) - permute(L4E,[1,2,3,4,8,6,7,5])\
        -permute(L4E,[1,2,3,4,5,7,6,8]) - permute(L4E,[1,2,3,4,5,8,7,6])\
        +permute(L4E,[1,2,3,4,7,8,5,6])
    L4E = L4E - permute(L4E,[3,2,1,4,5,6,7,8]) - permute(L4E,[4,2,3,1,5,6,7,8])\
        -permute(L4E,[1,3,2,4,5,6,7,8]) - permute(L4E,[1,4,3,2,5,6,7,8])\
        +permute(L4E,[3,4,1,2,5,6,7,8])

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return L4E

def triples_3body_diagonal(cc_t,ints,sys):

    print('\nCalculating 3-body triples diagonal... ')
    t_start = time.time()

    vA = ints['vA']
    vB = ints['vB']
    vC = ints['vC']
    t2a = cc_t['t2a']
    t2b = cc_t['t2b']
    t2c = cc_t['t2c']

    d3A_V = lambda a,i,b: -np.dot(vA['oovv'][i,:,a,b].T,t2a[a,b,i,:])
    d3A_O = lambda a,i,j:  np.dot(vA['oovv'][i,j,a,:].T,t2a[a,:,i,j])
    
    d3B_V = lambda a,i,c: -np.dot(vB['oovv'][i,:,a,c].T,t2b[a,c,i,:])
    d3B_O = lambda a,i,k:  np.dot(vB['oovv'][i,k,a,:].T,t2b[a,:,i,k])
    
    d3C_V = lambda a,k,c: -np.dot(vB['oovv'][:,k,a,c].T,t2b[a,c,:,k])
    d3C_O = lambda c,i,k:  np.dot(vB['oovv'][i,k,:,c].T,t2b[:,c,i,k])
    
    d3D_V = lambda a,i,b: -np.dot(vC['oovv'][i,:,a,b].T,t2c[a,b,i,:])
    d3D_O = lambda a,i,j:  np.dot(vC['oovv'][i,j,a,:].T,t2c[a,:,i,j])
    
    D3A_V = np.zeros((sys['Nunocc_a'],sys['Nocc_a'],sys['Nunocc_a']))
    D3A_O = np.zeros((sys['Nunocc_a'],sys['Nocc_a'],sys['Nocc_a']))
    D3B_V = np.zeros((sys['Nunocc_a'],sys['Nocc_a'],sys['Nunocc_b']))
    D3B_O = np.zeros((sys['Nunocc_a'],sys['Nocc_a'],sys['Nocc_b']))
    D3C_V = np.zeros((sys['Nunocc_a'],sys['Nocc_b'],sys['Nunocc_b']))
    D3C_O = np.zeros((sys['Nunocc_b'],sys['Nocc_a'],sys['Nocc_b']))
    D3D_V = np.zeros((sys['Nunocc_b'],sys['Nocc_b'],sys['Nunocc_b']))
    D3D_O = np.zeros((sys['Nunocc_b'],sys['Nocc_b'],sys['Nocc_b']))

    # A diagonal
    for a in range(sys['Nunocc_a']):
        for i in range(sys['Nocc_a']):
            for b in range(sys['Nunocc_a']):
                D3A_V[a,i,b] = d3A_V(a,i,b)
    for a in range(sys['Nunocc_a']):
        for i in range(sys['Nocc_a']):
            for j in range(sys['Nocc_a']):
                D3A_O[a,i,j] = d3A_O(a,i,j)
    
    # B diagonal
    for a in range(sys['Nunocc_a']):
        for i in range(sys['Nocc_a']):
            for c in range(sys['Nunocc_b']):
                D3B_V[a,i,c] = d3B_V(a,i,c)
    for a in range(sys['Nunocc_a']):
        for i in range(sys['Nocc_a']):
            for k in range(sys['Nocc_b']):
                D3B_O[a,i,k] = d3B_O(a,i,k)
    
   # C diagonal 
    for a in range(sys['Nunocc_a']):
        for k in range(sys['Nocc_b']):
            for c in range(sys['Nunocc_b']):
                D3C_V[a,k,c] = d3C_V(a,k,c)
    for c in range(sys['Nunocc_b']):
        for i in range(sys['Nocc_a']):
            for k in range(sys['Nocc_b']):
                D3C_O[c,i,k] = d3C_O(c,i,k)
    
    # D diagonal 
    for a in range(sys['Nunocc_b']):
        for i in range(sys['Nocc_b']):
            for b in range(sys['Nunocc_b']):
                D3D_V[a,i,b] = d3D_V(a,i,b)
    for a in range(sys['Nunocc_b']):
        for i in range(sys['Nocc_b']):
            for j in range(sys['Nocc_b']):
                D3D_O[a,i,j] = d3D_O(a,i,j)

    D3A = {'O' : D3A_O, 'V' : D3A_V}
    D3B = {'O' : D3B_O, 'V' : D3B_V}
    D3C = {'O' : D3C_O, 'V' : D3C_V}
    D3D = {'O' : D3D_O, 'V' : D3D_V}
    
    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))
    
    return D3A, D3B, D3C, D3D

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

def permute(x,perm_list):
    return x.transpose(tuple([y-1 for y in perm_list]))


def test_updates(matfile,ints,sys):

    from scipy.io import loadmat
    from HBar_module import HBar_CCSD

    data_dict = loadmat(matfile)
    cc_t = data_dict['cc_t']

    t1a = cc_t['t1a'][0,0]
    t1b = cc_t['t1b'][0,0]
    t2a = cc_t['t2a'][0,0]
    t2b = cc_t['t2b'][0,0]
    t2c = cc_t['t2c'][0,0]

    l1a = data_dict['l1a']
    l1b = data_dict['l1b']
    l2a = data_dict['l2a']
    l2b = data_dict['l2b']
    l2c = data_dict['l2c']

    cc_t = {'t1a' : t1a, 't1b' : t1b, 't2a' : t2a, 't2b' : t2b, 't2c' : t2c,
            'l1a' : l1a, 'l1b' : l1b, 'l2a' : l2a, 'l2b' : l2b, 'l2c' : l2c}

    H1A,H1B,H2A,H2B,H2C = HBar_CCSD(cc_t,ints,sys)

    # # test MM24A update
    # MM24A = build_MM24A(cc_t,H2A,sys)
    # print('|MM24A| = {}'.format(np.linalg.norm(MM24A)))

    # # test MM24B update
    # MM24B = build_MM24B(cc_t,H2A,H2B,sys)
    # print('|MM24B| = {}'.format(np.linalg.norm(MM24B)))

    # # test MM24C update
    # MM24C = build_MM24C(cc_t,H2A,H2B,H2C,sys)
    # print('|MM24C| = {}'.format(np.linalg.norm(MM24C)))

    # # test MM24D update
    # MM24D = build_MM24D(cc_t,H2B,H2C,sys)
    # print('|MM24D| = {}'.format(np.linalg.norm(MM24D)))

    # # test MM24D update
    # MM24E = build_MM24E(cc_t,H2C,sys)
    # print('|MM24E| = {}'.format(np.linalg.norm(MM24E)))

    # test L4A update
    L4A = build_L4A(cc_t,ints,sys)
    print('|L4A| = {}'.format(np.linalg.norm(L4A)))

    # # test L4B update
    # L4B = build_L4B(cc_t,ints,sys)
    # print('|L4B| = {}'.format(np.linalg.norm(L4B)))

    # # test L4C update
    # L4C = build_L4C(cc_t,ints,sys)
    # print('|L4C| = {}'.format(np.linalg.norm(L4C)))

    # # test L4D update
    # L4D = build_L4D(cc_t,ints,sys)
    # print('|L4D| = {}'.format(np.linalg.norm(L4D)))

    # test L4E update
    L4E = build_L4A(cc_t,ints,sys)
    print('|L4E| = {}'.format(np.linalg.norm(L4E)))

    # test 3-body diagonal
    D3A,D3B,D3C,D3D = triples_3body_diagonal(cc_t,ints,sys)
    print('|D3A_O| = {}'.format(np.linalg.norm(D3A['O'])))
    print('|D3A_V| = {}'.format(np.linalg.norm(D3A['V'])))
    print('|D3B_O| = {}'.format(np.linalg.norm(D3B['O'])))
    print('|D3B_V| = {}'.format(np.linalg.norm(D3B['V'])))
    print('|D3C_O| = {}'.format(np.linalg.norm(D3C['O'])))
    print('|D3C_V| = {}'.format(np.linalg.norm(D3C['V'])))
    print('|D3D_O| = {}'.format(np.linalg.norm(D3D['O'])))
    print('|D3D_V| = {}'.format(np.linalg.norm(D3D['V'])))

    return