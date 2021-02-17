import numpy as np

def parse_onebody(filename,sys):
    """This function reads the onebody.inp file from GAMESS
    and returns a numpy matrix"""

    Norb = sys['Norb']
    
    e1int = np.zeros((Norb,Norb))
    
    try:
        print('     onebody file : {}'.format(filename))
        with open(filename) as f_in:
            lines = f_in.readlines()
            ct = 0
            for i in range(Norb):
                for j in range(i+1):
                    val = float(lines[ct].split()[0])
                    e1int[i,j] = val
                    e1int[j,i] = val
                    ct += 1
    except IOError:
        print('Error: {} does not appear to exist.'.format(filename))
        sys.exit()

    return e1int
            
    
def parse_twobody(filename,sys):
    """This function reads the twobody.inp file from GAMESS
    and returns a numpy matrix"""

    try:
        print('     twobody file : {}'.format(filename))

        Norb = sys['Norb']

        # initialize numpy array
        e2int = np.zeros((Norb, Norb, Norb, Norb))

        # open file
        with open(filename) as f_in:
            
            # loop over lines
            for line in f_in:
                
                # split fields and parse
                fields = line.split()
                indices = tuple(map(int, fields[:4]))
                val = float(fields[4])
                
                # check whether value is nuclear repulsion
                # fill matrix otherwise
                if sum(indices) == 0:
                    e_nn = val
                else:
                    indices = tuple(i - 1 for i in indices)
                    e2int[indices] = val
                    
        # convert e2int from chemist notation (ia|jb) to
        # physicist notation <ij|ab>
        e2int = np.einsum('iajb->ijab', e2int)

    except IOError:
        print('Error: {} does not appear to exist.'.format(filename))
        sys.exit()

    return e_nn, e2int

def build_v(e2int):
    """This function generates the antisymmetrized version of the
    twobody matrix"""

    v_aa = e2int - np.einsum("pqrs->pqsr", e2int)
    v_ab = e2int
    v_bb = e2int - np.einsum('pqrs->pqsr',e2int)

    v = {
        "A" : v_aa,
        "B" : v_ab,
        "C" : v_bb
    }

    return v

def build_f(e1int,v,sys):
    """This function generates the Fock matrix using the formula
       F = Z + G where G is \sum_{i} <pi|v|qi>_A split for different
       spin cases"""

    Nocc_a = sys['Nocc_a']+sys['Nfrz_a']
    Nocc_b = sys['Nocc_b']+sys['Nfrz_b']


    # <p|f|q> = <p|z|q> + <pi|v|qi> + <pi~|v|qi~>
    f_a = e1int + np.einsum('piqi->pq',v['A'][:,:Nocc_a,:,:Nocc_a]) + np.einsum('piqi->pq',v['B'][:,:Nocc_b,:,:Nocc_b])

    # <p~|f|q~> = <p~|z|q~> + <p~i~|v|q~i~> + <ip~|v|iq~>
    f_b = e1int + np.einsum('piqi->pq',v['C'][:,:Nocc_b,:,:Nocc_b]) + np.einsum('ipiq->pq',v['B'][:Nocc_a,:,:Nocc_a,:])

    f = {
        "A" : f_a,
        "B" : f_b
    }

    return f

def slice_integrals(f,v,sys):

    Nocc_a = sys['Nocc_a']
    Nocc_b = sys['Nocc_b']
    Nunocc_a = sys['Nunocc_a']
    Nunocc_b = sys['Nunocc_b']

    oa = slice(sys['Nfrz_a'],sys['Nocc_a']+sys['Nfrz_a'])
    ob = slice(sys['Nfrz_b'],sys['Nocc_b']+sys['Nfrz_b'])
    ua = slice(sys['Nocc_a']+sys['Nfrz_a'],sys['Norb'])
    ub = slice(sys['Nocc_b']+sys['Nfrz_b'],sys['Norb'])

    fA_oo = f['A'][oa,oa]
    fA_ov = f['A'][oa,ua]
    fA_vo = f['A'][ua,oa]
    fA_vv = f['A'][ua,ua]

    fB_oo = f['B'][ob,ob]
    fB_ov = f['B'][ob,ub]
    fB_vo = f['B'][ub,ob]
    fB_vv = f['B'][ub,ub]

    vA_oooo = v['A'][oa,oa,oa,oa]
    vA_ooov = v['A'][oa,oa,oa,ua]
    vA_oovo = v['A'][oa,oa,ua,oa]
    vA_ovoo = v['A'][oa,ua,oa,oa]
    vA_vooo = v['A'][ua,oa,oa,oa]
    vA_oovv = v['A'][oa,oa,ua,ua]
    vA_ovov = v['A'][oa,ua,oa,ua]
    vA_voov = v['A'][ua,oa,oa,ua]
    vA_ovvo = v['A'][oa,ua,ua,oa]
    vA_vovo = v['A'][ua,oa,ua,oa]
    vA_vvoo = v['A'][ua,ua,oa,oa]
    vA_ovvv = v['A'][oa,ua,ua,ua]
    vA_vovv = v['A'][ua,oa,ua,ua]
    vA_vvov = v['A'][ua,ua,oa,ua]
    vA_vvvo = v['A'][ua,ua,ua,oa]
    vA_vvvv = v['A'][ua,ua,ua,ua]

    vB_oooo = v['B'][oa,ob,oa,ob]
    vB_ooov = v['B'][oa,ob,oa,ub]
    vB_oovo = v['B'][oa,ob,ua,ob]
    vB_ovoo = v['B'][oa,ub,oa,ob]
    vB_vooo = v['B'][ua,ob,oa,ob]
    vB_oovv = v['B'][oa,ob,ua,ub]
    vB_ovov = v['B'][oa,ub,oa,ub]
    vB_voov = v['B'][ua,ob,oa,ub]
    vB_ovvo = v['B'][oa,ub,ua,ob]
    vB_vovo = v['B'][ua,ob,ua,ob]
    vB_vvoo = v['B'][ua,ub,oa,ob]
    vB_vvvo = v['B'][ua,ub,ua,ob]
    vB_vvov = v['B'][ua,ub,oa,ub]
    vB_vovv = v['B'][ua,ob,ua,ub]
    vB_ovvv = v['B'][oa,ub,ua,ub]
    vB_vvvv = v['B'][ua,ub,ua,ub]

    vC_oooo = v['C'][ob,ob,ob,ob]
    vC_ooov = v['C'][ob,ob,ob,ub]
    vC_oovo = v['C'][ob,ob,ub,ob]
    vC_ovoo = v['C'][ob,ub,ob,ob]
    vC_vooo = v['C'][ub,ob,ob,ob]
    vC_oovv = v['C'][ob,ob,ub,ub]
    vC_ovov = v['C'][ob,ub,ob,ub]
    vC_voov = v['C'][ub,ob,ob,ub]
    vC_ovvo = v['C'][ob,ub,ub,ob]
    vC_vovo = v['C'][ub,ob,ub,ob]
    vC_vvoo = v['C'][ub,ub,ob,ob]
    vC_ovvv = v['C'][ob,ub,ub,ub]
    vC_vovv = v['C'][ub,ob,ub,ub]
    vC_vvov = v['C'][ub,ub,ob,ub]
    vC_vvvo = v['C'][ub,ub,ub,ob]
    vC_vvvv = v['C'][ub,ub,ub,ub]

    fA = {'oo' : fA_oo, 'ov' : fA_ov, 'vo' : fA_vo, 'vv' : fA_vv}
    fB = {'oo' : fB_oo, 'ov' : fB_ov, 'vo' : fB_vo, 'vv' : fB_vv}
    vA = {'oooo' : vA_oooo, 'ooov' : vA_ooov, 'oovo' : vA_oovo, 'ovoo' : vA_ovoo, \
          'vooo' : vA_vooo, 'oovv' : vA_oovv, 'ovov' : vA_ovov, 'ovvo' : vA_ovvo, \
          'vovo' : vA_vovo, 'voov' : vA_voov, 'vvoo' : vA_vvoo, 'vvvo' : vA_vvvo, \
          'vvov' : vA_vvov, 'vovv' : vA_vovv, 'ovvv' : vA_ovvv, 'vvvv' : vA_vvvv}
    vB = {'oooo' : vB_oooo, 'ooov' : vB_ooov, 'oovo' : vB_oovo, 'ovoo' : vB_ovoo, \
          'vooo' : vB_vooo, 'oovv' : vB_oovv, 'ovov' : vB_ovov, 'ovvo' : vB_ovvo, \
          'vovo' : vB_vovo, 'voov' : vB_voov, 'vvoo' : vB_vvoo, 'vvvo' : vB_vvvo, \
          'vvov' : vB_vvov, 'vovv' : vB_vovv, 'ovvv' : vB_ovvv, 'vvvv' : vB_vvvv}
    vC = {'oooo' : vC_oooo, 'ooov' : vC_ooov, 'oovo' : vC_oovo, 'ovoo' : vC_ovoo, \
          'vooo' : vC_vooo, 'oovv' : vC_oovv, 'ovov' : vC_ovov, 'ovvo' : vC_ovvo, \
          'vovo' : vC_vovo, 'voov' : vC_voov, 'vvoo' : vC_vvoo, 'vvvo' : vC_vvvo, \
          'vvov' : vC_vvov, 'vovv' : vC_vovv, 'ovvv' : vC_ovvv, 'vvvv' : vC_vvvv}

    return fA, fB, vA, vB, vC

def get_integrals(onebody_file,twobody_file,sys):

    print('')
    print('  Reading integrals...')

    e1int = parse_onebody(onebody_file,sys)
    e_nn, e2int = parse_twobody(twobody_file,sys)

    print('  Integrals read successfully!')

    Escf = e_nn
    for i in range(sys['Nocc_a']+sys['Nfrz_a']):
        Escf += e1int[i,i]
    for i in range(sys['Nocc_b']+sys['Nfrz_b']):
        Escf += e1int[i,i]
    for i in range(sys['Nocc_a']+sys['Nfrz_a']):
        for j in range(sys['Nocc_a']+sys['Nfrz_a']):
           Escf += 0.5*(e2int[i,j,i,j] - e2int[i,j,j,i])
    for i in range(sys['Nocc_a']+sys['Nfrz_a']):
        for j in range(sys['Nocc_b']+sys['Nfrz_b']):
           Escf += e2int[i,j,i,j]
    for i in range(sys['Nocc_b']+sys['Nfrz_b']):
        for j in range(sys['Nocc_b']+sys['Nfrz_b']):
           Escf += 0.5*(e2int[i,j,i,j] - e2int[i,j,j,i])

    v = build_v(e2int)
    f = build_f(e1int,v,sys)
    fA,fB,vA,vB,vC = slice_integrals(f,v,sys)

    # Escf = e_nn
    # Escf += np.einsum('ii->',fA['oo'],optimize=True)
    # Escf += np.einsum('ii->',fB['oo'],optimize=True)
    # Escf -= 0.5*np.einsum('ijij->',vA['oooo'],optimize=True)
    # Escf -= 0.5*np.einsum('ijij->',vC['oooo'],optimize=True)
    # Escf -= np.einsum('ijij->',vB['oooo'],optimize=True)


    ints = {'fA' : fA, 'fB' : fB, 'vA' : vA, 'vB' : vB, 'vC' : vC, 'Vnuc' : e_nn, 'Escf' : Escf}

    return ints
