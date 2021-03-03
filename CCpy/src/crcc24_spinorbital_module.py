import numpy as np 
import time 

def crcc24(cc_t,H1,H2,ints,sys):

    print('\n==================================++Entering CR-CC(2,4) Routine++=============================')
    print('****SPINORBITAL ROUTINE****')


    # Get Fock matrix elements
    f = ints['f']
    v = ints['v']

    # Get T2
    t2 = cc_t['t2']

    # Get optimum paths for diagram contraction functions
    #path_voov_1 = np.einsum_path

    # Spin case slicing vectors (for test spin-integrated module)
    oa = slice(0,sys['Nocc']-1,2)
    ob = slice(1,sys['Nocc'],2)
    ua = slice(0,sys['Nunocc']-1,2)
    ub = slice(1,sys['Nunocc'],2)

    # Get triples diagonal
    D3 = triples_diagonal(v['oovv'],cc_t['t2'],sys)

    # MM correction containers
    deltaA = 0.0
    deltaB = 0.0
    deltaC = 0.0
    deltaD = 0.0

    # Calculate <ijklabcd|H(2)|Phi>
    #MM24, D_voov, D_oooo, D_vvvv = build_MM24(cc_t,H2,sys)
    MM24 = np.zeros((sys['Nunocc'],sys['Nunocc'],sys['Nunocc'],sys['Nunocc'],\
                    sys['Nocc'],sys['Nocc'],sys['Nocc'],sys['Nocc']))

    # Calculate <Phi|(1+L1+L2)H(2)|ijklabcd>
    L4 = build_L4(cc_t,ints,sys)
    #L4 = np.zeros((sys['Nunocc'],sys['Nunocc'],sys['Nunocc'],sys['Nunocc'],\
    #                sys['Nocc'],sys['Nocc'],sys['Nocc'],sys['Nocc']))

    # debugging D denominator for each spin case
    denom_sum_AAAA_O = 0.0
    denom_sum_AAAA_V = 0.0
    denom_sum_AAAB_O = 0.0
    denom_sum_AAAB_V = 0.0
    denom_sum_AABB_O = 0.0
    denom_sum_AABB_V = 0.0

    # Full contraction loop
    # delta(2,4) = L4(abcdijkl)*MM24(abcdijkl) / D(abcdijkl)
    for a in range(sys['Nunocc']):
        for b in range(a+1,sys['Nunocc']):
            for c in range(b+1,sys['Nunocc']):
                for d in range(c+1,sys['Nunocc']):
                    for i in range(sys['Nocc']):
                        for j in range(i+1,sys['Nocc']):
                            for k in range(j+1,sys['Nocc']):
                                for l in range(k+1,sys['Nocc']):

                                    nalpha_vir = sum([np.mod(x+1,2) for x in list((a,b,c,d))])
                                    nalpha_occ = sum([np.mod(x+1,2) for x in list((i,j,k,l))])

                                    MM24_val = build_antisym_diagram(a,b,c,d,i,j,k,l,val_1,H2['voov'],t2)\
                                                    +build_antisym_diagram(a,b,c,d,i,j,k,l,val_2,H2['oooo'],t2)\
                                                    +build_antisym_diagram(a,b,c,d,i,j,k,l,val_3,H2['vvvv'],t2)

                                    LM = L4[a,b,c,d,i,j,k,l] * MM24_val
                                    #LM = L4[a,b,c,d,i,j,k,l] * MM24[a,b,c,d,i,j,k,l]

                                    denom_mp = f['oo'][i,i] + f['oo'][j,j] + f['oo'][k,k] + f['oo'][l,l]\
                                                    -f['vv'][a,a] - f['vv'][b,b] - f['vv'][c,c] - f['vv'][d,d]

                                    denom_1 = H1['oo'][i,i] + H1['oo'][j,j] + H1['oo'][k,k] + H1['oo'][l,l]\
                                                    -H1['vv'][a,a] - H1['vv'][b,b] - H1['vv'][c,c] - H1['vv'][d,d]
            
                                    denom_2 = -H2['oooo'][j,i,j,i]-H2['oooo'][k,i,k,i]-H2['oooo'][l,i,l,i]-H2['oooo'][k,j,k,j]\
                                                    -H2['oooo'][l,j,l,j]-H2['oooo'][l,k,l,k]-H2['voov'][a,i,i,a]-H2['voov'][a,j,j,a]\
                                                    -H2['voov'][a,k,k,a]-H2['voov'][a,l,l,a]-H2['voov'][b,i,i,b]-H2['voov'][b,j,j,b]\
                                                    -H2['voov'][b,k,k,b]-H2['voov'][b,l,l,b]-H2['voov'][c,i,i,c]-H2['voov'][c,j,j,c]\
                                                    -H2['voov'][c,k,k,c]-H2['voov'][c,l,l,c]-H2['voov'][d,i,i,d]-H2['voov'][d,j,j,d]\
                                                    -H2['voov'][d,k,k,d]-H2['voov'][d,l,l,d]-H2['vvvv'][a,b,a,b]-H2['vvvv'][a,c,a,c]\
                                                    -H2['vvvv'][a,d,a,d]-H2['vvvv'][b,c,b,c]-H2['vvvv'][b,d,b,d]-H2['vvvv'][c,d,c,d]

                                    denom_3 = +D3['O'][a,i,j]+D3['O'][a,i,k]+D3['O'][a,i,l]+D3['O'][a,j,k]\
                                                    +D3['O'][a,j,l]+D3['O'][a,k,l]+D3['O'][b,i,j]+D3['O'][b,i,k]\
                                                    +D3['O'][b,i,l]+D3['O'][b,j,k]+D3['O'][b,j,l]+D3['O'][b,k,l]\
                                                    +D3['O'][c,i,j]+D3['O'][c,i,k]+D3['O'][c,i,l]+D3['O'][c,j,k]\
                                                    +D3['O'][c,j,l]+D3['O'][c,k,l]+D3['O'][d,i,j]+D3['O'][d,i,k]\
                                                    +D3['O'][d,i,l]+D3['O'][d,j,k]+D3['O'][d,j,l]+D3['O'][d,k,l]\
                                                    -D3['V'][a,i,b]-D3['V'][a,j,b]-D3['V'][a,k,b]-D3['V'][a,l,b]\
                                                    -D3['V'][a,i,c]-D3['V'][a,j,c]-D3['V'][a,k,c]-D3['V'][a,l,c]\
                                                    -D3['V'][a,i,d]-D3['V'][a,j,d]-D3['V'][a,k,d]-D3['V'][a,l,d]\
                                                    -D3['V'][b,i,c]-D3['V'][b,j,c]-D3['V'][b,k,c]-D3['V'][b,l,c]\
                                                    -D3['V'][b,i,d]-D3['V'][b,j,d]-D3['V'][b,k,d]-D3['V'][b,l,d]\
                                                    -D3['V'][c,i,d]-D3['V'][c,j,d]-D3['V'][c,k,d]-D3['V'][c,l,d]

                                    denom_3_O = +D3['O'][a,i,j]+D3['O'][a,i,k]+D3['O'][a,i,l]+D3['O'][a,j,k]\
                                                    +D3['O'][a,j,l]+D3['O'][a,k,l]+D3['O'][b,i,j]+D3['O'][b,i,k]\
                                                    +D3['O'][b,i,l]+D3['O'][b,j,k]+D3['O'][b,j,l]+D3['O'][b,k,l]\
                                                    +D3['O'][c,i,j]+D3['O'][c,i,k]+D3['O'][c,i,l]+D3['O'][c,j,k]\
                                                    +D3['O'][c,j,l]+D3['O'][c,k,l]+D3['O'][d,i,j]+D3['O'][d,i,k]\
                                                    +D3['O'][d,i,l]+D3['O'][d,j,k]+D3['O'][d,j,l]+D3['O'][d,k,l]

                                    denom_3_V = -D3['V'][a,i,b]-D3['V'][a,j,b]-D3['V'][a,k,b]-D3['V'][a,l,b]\
                                                    -D3['V'][a,i,c]-D3['V'][a,j,c]-D3['V'][a,k,c]-D3['V'][a,l,c]\
                                                    -D3['V'][a,i,d]-D3['V'][a,j,d]-D3['V'][a,k,d]-D3['V'][a,l,d]\
                                                    -D3['V'][b,i,c]-D3['V'][b,j,c]-D3['V'][b,k,c]-D3['V'][b,l,c]\
                                                    -D3['V'][b,i,d]-D3['V'][b,j,d]-D3['V'][b,k,d]-D3['V'][b,l,d]\
                                                    -D3['V'][c,i,d]-D3['V'][c,j,d]-D3['V'][c,k,d]-D3['V'][c,l,d]

                                    deltaA += LM/denom_mp
                                    deltaB += LM/denom_1
                                    deltaC += LM/(denom_1+denom_2)
                                    deltaD += LM/(denom_1+denom_2+denom_3)

                                    if nalpha_occ == nalpha_vir:

                                        if nalpha_occ == 4:
                                            denom_sum_AAAA_O += (denom_3_O)
                                            denom_sum_AAAA_V += (denom_3_V)
                                        if nalpha_occ == 3: 
                                            denom_sum_AAAB_O += (denom_3_O)
                                            denom_sum_AAAB_V += (denom_3_V)
                                        if nalpha_occ == 2:
                                            denom_sum_AABB_O += (denom_3_O)
                                            denom_sum_AABB_V += (denom_3_V)


    print('D denominator sums:')
    print('Spin case A : D3_O = {}    D3_V = {}'.format(denom_sum_AAAA_O,denom_sum_AAAA_V))
    print('Spin case B : D3_O = {}    D3_V = {}'.format(denom_sum_AAAB_O,denom_sum_AAAB_V))
    print('Spin case C : D3_O = {}    D3_V = {}'.format(denom_sum_AABB_O,denom_sum_AABB_V))


    Ecorr = cc_energy(cc_t,ints)
    EcorrA = Ecorr + deltaA; EcorrB = Ecorr + deltaB; EcorrC = Ecorr + deltaC; EcorrD = Ecorr + deltaD
    E24A = ints['Escf'] + EcorrA
    E24B = ints['Escf'] + EcorrB
    E24C = ints['Escf'] + EcorrC
    E24D = ints['Escf'] + EcorrD
 
    print('CR-CC(2,4)_A = {} Eh     Ecorr_A = {} Eh     Delta_A = {} Eh'.format(E24A,EcorrA,deltaA))
    print('CR-CC(2,4)_B = {} Eh     Ecorr_B = {} Eh     Delta_B = {} Eh'.format(E24B,EcorrB,deltaB))
    print('CR-CC(2,4)_C = {} Eh     Ecorr_C = {} Eh     Delta_C = {} Eh'.format(E24C,EcorrC,deltaC))
    print('CR-CC(2,4)_D = {} Eh     Ecorr_D = {} Eh     Delta_D = {} Eh'.format(E24D,EcorrD,deltaD))

    Ecrcc24 = {'A' : E24A, 'B' : E24B, 'C' : E24C, 'D' : E24D}
    delta24 = {'A' : deltaA, 'B' : deltaB, 'C' : deltaC, 'D' : deltaD}

    print('')
    print('SPIN-INTEGRAGTED DEBUG')
    print('-----------------------------------------')
    MM24A = MM24[ua,ua,ua,ua,oa,oa,oa,oa]
    print('|MM(2,4)A| = {}'.format(np.linalg.norm(MM24A)))

    MM24B = MM24[ua,ua,ua,ub,oa,oa,oa,ob]
    DB_voov = D_voov[ua,ua,ua,ub,oa,oa,oa,ob]
    DB_oooo = D_oooo[ua,ua,ua,ub,oa,oa,oa,ob]
    DB_vvvv = D_vvvv[ua,ua,ua,ub,oa,oa,oa,ob]
    print('|MM(2,4)B| = {}'.format(np.linalg.norm(MM24B)))
    print('|voov| = {}'.format(np.linalg.norm(DB_voov)))
    print('|oooo| = {}'.format(np.linalg.norm(DB_oooo)))
    print('|vvvv| = {}'.format(np.linalg.norm(DB_vvvv)))

    MM24C = MM24[ua,ua,ub,ub,oa,oa,ob,ob]
    DC_voov = D_voov[ua,ua,ub,ub,oa,oa,ob,ob]
    DC_oooo = D_oooo[ua,ua,ub,ub,oa,oa,ob,ob]
    DC_vvvv = D_vvvv[ua,ua,ub,ub,oa,oa,ob,ob]
    print('|MM(2,4)C| = {}'.format(np.linalg.norm(MM24C)))
    print('|voov| = {}'.format(np.linalg.norm(DC_voov)))
    print('|oooo| = {}'.format(np.linalg.norm(DC_oooo)))
    print('|vvvv| = {}'.format(np.linalg.norm(DC_vvvv)))

    MM24D = MM24[ua,ub,ub,ub,oa,ob,ob,ob]
    print('|MM(2,4)D| = {}'.format(np.linalg.norm(MM24D)))

    MM24E = MM24[ub,ub,ub,ub,ob,ob,ob,ob]
    print('|MM(2,4)E| = {}'.format(np.linalg.norm(MM24E)))

    L4A = L4[ua,ua,ua,ua,oa,oa,oa,oa]
    print('|L4A| = {}'.format(np.linalg.norm(L4A)))

    L4B = L4[ua,ua,ua,ub,oa,oa,oa,ob]
    print('|L4B| = {}'.format(np.linalg.norm(L4B)))

    L4C = L4[ua,ua,ub,ub,oa,oa,ob,ob]
    print('|L4C| = {}'.format(np.linalg.norm(L4C)))

    L4D = L4[ua,ub,ub,ub,oa,ob,ob,ob]
    print('|L4D| = {}'.format(np.linalg.norm(L4D)))

    L4E = L4[ub,ub,ub,ub,ob,ob,ob,ob]
    print('|L4E| = {}'.format(np.linalg.norm(L4E)))

    print('-----------------------------------------')

    return Ecrcc24, delta24

def val_1(a,b,c,d,i,j,k,l,h_voov,t2):
    x = -np.einsum('me,m,e->',h_voov[c,:,k,:],t2[a,d,:,l],t2[:,b,i,j],optimize=True)
    x -= np.einsum('me,m,e->',h_voov[c,:,l,:],t2[a,d,:,k],t2[:,b,i,j],optimize=True)
    x -= np.einsum('me,m,e->',h_voov[b,:,k,:],t2[a,d,:,l],t2[:,c,i,j],optimize=True)
    x += np.einsum('me,m,e->',h_voov[b,:,l,:],t2[a,d,:,k],t2[:,c,i,j],optimize=True)
    return x 

def val_2(a,b,c,d,i,j,k,l,h_oooo,t2):
    return np.einsum('mn,m,n->',h_oooo[:,:,i,j],t2[a,d,:,l],t2[b,c,:,k],optimize=True)

def val_3(a,b,c,d,i,j,k,l,h_vvvv,t2):
    return np.einsum('ef,e,f->',h_vvvv[a,d,:,:],t2[:,b,i,j],t2[c,:,k,l],optimize=True)

def build_antisym_diagram(a,b,c,d,i,j,k,l,val,H,t2):

    x =   val(a,b,c,d,i,j,k,l,H,t2) # (1)
    x -=  val(b,a,c,d,i,j,k,l,H,t2) # -(ab)
    x -=  val(c,b,a,d,i,j,k,l,H,t2) # -(ac)
    x -=  val(a,d,c,b,i,j,k,l,H,t2) # -(bd)
    x -=  val(a,b,d,c,i,j,k,l,H,t2) # -(cd)
    x +=  val(b,a,d,c,i,j,k,l,H,t2) # +(ab)(cd)

    x -=  val(a,b,c,d,k,j,i,l,H,t2) # -(ik)
    x +=  val(b,a,c,d,k,j,i,l,H,t2) # +(ik)(ab)
    x +=  val(c,b,a,d,k,j,i,l,H,t2) # +(ik)(ac)
    x +=  val(a,d,c,b,k,j,i,l,H,t2) # +(ik)(bd)
    x +=  val(a,b,d,c,k,j,i,l,H,t2) # +(ik)(cd)
    x -=  val(b,a,d,c,k,j,i,l,H,t2) # -(ik)(ab)(cd)

    x -=  val(a,b,c,d,l,j,k,i,H,t2) # -(il)
    x +=  val(b,a,c,d,l,j,k,i,H,t2) # +(il)(ab)
    x +=  val(c,b,a,d,l,j,k,i,H,t2) # +(il)(ac)
    x +=  val(a,d,c,b,l,j,k,i,H,t2) # +(il)(bd)
    x +=  val(a,b,d,c,l,j,k,i,H,t2) # +(il)(cd)
    x -=  val(b,a,d,c,l,j,k,i,H,t2) # -(il)(ab)(cd)

    x -=  val(a,b,c,d,i,k,j,l,H,t2) # -(jk)
    x +=  val(b,a,c,d,i,k,j,l,H,t2) # +(jk)(ab)
    x +=  val(c,b,a,d,i,k,j,l,H,t2) # +(jk)(ac)
    x +=  val(a,d,c,b,i,k,j,l,H,t2) # +(jk)(bd)
    x +=  val(a,b,d,c,i,k,j,l,H,t2) # +(jk)(cd)
    x -=  val(b,a,d,c,i,k,j,l,H,t2) # -(jk)(ab)(cd)

    x -=  val(a,b,c,d,i,l,k,j,H,t2) # -(jl)
    x +=  val(b,a,c,d,i,l,k,j,H,t2) # +(jl)(ab)
    x +=  val(c,b,a,d,i,l,k,j,H,t2) # +(jl)(ac)
    x +=  val(a,d,c,b,i,l,k,j,H,t2) # +(jl)(bd)
    x +=  val(a,b,d,c,i,l,k,j,H,t2) # +(jl)(cd)
    x -=  val(b,a,d,c,i,l,k,j,H,t2) # -(jl)(ab)(cd)

    x +=  val(a,b,c,d,k,l,i,j,H,t2) # +(ik)(jl)
    x -=  val(b,a,c,d,k,l,i,j,H,t2) # -(ik)(jl)(ab)
    x -=  val(c,b,a,d,k,l,i,j,H,t2) # -(ik)(jl)(ac)
    x -=  val(a,d,c,b,k,l,i,j,H,t2) # -(ik)(jl)(bd)
    x -=  val(a,b,d,c,k,l,i,j,H,t2) # -(ik)(jl)(cd)
    x +=  val(b,a,d,c,k,l,i,j,H,t2) # +(ik)(jl)(ab)(cd)

    return x 

def build_MM24(cc_t,H2,sys):
    
    print('MM(2,4) construction... ')
    t_start = time.time()

    t2 = cc_t['t2']

    # A(jl/i/k)A(bc/a/d)
    D1 = -np.einsum('amie,bcmk,edjl->abcdijkl',H2['voov'],t2,t2,optimize=True)

    # Factorize the antisymmetrizer: A(pq/r/s) = A(pq/rs)A(rs)

    # A(jl/i/k) = A(jl/ik)A(ik)
    D1 += -permute(D1,[1,2,3,4,7,6,5,8])
    D1 += -permute(D1,[1,2,3,4,6,5,7,8]) - permute(D1,[1,2,3,4,5,7,6,8]) - permute(D1,[1,2,3,4,8,6,7,5])\
    -permute(D1,[1,2,3,4,5,6,8,7]) + permute(D1,[1,2,3,4,6,5,8,7])

    # A(bc/a/d) = A(bc/ad)A(ad)
    D1 += -permute(D1,[4,2,3,1,5,6,7,8])
    D1 += -permute(D1,[2,1,3,4,5,6,7,8]) - permute(D1,[3,2,1,4,5,6,7,8]) - permute(D1,[1,4,3,2,5,6,7,8])\
    -permute(D1,[1,2,4,3,5,6,7,8]) + permute(D1,[2,1,4,3,5,6,7,8])

    # (ij/kl)(bc/ad)
    D2 = np.einsum('mnij,adml,bcnk->abcdijkl',H2['oooo'],t2,t2,optimize=True)

    D2 += -permute(D2,[1,2,3,4,7,6,5,8]) - permute(D2,[1,2,3,4,8,6,7,5]) - permute(D2,[1,2,3,4,5,7,6,8])\
    -permute(D2,[1,2,3,4,5,8,7,6]) + permute(D2,[1,2,3,4,7,8,5,6])

    D2 += -permute(D2,[2,1,3,4,5,6,7,8]) - permute(D2,[3,2,1,4,5,6,7,8]) - permute(D2,[1,4,3,2,5,6,7,8])\
    -permute(D2,[1,2,4,3,5,6,7,8]) + permute(D2,[2,1,4,3,5,6,7,8])

    # (jk/il)(ab/cd)
    D3 = np.einsum('abef,fcjk,edil->abcdijkl',H2['vvvv'],t2,t2,optimize=True)

    D3 += -permute(D3,[1,2,3,4,6,5,7,8]) - permute(D3,[1,2,3,4,7,6,5,8]) - permute(D3,[1,2,3,4,5,6,8,7])\
    -permute(D3,[1,2,3,4,5,8,7,6]) + permute(D3,[1,2,3,4,6,5,8,7])

    D3 += -permute(D3,[3,2,1,4,5,6,7,8]) - permute(D3,[4,2,3,1,5,6,7,8]) - permute(D3,[1,3,2,4,5,6,7,8])\
    -permute(D3,[1,4,3,2,5,6,7,8]) + permute(D3,[3,4,1,2,5,6,7,8])

    MM24 = D1 + D2 + D3

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return MM24, D1, D2, D3

def build_L4(cc_t,ints,sys):

    print('Approximate L4 construction... ')
    
    t_start = time.time()

    V = ints['v']
    l2 = cc_t['l2']

    L4 = np.einsum('ijab,cdkl->abcdijkl',V['oovv'],l2,optimize=True)

    L4 += -permute(L4,[1,2,3,4,7,6,5,8]) - permute(L4,[1,2,3,4,8,6,7,5])\
        -permute(L4,[1,2,3,4,5,7,6,8]) - permute(L4,[1,2,3,4,5,8,7,6])\
        +permute(L4,[1,2,3,4,7,8,5,6]) # A(ij/kl)
    L4 += -permute(L4,[3,2,1,4,5,6,7,8]) - permute(L4,[4,2,3,1,5,6,7,8])\
        -permute(L4,[1,3,2,4,5,6,7,8]) - permute(L4,[1,4,3,2,5,6,7,8])\
        +permute(L4,[3,4,1,2,5,6,7,8]) # A(ab/cd)

    t_end = time.time()
    minutes, seconds = divmod(t_end-t_start, 60)
    print('finished in ({:0.2f}m  {:0.2f}s)'.format(minutes,seconds))

    return L4

def triples_diagonal(V_oovv,t2,sys):

    d3_V = lambda a,i,b: -np.dot(V_oovv[i,:,a,b].T,t2[a,b,i,:])
    d3_O = lambda a,i,j:  np.dot(V_oovv[i,j,a,:].T,t2[a,:,i,j])

    D3_V = np.zeros((sys['Nunocc'],sys['Nocc'],sys['Nunocc']))
    D3_O = np.zeros((sys['Nunocc'],sys['Nocc'],sys['Nocc']))

    for a in range(sys['Nunocc']):
        for i in range(sys['Nocc']):
            for b in range(sys['Nunocc']):
                D3_V[a,i,b] = d3_V(a,i,b)
    for a in range(sys['Nunocc']):
        for i in range(sys['Nocc']):
            for j in range(sys['Nocc']):
                D3_O[a,i,j] = d3_O(a,i,j)

    D3 = {'V' : D3_V, 'O' : D3_O}

    return D3

def cc_energy(cc_t,ints):

    f = ints['f']
    v = ints['v']
    t1 = cc_t['t1']
    t2 = cc_t['t2']

    Ecorr = 0.0
    Ecorr += np.einsum('me,em->',f['ov'],t1,optimize=True)
    Ecorr += 0.5*np.einsum('mnef,em,fn->',v['oovv'],t1,t1,optimize=True)
    Ecorr += 0.25*np.einsum('mnef,efmn->',v['oovv'],t2,optimize=True)

    return Ecorr


def permute(x,perm_list):
    str1 = ['a','b','c','d','i','j','k','l']
    str2 = ''.join([str1[x-1] for x in perm_list]) 
    str1 = ''.join(s for s in str1)
    contr = str1+'->'+str2
    #print(contr)
    return np.einsum(contr,x,optimize=True)
    #return x.transpose(tuple([y-1 for y in perm_list]))

def main(matfile):

    from scipy.io import loadmat

    # Load in integrals and T and Lambda vectors from spinorbital CCSD calculation
    print('')
    print('TEST SUBROUTINE:')
    print('Loading Matlab .mat file from {}'.format(matfile))
    print('Test case: H2O-2.0Re / DZ basis set')
    print('')

    data_dict = loadmat(matfile)
    t1 = data_dict['t1']
    t2 = data_dict['t2']
    l1 = data_dict['lambda1']
    l2 = data_dict['lambda2']
    h_oo = data_dict['H1_oo']
    h_vv = data_dict['H1_vv']
    h_voov = data_dict['H2_voov']
    h_oooo = data_dict['H2_oooo']
    h_vvvv = data_dict['H2_vvvv']
    V_oovv = data_dict['V_oovv']
    f_oo = data_dict['f_oo']
    f_vv = data_dict['f_vv']
    f_ov = data_dict['f_ov']

    Escf = data_dict['Escf'][0,0]
    Nocc = data_dict['Nocc'][0,0]
    Nunocc = data_dict['Nunocc'][0,0]
    Norb = Nocc + Nunocc

    v = {'oovv' : V_oovv}
    f = {'oo' : f_oo, 'vv' : f_vv, 'ov' : f_ov}
    H1 = {'oo' : h_oo, 'vv' : h_vv}
    H2 = {'voov' : h_voov, 'oooo' : h_oooo, 'vvvv' : h_vvvv}
    
    cc_t = {'t1' : t1, 't2' : t2, 'l1' : l1, 'l2' : l2}
    ints = {'v' : v, 'f' : f, 'Escf' : Escf}
    sys = {'Nocc' : Nocc, 'Nunocc' : Nunocc, 'Norb' : Norb}

    Ecorr = cc_energy(cc_t,ints)

    print('CCSD correlation energy = {} Eh'.format(Ecorr))

    Ecrcc24,delta24 = crcc24(cc_t,H1,H2,ints,sys)
    
    
    return

if __name__ == '__main__':
    matfile = '/home2/gururang/CCpy/tests/h2o-2.0-dz-spinorb.mat'
    main(matfile)
