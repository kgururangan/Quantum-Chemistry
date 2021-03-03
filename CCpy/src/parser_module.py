import numpy as np
from scipy.io import FortranFile

def parse_moefile(moefile):
	with FortranFile(moefile,'r') as f90_file:
		first_line_reals = f90_file.read_reals(dtype=np.float)
		tvec = f90_file.read_reals(dtype=np.float)
	#	with FortranFile(moefile,'r') as f90_file:
	#		first_line_ints = f90_file.read_ints(dtype=np.int64)
		corr_energy = first_line_reals[1]
	#	corr_energy_terms = first_line_reals[2:10]
	#	tvec_size = first_line_ints[10] 				
	#	assert len(tvec) == tvec_size, \
	#		   "Length of the T vector doesn't match!"
	f90_file.close()

#	return tvec, corr_energy, corr_energy_terms
	return tvec, corr_energy

def extract_t_vector(tvec,sys):

    noa = sys['Nocc_a']
    nua = sys['Nunocc_a']
    nob = sys['Nocc_b']
    nub = sys['Nunocc_b']

    t1a = np.zeros((nua,noa))
    t1b = np.zeros((nub,nob))
    t2a = np.zeros((nua,nua,noa,noa))
    t2b = np.zeros((nua,nub,noa,nob))
    t2c = np.zeros((nub,nub,nob,nob))

    rec = 0
    for i in range(noa):
        for a in range(nua):
            t1a[a,i] = tvec[rec]
            rec += 1

    for i in range(nob):
        for a in range(nub):
            t1b[a,i] = tvec[rec]
            rec += 1

    for i in range(noa):
        for j in range(noa):
            for a in range(nua):
                for b in range(nua):
                    t2a[a,b,i,j] = tvec[rec]
                    rec += 1

    for i in range(noa):
        for j in range(nob):
            for a in range(nua):
                for b in range(nub):
                    t2b[a,b,i,j] = tvec[rec]
                    rec += 1

    for i in range(nob):
        for j in range(nob):
            for a in range(nub):
                for b in range(nub):
                    t2c[a,b,i,j] = tvec[rec]
                    rec += 1

    cc_t = {'t1a' : t1a, 't1b' : t1b, 't2a' : t2a, 't2b' : t2b, 't2c' : t2c}

    return cc_t

def get_p_space_array(p_file,sys):

    print('Parsing p space file from {}'.format(p_file))
    
    nua = sys['Nunocc_a']
    noa = sys['Nocc_a']
    nub = sys['Nunocc_b']
    nob = sys['Nocc_b']

    p_space = np.zeros((nua,nua,nua,noa,noa,noa))
    with open(p_file,'r') as P:
        for line in P.readlines():
            p = [int(x) for x in line.split()]
            p_space[p[0]-nob-1,p[1]-nob-1,p[2]-nob-1,p[3]-1,p[4]-1,p[5]-1] = 1
    return p_space

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


def parse_extcorr_file(moefile,ints,sys):
	
    print('Parsing extcorr file: {}'.format(moefile))
    tvec, Ecorr_file = parse_moefile(moefile)
    cc_t = extract_t_vector(tvec,sys)
    Ecorr_parse = calc_cc_energy(cc_t,ints)
    print('Correlation energy from file = {}'.format(Ecorr_file))
    print('Correlation energy from parsed file = {}'.format(Ecorr_parse))

    return cc_t


