import argparse
import numpy as np
from system import build_system
from integrals import get_integrals
from ccsd_module import ccsd
from HBar_module import HBar_CCSD
from ccsdt_module import ccsdt
#from eomccsd_module import test_updates
from left_ccsd_module import left_ccsd
from crcc23_module import crcc23
from crcc24_module import test_updates


def main(args):

    gamess_fid = args.gamess
    nfroz = args.frozen

    #gamess_fid = 'H2O-cc-pVDZ.log'
    #gamess_fid = 'rectangle-d2h.log'
    
    work_dir = "/Users/karthik/Dropbox/Hartree Fock/hartree_fock/v4/CCpy/tests/"
    matfile = work_dir + 'matlab_results.mat'
    #HBar_matfile = work_dir + 'matlab_HBar.mat'

    gamess_file = work_dir + gamess_fid
   # onebody_file = work_dir + 'onebody.inp'
   # twobody_file = work_dir + 'twobody.inp'
    onebody_file = '/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/onebody.inp'
    twobody_file = '/Users/karthik/Desktop/CC_matlab_tests/rectangle-d2h-pvdz/twobody.inp'

    sys = build_system(gamess_file,nfroz)
    print('')
    print('System Information:')
    print('-------------------------------------------------')
    print('  Number of correlated electrons = {}'.format(sys['Nelec']))
    print('  Number of frozen electrons = {}'.format(2*nfroz))
    print('  Charge = {}'.format(sys['charge']))
    print('  Point group = {}'.format(sys['point_group']))
    print('  Spin multiplicity of reference = {}'.format(sys['multiplicity']))
    print('')
    print('    MO #    Energy (Eh)   Symmetry    Occupation')
    print('-------------------------------------------------')
    for i in range(sys['Norb']):
        print('     {}       {}          {}         {}'.format(i+1,sys['mo_energy'][i],sys['sym'][i],sys['mo_occ'][i]))


    ints = get_integrals(onebody_file,twobody_file,sys)
    print('')
    print('  Nuclear Repulsion Energy = {} Eh'.format(ints['Vnuc']))
    print('  Reference Energy = {} Eh'.format(ints['Escf']))

    cc_t, Eccsd = ccsd(sys,ints,tol=1e-011)
    H1A,H1B,H2A,H2B,H2C = HBar_CCSD(cc_t,ints,sys)
    cc_t = left_ccsd(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys,shift=2,tol=1e-011)
    Ecrcc23 = crcc23(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)

    cc_t, Eccsdt = ccsdt(sys,ints,tol=1e-011)

    print('CR-CC(2,3)_A = {} mEh'.format( (Ecrcc23['A'] - Eccsdt) * 1000))
    print('CR-CC(2,3)_D = {} mEh'.format( (Ecrcc23['D'] - Eccsdt) * 1000))

    #test_updates(matfile,ints,sys)


    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='parser for Python CC implementation')
    parser.add_argument('gamess',type=str,help='GAMESS log file')
    parser.add_argument('-f','--frozen',type=int,help='Number of frozen spatial orbitals',default=0)

    args = parser.parse_args()
    main(args)

