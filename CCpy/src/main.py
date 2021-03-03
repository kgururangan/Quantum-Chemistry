import os
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
from ccp3_module import ccp3
from crcc24_module import crcc24, test_updates
from parser_module import parse_extcorr_file, get_p_space_array 

def get_file_names_with_strings(path,identifier):
    full_list = os.listdir(path)
    final_list = [str(x) for x in full_list if identifier in x]
    if len(final_list) > 0:
        return final_list[0]
    else: return ''

def main(args):

    work_dir = os.path.abspath(args.directory)+'/'
    nfroz = args.frozen

    gamess_file = work_dir + get_file_names_with_strings(work_dir,'log')
    onebody_file = work_dir + get_file_names_with_strings(work_dir,'onebody')
    twobody_file = work_dir + get_file_names_with_strings(work_dir,'twobody')
    extcorr_file = work_dir + get_file_names_with_strings(work_dir,'moe')
    pspace_file = work_dir + get_file_names_with_strings(work_dir,'.dat-p')

    sys = build_system(gamess_file,nfroz)
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

    #print('')
    #print('External Correction:')
    #print('=========================================')
    #cc_t = parse_extcorr_file(extcorr_file,ints,sys)
    #H1A,H1B,H2A,H2B,H2C = HBar_CCSD(cc_t,ints,sys)
    #cc_t = left_ccsd(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys,tol=1e-11,maxit=500,shift=0.3)
    #p_space = get_p_space_array(pspace_file,sys)
    #Ecrcc23,_ = ccp3(cc_t,p_space,H1A,H1B,H2A,H2B,H2C,ints,sys)

    #print('')
    #Efci = -75.9119458469
    #print('Error rel. FCI = {} mEh'.format( (Ecrcc23['D']-Efci)*1000 ))

    cc_t, Eccsd = ccsd(sys,ints,shift=0.0)
    H1A,H1B,H2A,H2B,H2C = HBar_CCSD(cc_t,ints,sys)
    cc_t = left_ccsd(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
   # test_updates(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys)
    _,E23 = crcc23(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys,flag_RHF=True)
    _,E24 = crcc24(cc_t,H1A,H1B,H2A,H2B,H2C,ints,sys,flag_RHF=True)

    # Expected CR-CC(2,4) results
    #print('')
    #print('Expected A correction : CR-CC(2,3) =  -0.0088261909 Eh  CR-CC(2,4) = -0.0014377831 Eh')
    #print('Expected B correction : CR-CC(2,3) =  -0.0076279185 Eh  CR-CC(2,4) = -0.0014636650 Eh')
    #print('Expected C correction : CR-CC(2,3) =  -0.0124386886 Eh  CR-CC(2,4) = +0.0009754553 Eh')
    #print('Expected D correction : CR-CC(2,3) =  -0.0117381126 Eh  CR-CC(2,4) = -0.0017982619 Eh')

    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='parser for Python CC implementation')
    parser.add_argument('directory',type=str,help='Path to working directory containing GAMESS log file and integrals')
    parser.add_argument('-f','--frozen',type=int,help='Number of frozen spatial orbitals',default=0)

    args = parser.parse_args()
    main(args)

