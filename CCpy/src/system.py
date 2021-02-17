import argparse
#import resultsFile

def build_system(gamess_file,Nfroz):
    import resultsFile

    gfile = resultsFile.getFile(gamess_file)
    print_gamess_info(gfile,gamess_file)

    gfile_options = gfile.options

    Nelec = gfile.num_elec
    Nocc_a = gfile.num_alpha
    Nocc_b = gfile.num_beta
    methods = gfile.methods
    scftyp = methods[0]
    mos = gfile.mo_sets[scftyp]
    point_group = gfile.point_group
    charge = gfile.charge
    multiplicity = gfile.multiplicity
    basis = gfile.basis
    occ_num = gfile.occ_num
    mo_sets = gfile.mo_sets
    closed_mos = gfile.closed_mos
    active_mos = gfile.active_mos
    virtual_mos = gfile.virtual_mos


    Norb = 0
    sym = []
    mo_energies = []
    mo_vecs = []
    for orb in mos:
        Norb += 1
        sym.append(orb.sym)
        mo_energies.append(orb.eigenvalue)
        mo_vecs.append(orb.vector)

    occ_num = [None] * Norb
    for i in range(Norb):
        if i in closed_mos:
            occ_num[i] = 2.0
        if i in active_mos:
            occ_num[i] = 1.0
        if i in virtual_mos:
            occ_num[i] = 0.0

    Nunocc_a = Norb  - Nocc_a
    Nunocc_b = Norb  - Nocc_b

    sys_t = {'Nelec' : Nelec-2*Nfroz, \
             'Nocc_a' : Nocc_a-Nfroz, \
             'Nocc_b' : Nocc_b-Nfroz, \
             'Nunocc_a' : Nunocc_a, \
             'Nunocc_b' : Nunocc_b, \
             'Norb' : Norb, \
             'Nfrz_a' : Nfroz, \
             'Nfrz_b' : Nfroz, \
             'sym' : sym, \
             'mo_energy' : mo_energies, \
             'mo_vector' : mo_vecs, \
             'point_group' : point_group, \
             'charge' : charge, \
             'multiplicity' : multiplicity, \
             'basis' : basis, \
             'mo_occ' : occ_num}

    return sys_t

def print_gamess_info(gfile,gamess_file):

    gfile_options = gfile.options

    if 'GBASIS' in gfile_options:
        basis_name = gfile_options['GBASIS']
    else:
        basis_name = 'User-defined'

    if gfile_options['ISPHER'] == '1':
        ispher = 'Spherical'
    else:
        ispher = 'Cartesian'

    print('')
    print('GAMESS Run Information:')
    print('-----------------------------------------------')
    print('  GAMESS file location : {}'.format(gamess_file))
    print('  Title of run : {}'.format(gfile.title))
    print('  GAMESS calculation performed by {} on {}'.format(gfile.author,gfile.date))
    print('  Calculation run on {}'.format(gfile.machine))
    print('  SCF Type : {}'.format(gfile_options['SCFTYP']))
    print('  Basis : {} ({})'.format(basis_name,ispher))
    print('')
    return

def main(args):

    sys_t = build_system(args.gamess_file,args.frozen)
    for key, value in sys_t.items():
        print(key,'->',value)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gamess_file',type=str,help='Path to GAMESS log file containing SCF calculation')
    parser.add_argument('-f','--frozen',type=int,default=0,help='Number of frozen spatial orbitals')
    args = parser.parse_args()
    main(args)
