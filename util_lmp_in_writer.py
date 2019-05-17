import random


def write(fname='py_lmp_in', **kwargs):
    # analyze datafile
    lines = open(kwargs['data_fname']).readlines()
    for idx, line in enumerate(lines):
        if line.endswith('atoms\n'):
            atoms_count = int(line.split()[0])
        elif line.endswith('atom types\n'):
            atom_types_count = int(line.split()[0])
        elif line.endswith('bonds\n'):
            bonds_count = int(line.split()[0])
        elif line.endswith('bond types\n'):
            bond_types_count = int(line.split()[0])
        elif line.startswith('Atoms'):
            atoms_idx = idx + 2

    atom_lines = lines[atoms_idx : atoms_idx + atoms_count]
    atom_types = [int(line.split()[2]) for line in atom_lines]
    mmt_atoms = atom_types.count(1)
    mod_atoms = atom_types.count(2) + atom_types.count(3)
    poly_atoms = atom_types.count(4)

    assert(atom_types_count == 4)
    assert(bond_types_count == 5)
    assert(mmt_atoms + mod_atoms + poly_atoms == atoms_count)

    print('ok,', atoms_count, 'atoms')

    # write out lmp script content
    f = open(fname, 'w')
    tmppr = lambda *x : print(*x, file=f)

    tmppr('# Test created structure of MMT and modifier')
    tmppr('')
    tmppr('units         lj')
    tmppr('comm_modify   vel yes')
    tmppr('newton        on')
    tmppr('special_bonds lj/coul 1 1 1')
    tmppr('atom_style    full')
    tmppr('bond_style    harmonic')
    tmppr('neighbor      3 bin')
    tmppr('neigh_modify  delay 0 every 1 check no page 500000 one 50000')
    tmppr('')
    tmppr('read_data    ', kwargs['data_fname'])
    tmppr('mass          1 1')
    tmppr('mass          2 1.74')
    tmppr('mass          3 1.74')
    tmppr('mass          4 1.74')
    tmppr('')
    tmppr('pair_style    hybrid/overlay coul/cut {0} dpd 1 {1} {2}'.format(
        kwargs['pair_radius'], kwargs['dpd_cutoff'], random.randint(1, 100000000)))
    tmppr('')
    # a_ij,  gamma, r_ij
    tmppr('pair_coeff    1 1 dpd {0} {1}'.format(kwargs['a_cl_cl'],
        kwargs['gamma']))
    tmppr('pair_coeff    1 2 dpd {0} {1} {2}'.format(kwargs['a_cl_mh'],
        kwargs['gamma'], (kwargs['r_cl'] + kwargs['r_mh'])/2))
    tmppr('pair_coeff    1 3 dpd {0} {1} {2}'.format(kwargs['a_cl_mt'],
        kwargs['gamma'], (kwargs['r_cl'] + kwargs['r_mt'])/2))
    tmppr('pair_coeff    1 4 dpd {0} {1} {2}'.format(kwargs['a_cl_po'],
        kwargs['gamma'], (kwargs['r_cl'] + kwargs['r_po'])/2))
    tmppr('pair_coeff    2 2 dpd {0} {1} {2}'.format(kwargs['a_mh_mh'],
        kwargs['gamma'], kwargs['r_mh']))
    tmppr('pair_coeff    2 3 dpd {0} {1} {2}'.format(kwargs['a_mh_mt'],
        kwargs['gamma'], (kwargs['r_mh'] + kwargs['r_mt'])/2))
    tmppr('pair_coeff    2 4 dpd {0} {1} {2}'.format(kwargs['a_mh_po'],
        kwargs['gamma'], (kwargs['r_mh'] + kwargs['r_po'])/2))
    tmppr('pair_coeff    3 3 dpd {0} {1} {2}'.format(kwargs['a_mt_mt'],
        kwargs['gamma'], kwargs['r_mt']))
    tmppr('pair_coeff    3 4 dpd {0} {1} {2}'.format(kwargs['a_mt_po'],
        kwargs['gamma'], (kwargs['r_mt'] + kwargs['r_po'])/2))
    tmppr('pair_coeff    4 4 dpd {0} {1} {2}'.format(kwargs['a_po_po'],
        kwargs['gamma'], kwargs['r_po']))

    tmppr('')
    tmppr('pair_coeff    * * coul/cut')
    tmppr('bond_coeff    1 {0} {1}'.format(kwargs['cl_ed_k'], kwargs['cl_ed_l']))
    tmppr('bond_coeff    2 {0} {1}'.format(kwargs['cl_di_k'], kwargs['cl_di_l']))
    tmppr('bond_coeff    3 {0} {1}'.format(kwargs['mht_k'], kwargs['mht_l']))
    tmppr('bond_coeff    4 {0} {1}'.format(kwargs['mtt_k'], kwargs['mtt_l']))
    tmppr('bond_coeff    5 {0} {1}'.format(kwargs['po_l'], kwargs['po_l']))
    tmppr('')
    tmppr('thermo        10')
    tmppr('thermo_style  custom step density press pe ebond ecoul evdwl lz')
    tmppr('')
    tmppr('group         charged id 1:{0}'.format(mmt_atoms + mod_atoms))
    tmppr('group         soft_rared id {0}:{1}:100'.format(
        mmt_atoms + mod_atoms + 1, atoms_count))
    tmppr('group         imagable union charged soft_rared')
    tmppr('dump         ', 
          'd1 imagable image 1000 images_dump/*.jpg type type view 90 0')
    tmppr('')
    tmppr('timestep      0.0001')
    tmppr('fix          ',
          '1 all npt x 0 0 10000 y 0 0 10000 z 1000 1000 10000 temp 1 1 10000')
    tmppr('')
    tmppr('label loop')
    tmppr('variable a loop 500')
    tmppr('    run 5000 start 0 stop 2500000')
    tmppr('    write_data dpd.*.data')
    tmppr('    next a')
    tmppr('jump SELF loop')


if __name__ == '__main__':
    ###############
    # some test parameters

    big_a = 100
    small_a = 1

    big_k = 100
    small_k = 10

    pair_radius = 7
    dpd_cutoff = 2

    big_r = 1.72
    small_r = 1.72

    ###############

    kwargs = {
        'data_fname':   'periodic_composite.data',

        'pair_radius':  pair_radius,
        'dpd_cutoff':   dpd_cutoff,

        'a_cl_cl':      0,
        'a_cl_mh':      big_a,
        'a_cl_mt':      big_a,
        'a_cl_po':      big_a,
        'a_mh_mh':      small_a,
        'a_mh_mt':      small_a,
        'a_mh_po':      small_a,
        'a_mt_mt':      small_a,
        'a_mt_po':      small_a,
        'a_po_po':      small_a,

        'cl_ed_k':      big_k,
        'cl_ed_l':      2,
        'cl_di_k':      big_k,
        'cl_di_l':      2 * 3**0.5,
        'mht_k':        small_k,
        'mht_l':        2,
        'mtt_k':        small_k,
        'mtt_l':        2,
        'po_k':         small_k,
        'po_l':         2,

        'gamma':        5,

        'r_cl':         big_r,
        'r_mh':         small_r,
        'r_mt':         small_r,
        'r_po':         small_r,
    }

    write(**kwargs)
