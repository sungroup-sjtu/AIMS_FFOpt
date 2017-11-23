#!/usr/bin/env python3
# coding=utf-8

import os, sys

sys.path.append('..')

from ffoptimizer import Optimizer

if __name__ == '__main__':
    optimizer = Optimizer(db_file='ffoptimizer.db')

    cmd = sys.argv[1]

    if cmd == 'init':
        if len(sys.argv) == 6:
            work_dir = sys.argv[5]
        else:
            work_dir = os.getcwd()

        optimizer.init_task(task_name=sys.argv[2],
                            data_file=sys.argv[3],
                            ppf_file=sys.argv[4],
                            work_dir=work_dir)

    elif cmd == 'list':
        optimizer.list_task()

    elif cmd == 'reset':
        ppf_file = None
        if len(sys.argv) > 3:
            ppf_file = sys.argv[3]
        optimizer.reset_task(sys.argv[2], ppf_file)

    elif cmd == 'remove':
        optimizer.remove_task(sys.argv[2])

    elif cmd == 'optimize':
        if len(sys.argv) > 3:
            optimizer.n_parallel = int(sys.argv[3])

        # torsion parameters to be optimized
        torsions = [
            ('C2-MP2.qmd', 'C2.msd', 'TORS H5 C1 C2 H6 500.0 60.0 15.0 12', 'h_1, c_4, c_4, h_1'),
            ('C6-MP2.qmd', 'C6.msd', 'TORS C5 C8 C11 C14 500.0 180.0 10.0 18', 'c_4, c_4, c_4, c_4'),
        ]

        # torsion parameters to be modified
        modify_torsions = [
            # ('c_4, c_4, c_4, c_4', 2, 0.1),
        ]

        # temperature dependent LJ parameters
        optimizer.drde_dict = {
            'h_1_dr': -0.01,
            'h_1_de': 0.045,
            'c_4_dl': 0.020,
            'o_2_d2': 0.027,
        }

        optimizer.optimize(task_name=sys.argv[2],
                           torsions=torsions,
                           modify_torsions=modify_torsions,
                           #weight_expansivity=0.01,
                           penalty_sigma=2,
                           penalty_epsilon=0.4,
                           penalty_charge=0,
                           drde_atoms={
                               # 'all_dr': -0.01,
                               # 'all_de': 0.056
                           })

    else:
        print('Unknown command')
