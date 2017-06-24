#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')

from ffoptimizer import Optimizer

if __name__ == '__main__':
    optimizer = Optimizer(db_file='ffoptimizer.db')

    cmd = sys.argv[1]
    if cmd == 'init':
        optimizer.init_task(task_name=sys.argv[2],
                            data_file=sys.argv[3],
                            ppf_file=sys.argv[4],
                            work_dir=sys.argv[5])

    elif cmd == 'list':
        optimizer.list_task()

    elif cmd == 'reset':
        optimizer.reset_task(sys.argv[2])

    elif cmd == 'remove':
        optimizer.remove_task(sys.argv[2])

    elif cmd == 'optimize':
        # optimizer.optimize(task_name=sys.argv[2], wExpansivity=0)
        optimizer.optimize(task_name=sys.argv[2],
                           penalty={'r0': 1, 'e0': 1},
                           wExpansivity=0.15,
                           qmd='MP2-C6.qmd',
                           msd='C6.msd',
                           torsion='TORS C5 C8 C11 C14 500.0 180.0 10.0 18')

    elif cmd == 'plot':
        iterations = None
        if len(sys.argv) > 3:
            iterations = tuple(map(int, sys.argv[3:]))
        optimizer.plot(sys.argv[2], iterations)

    else:
        print('Unknown command')
