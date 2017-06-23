#!/usr/bin/env python3
# coding=utf-8

import os
import sys

sys.path.append('..')

from ffoptimizer import Optimizer

if __name__ == '__main__':
    optimizer = Optimizer(db_file='ffoptimizer.db')

    cmd = sys.argv[1]
    if cmd == 'init':
        optimizer.init_task(sys.argv[2], sys.argv[3], sys.argv[4],
                            sys.argv[5])  # task_name, data_file, ppf_file, work_dir

    elif cmd == 'optimize':
        task_name = sys.argv[2]
        # optimizer.optimize(ppf_file, 0.15)
        optimizer.optimize(task_name, 0.15, 'MP2-C6.qmd', 'C6.msd', 'TORS C5 C8 C11 C14 500.0 180.0 10.0 18')

    elif cmd == 'remove':
        optimizer.remove_task(sys.argv[2])

    elif cmd == 'list':
        optimizer.list_task()

    elif cmd == 'plot':
        optimizer.plot(sys.argv[2], tuple(map(int, sys.argv[3:])))

    else:
        print('Unknown command')
