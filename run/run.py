#!/usr/bin/env python3
# coding=utf-8

import os
import sys

sys.path.append('..')

from optimizer import Optimizer

if __name__ == '__main__':
    optimizer = Optimizer(db_file='ffoptimizer.db')

    cmd = sys.argv[1]
    if cmd == 'init':
        optimizer.init_db(sys.argv[2])

    elif cmd == 'optimize':
        ppf_file = os.path.abspath(sys.argv[2])
        # optimizer.optimize(ppf_file, 0.15)
        optimizer.optimize(ppf_file, 0.15, 'MP2-C6.qmd', 'C6.msd', 'TORS C5 C8 C11 C14 500.0 180.0 15.0 12')

    elif cmd == 'plot':
        ppfs = sys.argv[2:]
        optimizer.plot(ppfs, iteration=None)

    else:
        print('Unknown command')
