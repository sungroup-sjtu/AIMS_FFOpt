#!/usr/bin/env python3
# coding=utf-8

import os
import sys
from optimizer import Optimizer


if __name__ == '__main__':
    optimizer = Optimizer(db_file='ffoptimizer.db')

    cmd = sys.argv[1]
    if cmd == 'init':
        optimizer.init_db(sys.argv[2])

    elif cmd == 'optimize':
        ppf_file = os.path.abspath(sys.argv[2])
        params_out = optimizer.optimize_npt(ppf_file)

    elif cmd == 'npt':
        ppf_file = os.path.abspath(sys.argv[2])
        optimizer.npt(ppf_file)

    elif cmd == 'plot':
        ppfs = sys.argv[2:]
        optimizer.plot(ppfs)

    else:
        print('Unknown command')
