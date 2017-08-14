#!/usr/bin/env python3
# coding=utf-8

import math

from collections import OrderedDict

def load_expt_data(filename):
    exp_data = OrderedDict()
    with open(filename) as f:
        lines = f.read().splitlines()
    for line in lines[1:]:
        if line.strip() == '' or line.startswith('#'):
            continue

        words = line.split()
        name = words[0]
        smiles = words[1]
        paras = list(map(float,words[3:]))
        exp_data[smiles] = {'name':name, 'paras':paras}
    return exp_data

print('%-10s %-40s %4s %4s %8s %8s %8s %8s' %('#Name', 'SMILES', 'T', 'P', 'Density', 'Weight', 'Hvap', 'Weight'))

expt_data = load_expt_data('mols.txt')
for smiles, info in expt_data.items():
    name = info['name']
    paras = info['paras']

    t_min = int(math.ceil(float(paras[0])))
    t_max = int(math.floor(float(paras[1])))
    t_mid = int(math.ceil((t_min + t_max) / 2))
    dt = t_max - t_min

    t_list = [t_min + 25, t_max - 25]
    if dt < 50:
        print('# %s: dt < 50, not good' % name)
        t_list = [t_min + 5, t_max - 5]
    elif dt < 100:
        print('# %s: dt < 100, not good' % name)
        t_list = [t_min + 5, t_max - 5]

    exp_density = []
    exp_hvap = []

    A,B,C,n = paras[2:6]
    for t in t_list:
        val = A*B**(-(1-t/C)**n)  # kg/m^3
        exp_density.append(val)

    A,B,n = paras[6:9]
    for t in t_list:
        val = A*(1-t/B)**n
        exp_hvap.append(val)

    for i, t in enumerate(t_list):
        wDensity = 1
        wHvap = 0.2
        print('%-10s %-40s %4i %4i %8.3f %8.2f %8.1f %8.2f' %(name, smiles, t, 1, exp_density[i], wDensity, exp_hvap[i], wHvap))
