#!/usr/bin/env python3
# coding=utf-8

import os, sys, json
import math
import pylab
import numpy as np

from collections import OrderedDict

def load_expt_data(filename):
    exp_data = OrderedDict()
    with open(filename) as f:
        lines = f.read().splitlines()
    for line in lines[1:]:
        if line.strip() == '':
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
    if dt <= 20:
        raise Exception('dt <= 20 not good')
    t_list = [t_min + 10, t_max - 10]

    if abs(t_min - 298) > abs(t_max - 298):
        weight = [1, 2]
    elif abs(t_min - 298) < abs(t_max - 298):
        weight = [2, 1]
    else:
        weight = [1.5, 1.5]

    #t_list = [t_min + 5, t_mid, t_max - 5]
    #t_list = list(range(t_min, t_max+1, int(dt / 4)))

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
        # wDensity = math.exp(-abs(t-298)/200) + 1
        wDensity = 1
        wHvap = 0.2
        print('%-10s %-40s %4i %4i %8.3f %8.2f %8.1f %8.2f' %(name, smiles, t, 1, exp_density[i], wDensity, exp_hvap[i], wHvap))

