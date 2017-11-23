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
        if line.strip() == '' or line[0] == '#':
            continue

        words = line.split()
        name = words[0]
        smiles = words[1]
        paras = words[3:]
        exp_data[smiles] = {'name':name, 'paras':paras}
    return exp_data

print('%-10s %-40s %4s %4s %8s %8s %8s %8s' %('#Name', 'SMILES', 'T', 'P', 'Density', 'Weight', 'Hvap', 'Weight'))

expt_data = load_expt_data(sys.argv[1])
for smiles, info in expt_data.items():
    name = info['name']
    paras = info['paras']

    t_fus = int(math.ceil(float(paras[0])))
    t_vap = int(math.floor(float(paras[1])))

    try:
        t_c = int(math.floor(float(paras[2])))
    except:
        t_c = t_vap + 200

    t_span = t_vap - t_fus

    t_list = [t_fus + 25, t_vap, t_c * 0.8]

    if t_span < 50:
        print('# %s: t_span < 50, not good' % name)
        t_list[0] = t_fus + 10

    if t_list[0] < 200:
        #print('# %s: t_fus < 200, not good' % name)
        t_list[0] = 200
        if t_list[1] < 250:
            print('# %s: t_vap < 250, not good' % name)
    if t_list[-1] > 600:
        print('# %s: t_c > 700, not good' % name)
        t_list[-1] = 600

    if len(sys.argv) > 2:
        T = int(sys.argv[2])
        t_list = [T]
        if t_fus > T or t_vap < T:
            print('# %s: %i not in range, not good' %(name, T))

    exp_density = []
    exp_hvap = []
    exp_p0 = []

    A,B,C,n = map(float, paras[3:7])
    for t in t_list:
        val = A*B**(-(1-t/C)**n)  # kg/m^3
        exp_density.append(val)

    A,B,n = map(float, paras[12:15])
    for t in t_list:
        val = A*(1-t/B)**n
        exp_hvap.append(val)

    A,B,C = map(float, paras[20:23])
    for t in t_list:
        val = 10**(A-B/(t-273+C))*0.001333
        exp_p0.append(val)

    for i, t in enumerate(t_list):
        wDensity = 1
        wHvap = 0.3
        if i % 3 == 2:
            wHvap = 0
        #try:
        print('%-10s %-40s %4i %4i %8.3f %8.2f %8.1f %8.2f' %(name, smiles, t, round(exp_p0[i]), exp_density[i], wDensity, exp_hvap[i], wHvap))
        #except:
            #pass

