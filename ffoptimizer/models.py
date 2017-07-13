import copy
import math
import os
import sys
import shutil
from collections import OrderedDict

import numpy as np
import pybel
import panedr
from pandas import Series
from sqlalchemy import Column, Integer, Text, Float, String, ForeignKey
from sqlalchemy.orm import relationship

from functools import partial

NotNullColumn = partial(Column, nullable=False)

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata

from .ppf import PPF
from .config import Config

sys.path.append(Config.MS_TOOLS_DIR)

from mstools.utils import create_mol_from_smiles, cd_or_create_and_cd
from mstools.jobmanager import Local, Torque, Slurm
from mstools.simulation.gmx import Npt

if Config.JOB_MANAGER == 'local':
    jobmanager = Local(nprocs=Config.NPROC_PER_JOB)
elif Config.JOB_MANAGER == 'torque':
    jobmanager = Torque(queue_dict=Config.QUEUE_DICT)
elif Config.JOB_MANAGER == 'slurm':
    jobmanager = Slurm(queue_dict=Config.QUEUE_DICT)
else:
    raise Exception('Job manager not supported')

kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
          'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}
npt = Npt(**kwargs)


class Task(Base):
    __tablename__ = 'task'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200), unique=True)
    ppf = NotNullColumn(Text)
    cwd = NotNullColumn(String(200))
    iteration = NotNullColumn(Integer, default=0)

    targets = relationship('Target', lazy='dynamic')
    results = relationship('Result', lazy='dynamic')

    @property
    def dir(self):
        return os.path.join(self.cwd, self.name)

    def npt_started(self):
        for target in self.targets:
            if not target.npt_started():
                return False
        return True

    def npt_finished(self):
        for target in self.targets:
            if not target.npt_finished():
                return False
        return True


class Target(Base):
    __tablename__ = 'target'
    id = NotNullColumn(Integer, primary_key=True)
    task_id = NotNullColumn(Integer, ForeignKey(Task.id))
    name = NotNullColumn(String(200))
    smiles = NotNullColumn(Text)
    n_mol = Column(Integer, nullable=True)
    T = NotNullColumn(Integer)
    P = NotNullColumn(Integer)
    density = NotNullColumn(Float)
    hvap = NotNullColumn(Float)
    wDens = NotNullColumn(Float)
    wHvap = NotNullColumn(Float)

    task = relationship(Task)

    def __repr__(self):
        return '<Target: %s %s %i %i>' % (self.name, self.smiles, self.T, self.P)

    def calc_n_mol(self, n_atoms=3000, n_mol=0):
        py_mol = pybel.readstring('smi', self.smiles)
        py_mol.addh()
        self.n_mol = math.ceil(n_atoms / len(py_mol.atoms))
        if self.n_mol < n_mol:
            self.n_mol = n_mol

    @property
    def RT(self):
        return 8.314 * self.T / 1000  # kJ/mol

    @property
    def dir_base_npt(self):
        return os.path.join(self.task.dir, 'NPT-%s' % (self.name))

    @property
    def dir(self):
        return os.path.join(self.dir_base_npt, '%i-%i-%i' % (self.T, self.P, self.task.iteration))

    def run_npt(self, ppf_file=None, paras_diff: OrderedDict = None) -> [str]:
        cd_or_create_and_cd(self.dir_base_npt)

        if not os.path.exists('init.msd'):
            pdb = 'mol.pdb'
            mol2 = 'mol.mol2'
            py_mol = create_mol_from_smiles(self.smiles, pdb_out=pdb, mol2_out=mol2)
            mass = py_mol.molwt * self.n_mol
            length = (10 / 6.022 * mass / (self.density - 0.1)) ** (1 / 3)  # assume cubic box

            print('Build coordinates using Packmol: %s molecules ...' % self.n_mol)
            npt.packmol.build_box([pdb], [self.n_mol], 'init.pdb', length=length - 2, tolerance=1.7, silent=True)

            print('Create box using DFF ...')
            npt.dff.build_box_after_packmol([mol2], [self.n_mol], 'init.msd', mol_corr='init.pdb', length=length)

        cd_or_create_and_cd(self.dir)

        shutil.copy('../init.msd', npt.msd)

        ### temperature dependence of epsilon
        if paras_diff is not None:
            paras = copy.copy(paras_diff)
            for k, v in paras.items():
                if k.endswith('de'):
                    atype = k[:-3]
                    paras[atype + '_e0'] += v * (self.T - 298)
            ppf = PPF(ppf_file)
            ppf.set_nb_paras(paras)
            ppf_file = 'ff.ppf'
            ppf.write(ppf_file)
        npt.export(ppf=ppf_file, minimize=False)

        npt.jobmanager.refresh_preferred_queue()
        # TODO because of the float error in gmx edr file, MAKE SURE nst_edr equals nst_trr and nst_xtc
        commands = npt.prepare(T=self.T, P=self.P, jobname='NPT-%s-%i' % (self.name, self.T),
                               dt=0.002, nst_eq=int(3E5), nst_run=int(2E5), nst_edr=200, nst_trr=200, nst_xtc=200)

        commands.insert(0, 'touch _started_')

        if paras_diff is not None:
            npt.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='diff.mdp', nstxtcout=0, restart=True)
            commands.append('export GMX_MAXCONSTRWARN=-1')
            nprocs = npt.jobmanager.nprocs

            import multiprocessing
            def worker(k, return_dict):
                worker_commands = []
                for i in [-1, 1]:
                    basename = 'diff%i.%s' % (i, k)
                    ppf_diff = basename + '.ppf'
                    msd_diff = basename + '.msd'
                    gro_diff = basename + '.gro'
                    top_diff = basename + '.top'
                    top_diff_hvap = basename + '-hvap.top'

                    ### temperature dependence of epsilon
                    paras_delta = copy.copy(paras)
                    paras_delta[k] += PPF.get_delta_for_para(k) * i
                    ppf = PPF(ppf_file)
                    ppf.set_nb_paras(paras_delta)
                    ppf.write(ppf_diff)

                    shutil.copy(npt.msd, msd_diff)

                    npt.dff.set_charge([msd_diff], ppf_diff, dfi_name=basename)
                    npt.dff.export_gmx(msd_diff, ppf_diff, gro_out=gro_diff, top_out=top_diff, dfi_name=basename)

                    cmd = npt.gmx.grompp(mdp='diff.mdp', top=top_diff, tpr_out=basename + '.tpr', get_cmd=True)
                    worker_commands.append(cmd)
                    cmd = npt.gmx.mdrun(name=basename, nprocs=nprocs, rerun='npt.trr', get_cmd=True)
                    worker_commands.append(cmd)

                    npt.gmx.generate_top_for_hvap(top_diff, top_diff_hvap)

                    cmd = npt.gmx.grompp(mdp='diff.mdp', top=top_diff_hvap, tpr_out=basename + '-hvap.tpr',
                                         get_cmd=True)
                    worker_commands.append(cmd)
                    cmd = npt.gmx.mdrun(name=basename + '-hvap', nprocs=nprocs, rerun='npt.trr', get_cmd=True)
                    worker_commands.append(cmd)

                    os.remove(ppf_diff)
                    os.remove(msd_diff)
                    os.remove(gro_diff)
                    os.remove(basename + '.dfi')
                    os.remove(basename + '.dfo')

                return_dict[k] = worker_commands

            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            jobs = []
            for k in paras_diff.keys():
                ### temperature dependence of epsilon
                if k.endswith('de'):
                    continue
                p = multiprocessing.Process(target=worker, args=(k, return_dict))
                jobs.append(p)
                p.start()
            for p in jobs:
                p.join()
            for worker_commands in return_dict.values():
                commands += worker_commands

        commands.append('touch _finished_')
        npt.jobmanager.generate_sh(os.getcwd(), commands,
                                   name='NPT-%s-%i-%i' % (self.name, self.T, self.task.iteration))

        if npt.jobmanager.queue != 'gtx':
            npt.run()
            commands = []

        return commands

    def get_npt_result(self, iteration=None) -> (float, float):
        if iteration is None:
            iteration = self.task.iteration
        os.chdir(self.dir_base_npt)
        os.chdir('%i-%i-%i' % (self.T, self.P, iteration))
        print(os.getcwd())

        df = panedr.edr_to_df('npt.edr')
        density = df.Density.mean() / 1000  # convert to g/mL
        self.sim_dens = density  # save self.sim_dens for calculating thermal expansivity
        df = panedr.edr_to_df('hvap.edr')
        hvap = self.RT - df.Potential.mean() / self.n_mol

        return density, hvap

    def get_dDens_dHvap_list_from_paras(self, paras: OrderedDict):
        # read density and Hvap list
        os.chdir(self.dir)

        df = panedr.edr_to_df('npt.edr')
        self.dens_series_npt: Series = df.Density / 1000  # convert to g/mL
        df = panedr.edr_to_df('hvap.edr')
        self.hvap_series_npt: Series = self.RT - df.Potential / self.n_mol

        dDdp_list = []
        dHdp_list = []
        for k in paras.keys():
            ### temperature dependence of epsilon
            if k.endswith('de'):
                dDdp_list.append(dDdp_list[-1] * (self.T - 298))
                dHdp_list.append(dHdp_list[-1] * (self.T - 298))
            dDdp, dHdp = self.get_dDens_dHvap_from_para(k)
            dDdp_list.append(dDdp)
            dHdp_list.append(dHdp)
        self.dDdp_array = np.array(dDdp_list)  # save dDdp_array for calculating thermal expansivity

        return dDdp_list, dHdp_list

    def get_dDens_dHvap_from_para(self, k) -> (float, float):
        os.chdir(self.dir)

        # energy and Hvap after diff
        df = panedr.edr_to_df('diff1.%s.edr' % k)
        pene_array_diff_p = np.array(df.Potential)

        df = panedr.edr_to_df('diff1.%s-hvap.edr' % k)
        hvap_array_diff_p = np.array(self.RT - df.Potential / self.n_mol)

        df = panedr.edr_to_df('diff-1.%s.edr' % k)
        pene_array_diff_n = np.array(df.Potential)

        df = panedr.edr_to_df('diff-1.%s-hvap.edr' % k)
        hvap_array_diff_n = np.array(self.RT - df.Potential / self.n_mol)

        # calculate the derivative series dA/dp
        delta = PPF.get_delta_for_para(k)
        dPene_array = (pene_array_diff_p - pene_array_diff_n) / delta / 2
        dHvap_array = (hvap_array_diff_p - hvap_array_diff_n) / delta / 2

        # extract out the required density and hvap
        # TODO because of the float error in gmx edr file, the index in Series is errorous. Convert to array
        dens_array = np.array(self.dens_series_npt)
        hvap_array = np.array(self.hvap_series_npt)

        # calculate the derivative dA/dp according to ForceBalance
        densXdPene = dens_array * dPene_array
        hvapXdPene = hvap_array * dPene_array

        dDdp = -1 / self.RT * (densXdPene.mean() - dens_array.mean() * dPene_array.mean())
        dHdp = dHvap_array.mean() - 1 / self.RT * (hvapXdPene.mean() - hvap_array.mean() * dPene_array.mean())
        # !!! To accurately calculate the covariant, using dens_array.mean() instead of dens_series_npt.mean()

        return dDdp, dHdp

    def npt_finished(self) -> bool:
        log_finished = os.path.join(self.dir, '_finished_')
        if os.path.exists(log_finished):
            return True

        return False

    def npt_started(self) -> bool:
        log_started = os.path.join(self.dir, '_started_')
        if os.path.exists(log_started):
            return True

        return False

    def clear_npt_result(self):
        log_started = os.path.join(self.dir, '_started_')
        log_finished = os.path.join(self.dir, '_finished_')
        try:
            os.remove(log_started)
            os.remove(log_finished)
        except:
            pass


class Result(Base):
    __tablename__ = 'result'
    id = NotNullColumn(Integer, primary_key=True)
    task_id = NotNullColumn(Integer, ForeignKey(Task.id))
    iteration = NotNullColumn(Integer)
    ppf = NotNullColumn(String(200))
    parameter = Column(Text, nullable=True)
    residual = Column(Text, nullable=True)
    jacobian = Column(Text, nullable=True)

    task = relationship(Task)
