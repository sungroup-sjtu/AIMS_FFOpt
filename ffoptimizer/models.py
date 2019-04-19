import copy
import math
import os
import shutil
import sys
from collections import OrderedDict
from functools import partial

import numpy as np
import panedr
import pybel
from sqlalchemy import Column, Integer, Text, Float, String, ForeignKey
from sqlalchemy.orm import relationship

NotNullColumn = partial(Column, nullable=False)

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata

from .para_tool import get_delta_for_para
from .config import Config

sys.path.append(Config.MS_TOOLS_DIR)

from mstools.utils import cd_or_create_and_cd
from mstools.jobmanager import Local, Torque, Slurm
from mstools.simulation.gmx import Npt, NvtGas
from mstools.wrapper.ppf import PPF, delta_ppf

if Config.JOB_MANAGER == 'slurm':
    PBS = Slurm
elif Config.JOB_MANAGER == 'torque':
    PBS = Torque
else:
    raise Exception('Job manager not supported')

jobmanager = PBS(*Config.PBS_QUEUE, env_cmd=Config.PBS_ENV_CMD)
jobmanager.time = 1

kwargs = {'packmol_bin': Config.PACKMOL_BIN,
          'dff_root'   : Config.DFF_ROOT,
          'dff_table'  : Config.DFF_TABLE,
          'gmx_bin'    : Config.GMX_BIN,
          'gmx_mdrun'  : Config.GMX_MDRUN,
          'jobmanager' : jobmanager}
npt = Npt(**kwargs)
vacuum = NvtGas(**kwargs)


def wrapper_target(target_funcname_args):
    target, funcname, args = target_funcname_args
    func = getattr(target, funcname)
    return func(args)


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
            if target.need_npt and not target.npt_started():
                return False
        return True

    def npt_finished(self):
        for target in self.targets:
            if target.need_npt and not target.npt_finished():
                return False
        return True

    def vacuum_started(self):
        for target in self.targets:
            if target.need_vacuum and not target.vacuum_started():
                return False
        return True

    def vacuum_finished(self):
        for target in self.targets:
            if target.need_vacuum and not target.vacuum_finished():
                return False
        return True

    def check_need_hvap_files(self):
        for target in self.targets:
            if target.wHvap > 0 and not target.need_vacuum:
                self.need_hvap_files = True
                break
        else:
            self.need_hvap_files = False


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
    wDens = NotNullColumn(Float)
    hvap = NotNullColumn(Float)
    wHvap = NotNullColumn(Float)
    # TODO dielectric
    dielectric = NotNullColumn('st', Float, default=1)

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
    def dir_base_vacuum(self):
        return os.path.join(self.task.dir, 'VACUUM-%s' % (self.name))

    @property
    def dir_npt(self):
        return os.path.join(self.dir_base_npt, '%i-%i-%i' % (self.T, self.P, self.task.iteration))

    @property
    def dir_vacuum(self):
        return os.path.join(self.dir_base_vacuum, '%i-%i-%i' % (self.T, self.P, self.task.iteration))

    @property
    def need_npt(self) -> bool:
        return self.wDens > 0 or self.wHvap > 0

    @property
    def need_vacuum(self) -> bool:
        return self.wHvap > 0 and self.smiles.find('.') > -1

    def run_npt(self, ppf_file: str, paras_diff: OrderedDict, drde_dict: {}) -> [str]:
        cd_or_create_and_cd(self.dir_base_npt)

        if not os.path.exists('init.msd'):
            smiles_list = self.smiles.split('.')
            density = self.density if self.density > 0 else None
            npt.set_system(smiles_list, n_atoms=3000, n_mol_ratio=[1] * len(smiles_list), density=density)
            npt.build(export=False)

        cd_or_create_and_cd(self.dir_npt)

        shutil.copy('../init.msd', npt.msd)

        ### temperature dependence
        ppf_run = 'run.ppf'
        paras = copy.copy(paras_diff)
        for k, v in drde_dict.items():
            if k not in paras.keys():
                paras[k] = v
        delta_ppf(ppf_file, ppf_run, self.T, paras)
        ###

        npt.export(ppf=ppf_run)

        # TODO dielectric
        npt.gmx._DIELECTRIC = self.dielectric or 1

        # TODO Because of the float error in gmx edr file, MAKE SURE nst_edr equals nst_trr and nst_xtc
        commands = npt.prepare(T=self.T, P=self.P, jobname='NPT-%s-%i' % (self.name, self.T),
                               dt=0.002, nst_eq=int(3E5), nst_run=int(2E5), nst_edr=200, nst_trr=200, nst_xtc=200)
                               # dt=0.001, nst_eq=int(5E5), nst_run=int(5E5), nst_edr=500, nst_trr=500, nst_xtc=500)

        commands.insert(0, 'touch _started_')

        if paras_diff is not None:
            npt.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='diff.mdp', nstxtcout=0, restart=True)
            commands.append('export GMX_MAXCONSTRWARN=-1')
            nprocs = npt.jobmanager.nprocs

            import multiprocessing
            def worker(key, return_dict):
                worker_commands = []
                for i in [-1, 1]:
                    basename = 'diff%i.%s' % (i, key)
                    ppf_diff = basename + '.ppf'
                    msd_diff = basename + '.msd'
                    gro_diff = basename + '.gro'
                    top_diff = basename + '.top'
                    top_diff_hvap = basename + '-hvap.top'

                    ### temperature dependence
                    paras_delta = copy.copy(paras_diff)
                    paras_delta[key] += get_delta_for_para(key) * i
                    for fuck, v in drde_dict.items():
                        if fuck not in paras_delta.keys():
                            paras_delta[fuck] = v
                    delta_ppf(ppf_file, ppf_diff, self.T, paras_delta)
                    ###

                    shutil.copy(npt.msd, msd_diff)

                    npt.dff.set_charge([msd_diff], ppf_diff, dfi_name=basename)
                    npt.dff.export_gmx(msd_diff, ppf_diff, gro_out=gro_diff, top_out=top_diff, dfi_name=basename)

                    cmd = npt.gmx.grompp(mdp='diff.mdp', top=top_diff, tpr_out=basename + '.tpr', get_cmd=True)
                    worker_commands.append(cmd)
                    cmd = npt.gmx.mdrun(name=basename, nprocs=nprocs, rerun='npt.trr', get_cmd=True)
                    worker_commands.append(cmd)

                    if self.task.need_hvap_files:
                        npt.gmx.generate_top_for_hvap(top_diff, top_diff_hvap)

                        cmd = npt.gmx.grompp(mdp='diff.mdp', top=top_diff_hvap, tpr_out=basename + '-hvap.tpr',
                                             get_cmd=True)
                        worker_commands.append(cmd)
                        cmd = npt.gmx.mdrun(name=basename + '-hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.trr',
                                            get_cmd=True)
                        worker_commands.append(cmd)

                    os.remove(ppf_diff)
                    os.remove(msd_diff)
                    os.remove(gro_diff)
                    os.remove(basename + '.dfi')
                    os.remove(basename + '.dfo')

                return_dict[key] = worker_commands

            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            jobs = []
            for k in paras_diff.keys():
                p = multiprocessing.Process(target=worker, args=(k, return_dict))
                jobs.append(p)
                p.start()
            for p in jobs:
                p.join()
            for worker_commands in return_dict.values():
                commands += worker_commands

        commands.append('touch _finished_')
        jobmanager.generate_sh(os.getcwd(), commands, name='NPT-%s-%i-%i' % (self.name, self.T, self.task.iteration))

        if jobmanager.ngpu == 0:
            npt.run()
            commands = []

        return commands

    def run_vacuum(self, ppf_file: str, paras_diff: OrderedDict, drde_dict: {}) -> [str]:
        cd_or_create_and_cd(self.dir_base_vacuum)

        if not os.path.exists('init.msd'):
            smiles_list = self.smiles.split('.')
            vacuum.set_system(smiles_list, n_atoms=3000, n_mol_list=[1] * len(smiles_list), density=None)
            vacuum.length = 40
            vacuum.build(export=False)

        cd_or_create_and_cd(self.dir_vacuum)

        shutil.copy('../init.msd', vacuum.msd)

        ### temperature dependence
        ppf_run = 'run.ppf'
        paras = copy.copy(paras_diff)
        for k, v in drde_dict.items():
            if k not in paras.keys():
                paras[k] = v
        delta_ppf(ppf_file, ppf_run, self.T, paras)
        ###

        vacuum.export(ppf=ppf_run)

        # TODO dielectric
        vacuum.gmx._DIELECTRIC = self.dielectric or 1
        # TODO Because of the float error in gmx edr file, MAKE SURE nst_edr equals nst_trr
        commands = vacuum.prepare(T=self.T, jobname='NPT-%s-%i' % (self.name, self.T),
                                  dt=0.002, nst_eq=int(2E5), nst_run=int(5E5), nst_edr=100, nst_trr=100, nst_xtc=0)
                                  # dt=0.001, nst_eq=int(4E5), nst_run=int(10E5), nst_edr=500, nst_trr=500, nst_xtc=0)

        commands.insert(0, 'touch _started_')

        if paras_diff is not None:
            vacuum.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='diff.mdp', nstxtcout=0, restart=True)
            commands.append('export GMX_MAXCONSTRWARN=-1')
            nprocs = vacuum.jobmanager.nprocs

            import multiprocessing
            def worker(key, return_dict):
                worker_commands = []
                for i in [-1, 1]:
                    basename = 'diff%i.%s' % (i, key)
                    ppf_diff = basename + '.ppf'
                    msd_diff = basename + '.msd'
                    gro_diff = basename + '.gro'
                    top_diff = basename + '.top'

                    ### temperature dependence
                    paras_delta = copy.copy(paras_diff)
                    paras_delta[key] += get_delta_for_para(key) * i
                    for fuck, v in drde_dict.items():
                        if fuck not in paras_delta.keys():
                            paras_delta[fuck] = v
                    delta_ppf(ppf_file, ppf_diff, self.T, paras_delta)
                    ###

                    shutil.copy(vacuum.msd, msd_diff)

                    vacuum.dff.set_charge([msd_diff], ppf_diff, dfi_name=basename)
                    vacuum.dff.export_gmx(msd_diff, ppf_diff, gro_out=gro_diff, top_out=top_diff, dfi_name=basename)

                    cmd = vacuum.gmx.grompp(mdp='diff.mdp', top=top_diff, tpr_out=basename + '.tpr', get_cmd=True)
                    worker_commands.append(cmd)
                    cmd = vacuum.gmx.mdrun(name=basename, nprocs=nprocs, n_omp=nprocs, rerun='nvt.trr', get_cmd=True)
                    worker_commands.append(cmd)

                    os.remove(ppf_diff)
                    os.remove(msd_diff)
                    os.remove(gro_diff)
                    os.remove(basename + '.dfi')
                    os.remove(basename + '.dfo')

                return_dict[key] = worker_commands

            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            jobs = []
            for k in paras_diff.keys():
                p = multiprocessing.Process(target=worker, args=(k, return_dict))
                jobs.append(p)
                p.start()
            for p in jobs:
                p.join()
            for worker_commands in return_dict.values():
                commands += worker_commands

        commands.append('touch _finished_')
        jobmanager.generate_sh(os.getcwd(), commands, name='VACUUM-%s-%i-%i' % (self.name, self.T, self.task.iteration))

        if jobmanager.ngpu == 0:
            vacuum.run()
            commands = []

        return commands

    def get_density(self) -> float:
        os.chdir(self.dir_npt)
        print(os.getcwd())

        df = panedr.edr_to_df('npt.edr')
        density = df.Density.mean() / 1000  # convert to g/mL

        self.sim_dens = density  # save self.sim_dens for calculating expansivity
        return density

    def get_hvap(self) -> float:
        os.chdir(self.dir_npt)
        print(os.getcwd())

        if not self.need_vacuum:
            df = panedr.edr_to_df('hvap.edr')
            hvap = self.RT - df.Potential.mean() / self.n_mol
        else:
            df = panedr.edr_to_df('npt.edr')
            pe_liq = df.Potential.mean()
            os.chdir(self.dir_vacuum)
            print(os.getcwd())

            df = panedr.edr_to_df('nvt.edr')
            pe_gas = df.Potential.mean()
            hvap = self.RT + pe_gas - pe_liq / self.n_mol

        return hvap

    def get_dDens_list_from_paras(self, paras: OrderedDict):
        os.chdir(self.dir_npt)

        df = panedr.edr_to_df('npt.edr')
        # TODO Because of the float error in gmx edr file, the index in Series is erroneous. Convert to array
        self.dens_array = np.array(df.Density) / 1000  # convert to g/mL

        # dDdp_list = [self.get_dDens_from_para(k) for k in paras.keys()]
        from multiprocessing import Pool
        with Pool(len(paras)) as p:
            dDdp_list = p.map(wrapper_target, [(self, 'get_dDens_from_para', k) for k in paras.keys()])

        self.dDdp_array = np.array(dDdp_list)  # save dDdp_array for calculating expansivity
        return dDdp_list

    def get_dDens_from_para(self, k) -> (float, float):
        os.chdir(self.dir_npt)

        # energy and Hvap after diff
        try:
            df = panedr.edr_to_df('diff1.%s.edr' % k)
        except:
            raise Exception('File not exist: ' + os.path.abspath('diff1.%s.edr' % k))
        pene_array_diff_p = np.array(df.Potential)

        try:
            df = panedr.edr_to_df('diff-1.%s.edr' % k)
        except:
            raise Exception('File not exist: ' + os.path.abspath('diff-1.%s.edr' % k))
        pene_array_diff_n = np.array(df.Potential)

        # calculate the derivative series dA/dp
        delta = get_delta_for_para(k)
        dPene_array = (pene_array_diff_p - pene_array_diff_n) / delta / 2

        # calculate the derivative dA/dp according to ForceBalance
        # TODO To accurately calculate the covariant, using dens_array.mean() instead of dens_series.mean()
        dDdp = -1 / self.RT * ((self.dens_array * dPene_array).mean() - self.dens_array.mean() * dPene_array.mean())

        return dDdp

    def get_dHvap_list_from_paras(self, paras: OrderedDict):
        os.chdir(self.dir_npt)

        if not self.need_vacuum:
            df = panedr.edr_to_df('hvap.edr')
            self.hvap_array = self.RT - np.array(df.Potential) / self.n_mol
        else:
            df = panedr.edr_to_df('npt.edr')
            self.pe_liq_array = np.array(df.Potential)

            os.chdir(self.dir_vacuum)
            df = panedr.edr_to_df('nvt.edr')
            self.pe_gas_array = np.array(df.Potential)

        # dHdp_list = [self.get_dHvap_from_para(k) for k in paras.keys()]
        from multiprocessing import Pool
        with Pool(len(paras)) as p:
            dHdp_list = p.map(wrapper_target, [(self, 'get_dHvap_from_para', k) for k in paras.keys()])

        return dHdp_list

    def get_dHvap_from_para(self, k) -> (float, float):
        os.chdir(self.dir_npt)

        # energy and Hvap after diff
        try:
            df = panedr.edr_to_df('diff1.%s.edr' % k)
        except:
            raise Exception('File not exist: ' + os.path.abspath('diff1.%s.edr' % k))
        pene_array_diff_p = np.array(df.Potential)

        try:
            df = panedr.edr_to_df('diff-1.%s.edr' % k)
        except:
            raise Exception('File not exist: ' + os.path.abspath('diff-1.%s.edr' % k))
        pene_array_diff_n = np.array(df.Potential)

        # calculate the derivative series dA/dp
        delta = get_delta_for_para(k)
        dPene_array = (pene_array_diff_p - pene_array_diff_n) / delta / 2

        if not self.need_vacuum:
            try:
                df = panedr.edr_to_df('diff1.%s-hvap.edr' % k)
            except:
                raise Exception('File not exist: ' + os.path.abspath('diff1.%s-hvap.edr' % k))
            hvap_array_diff_p = self.RT - np.array(df.Potential) / self.n_mol

            try:
                df = panedr.edr_to_df('diff-1.%s-hvap.edr' % k)
            except:
                raise Exception('File not exist: ' + os.path.abspath('diff-1.%s-hvap.edr' % k))
            hvap_array_diff_n = self.RT - np.array(df.Potential) / self.n_mol

            dHvap_array = (hvap_array_diff_p - hvap_array_diff_n) / delta / 2

            dHdp = dHvap_array.mean() - 1 / self.RT * (
                    (self.hvap_array * dPene_array).mean() - self.hvap_array.mean() * dPene_array.mean())
        else:
            dELIQdp = dPene_array.mean() - 1 / self.RT * (
                    (self.pe_liq_array * dPene_array).mean() - self.pe_liq_array.mean() * dPene_array.mean())

            os.chdir(self.dir_vacuum)

            try:
                df = panedr.edr_to_df('diff1.%s.edr' % k)
            except:
                raise Exception('File not exist: ' + os.path.abspath('diff1.%s.edr' % k))
            pene_array_diff_p = np.array(df.Potential)

            try:
                df = panedr.edr_to_df('diff-1.%s.edr' % k)
            except:
                raise Exception('File not exist: ' + os.path.abspath('diff-1.%s.edr' % k))
            pene_array_diff_n = np.array(df.Potential)
            dPene_array = (pene_array_diff_p - pene_array_diff_n) / delta / 2

            dEGASdp = dPene_array.mean() - 1 / self.RT * (
                    (self.pe_gas_array * dPene_array).mean() - self.pe_gas_array.mean() * dPene_array.mean())

            dHdp = dEGASdp - dELIQdp / self.n_mol

        return dHdp

    def npt_started(self) -> bool:
        log_started = os.path.join(self.dir_npt, '_started_')
        if os.path.exists(log_started):
            return True

        return False

    def npt_finished(self) -> bool:
        log_finished = os.path.join(self.dir_npt, '_finished_')
        if os.path.exists(log_finished):
            return True

        return False

    def vacuum_started(self) -> bool:
        log_started = os.path.join(self.dir_vacuum, '_started_')
        if os.path.exists(log_started):
            return True

        return False

    def vacuum_finished(self) -> bool:
        log_finished = os.path.join(self.dir_vacuum, '_finished_')
        if os.path.exists(log_finished):
            return True

        return False


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
