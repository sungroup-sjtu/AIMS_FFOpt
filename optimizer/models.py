import copy
import math
import os
import sys
import shutil
from collections import OrderedDict

import pybel
import panedr
from pandas import Series
from sqlalchemy import Column, Integer, Text, Float, String

from functools import partial

NotNullColumn = partial(Column, nullable=False)

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata

from .ppf import PPF, get_delta_for_para
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


class Target(Base):
    __tablename__ = 'target'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200))
    smiles = NotNullColumn(Text)
    n_mol = Column(Integer, nullable=True)
    T = NotNullColumn(Integer)
    P = NotNullColumn(Integer)
    density = NotNullColumn(Float)
    hvap = NotNullColumn(Float)
    wDensity = NotNullColumn(Float)
    wHvap = NotNullColumn(Float)
    iteration = NotNullColumn(Integer, default=0)

    def __repr__(self):
        return '<Target: %s %s %i>' % (self.name, self.smiles, self.T)

    def calc_n_mol(self, n_atoms=3000, n_mol=0):
        py_mol = pybel.readstring('smi', self.smiles)
        py_mol.addh()
        self.n_mol = math.ceil(n_atoms / len(py_mol.atoms))
        if self.n_mol < n_mol:
            self.n_mol = n_mol

    @property
    def dir_base_npt(self):
        return os.path.join(Config.WORK_DIR, 'NPT-%s' % (self.name))

    @property
    def dir_child(self):
        return '%i-%i-%i' % (self.T, self.P, self.iteration)

    def run_npt(self, ppf_file=None, paras_diff: OrderedDict = None):
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

        cd_or_create_and_cd(os.path.basename('%s' % ppf_file)[:-4])

        npt.msd = '../init.msd'
        npt.export(ppf=ppf_file, minimize=True)

        cd_or_create_and_cd(self.dir_child)

        npt.jobmanager.refresh_preferred_queue()
        commands = npt.prepare(model_dir='..', T=self.T, P=self.P, jobname='NPT-%s-%i' % (self.name, self.T),
                               dt=0.002, nst_eq=int(3E5), nst_run=int(2E5), nst_trr=250, nst_xtc=250)

        nprocs = npt.jobmanager.nprocs
        if paras_diff is not None:
            commands.append('export GMX_MAXCONSTRWARN=-1')
            for k in paras_diff.keys():
                msd_out = '_tmp.msd'
                ppf_out = '_tmp.ppf'
                for i in [-1, 1]:
                    basename = 'diff-%s.%i' % (k, i)
                    top_out = basename + '.top'
                    top_out_hvap = basename + '-hvap.top'

                    paras = copy.copy(paras_diff)
                    paras[k] += get_delta_for_para(k) * i
                    ppf = PPF(ppf_file)
                    ppf.set_lj_para(paras)
                    ppf.write(ppf_out)

                    shutil.copy('../../init.msd', msd_out)
                    npt.dff.set_charge([msd_out], ppf_out)
                    npt.dff.export_gmx(msd_out, ppf_out, gro_out='_tmp.gro', top_out=top_out)

                    npt.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='diff.mdp', nstxtcout=0)
                    cmd = npt.gmx.grompp(mdp='diff.mdp', top=top_out, tpr_out=basename + '.tpr', get_cmd=True)
                    commands.append(cmd)
                    cmd = npt.gmx.mdrun(name=basename, nprocs=nprocs, rerun='npt.trr', get_cmd=True)
                    commands.append(cmd)

                    npt.gmx.generate_top_for_hvap(top_out, top_out_hvap)

                    cmd = npt.gmx.grompp(mdp='diff.mdp', top=top_out_hvap, tpr_out=basename + '-hvap.tpr', get_cmd=True)
                    commands.append(cmd)
                    cmd = npt.gmx.mdrun(name=basename + '-hvap', nprocs=nprocs, rerun='npt.trr', get_cmd=True)
                    commands.append(cmd)

        commands.append('touch _finished_')
        npt.jobmanager.generate_sh(os.getcwd(), commands, name='NPT-%s-%i-%i' % (self.name, self.T, self.iteration))
        npt.run()

    def get_npt_result(self, subdir, iteration=None) -> (float, float):
        if iteration is None:
            iteration = self.iteration
        os.chdir(self.dir_base_npt)
        os.chdir(subdir)
        os.chdir('%i-%i-%i' % (self.T, self.P, iteration))
        print(os.getcwd())

        df = panedr.edr_to_df('npt.edr')
        density = df.Density.mean() / 1000
        df = panedr.edr_to_df('hvap.edr')
        hvap = 8.314 * self.T / 1000 - df.Potential.mean() / self.n_mol
        return density, hvap

    def get_dDens_dHvap_from_paras(self, ppf_file, paras: OrderedDict):
        # read density and Hvap series
        df = panedr.edr_to_df('npt.edr')
        self.dens_series_npt: Series = df.Density
        df = panedr.edr_to_df('hvap.edr')
        self.hvap_series_npt: Series = 8.314 * self.T / 1000 - df.Potential / self.n_mol

        dDens = []
        dHvap = []
        for k in paras.keys():
            dD, dH = self.get_dDens_dHvap_from_para(ppf_file, k)
            dDens.append(dD)
            dHvap.append(dH)
        return dDens, dHvap

    def get_dDens_dHvap_from_para(self, ppf_file, k) -> (float, float):
        os.chdir(self.dir_base_npt)
        subdir = os.path.basename(ppf_file)[:-4]
        os.chdir(subdir)
        os.chdir(self.dir_child)

        # energy and Hvap after diff
        df = panedr.edr_to_df('diff-%s.1.edr' % k)
        pene_series_diff_p = df.Potential

        df = panedr.edr_to_df('diff-%s.1-hvap.edr' % k)
        hvap_series_diff_p = 8.314 * self.T / 1000 - df.Potential / self.n_mol

        df = panedr.edr_to_df('diff-%s.-1.edr' % k)
        pene_series_diff_n = df.Potential

        df = panedr.edr_to_df('diff-%s.-1-hvap.edr' % k)
        hvap_series_diff_n = 8.314 * self.T / 1000 - df.Potential / self.n_mol

        # calculate the derivative series dA/dp
        delta = get_delta_for_para(k)
        dPene_series: Series = (pene_series_diff_p - pene_series_diff_n) / delta / 2
        dHvap_series: Series = (hvap_series_diff_p - hvap_series_diff_n) / delta / 2

        # extract out the required density and hvap
        dens_series = self.dens_series_npt.loc[dPene_series.index]
        hvap_series = self.hvap_series_npt.loc[dPene_series.index]

        # calculate the derivative dA/dp according to ForceBalance
        densXdPene = dens_series * dPene_series
        hvapXdPene = hvap_series * dPene_series

        dDdp = -1 / 8.314 / self.T * (densXdPene.mean() - dens_series.mean() * dPene_series.mean())
        dHdp = dHvap_series.mean() - 1 / 8.314 / self.T * (hvapXdPene.mean() - hvap_series.mean() * dPene_series.mean())
        return dDdp, dHdp

    def npt_finished(self, ppf_file) -> bool:
        subdir = os.path.basename(ppf_file)[:-4]
        log_finished = os.path.join(self.dir_base_npt, subdir, self.dir_child, '_finished_')
        if os.path.exists(log_finished):
            return True

        return False

    def npt_started(self, ppf_file) -> bool:
        subdir = os.path.basename(ppf_file)[:-4]
        sh_job = os.path.join(self.dir_base_npt, subdir, self.dir_child, jobmanager.sh)
        if os.path.exists(sh_job):
            return True

        return False

    def clear_npt_result(self, ppf_file):
        subdir = os.path.basename(ppf_file)[:-4]
        dir = os.path.join(self.dir_base_npt, subdir, self.dir_child)
        log_finished = os.path.join(dir, '_finished_')
        sh_job = os.path.join(dir, jobmanager.sh)
        try:
            os.remove(log_finished)
            shutil.move(sh_job, sh_job + '.bak')
        except:
            pass


class Result(Base):
    __tablename__ = 'result'
    id = NotNullColumn(Integer, primary_key=True)
    ppf = NotNullColumn(String(200))
    parameter = Column(Text, nullable=True)
    residual = Column(Text, nullable=True)
    jacobian = Column(Text, nullable=True)
