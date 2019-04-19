import json
import os
import shutil
import sys
import time
import copy
from collections import OrderedDict

import numpy as np
from lmfit import Parameters, Minimizer

from .db import DB
from .models import Task, Target, Result, PPF
from .para_tool import get_bound_for_para, replace_para_name, restore_para_name


class Optimizer():
    def __init__(self, db_file):
        self.db = DB(db_file)
        self.db.conn()
        self.CWD = os.getcwd()
        self.n_parallel = 8
        self.drde_dict = {}
        self.max_iter = 9

    def init_task(self, task_name, data_file, ppf_file, work_dir):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is not None:
            print('Error: Task %s already exist' % task_name)
            sys.exit(0)

        task = Task()
        task.name = task_name
        task.cwd = os.path.abspath(work_dir)
        ppf = PPF(ppf_file)
        task.ppf = str(ppf)

        self.db.session.add(task)
        self.db.session.flush()

        with open(data_file) as f:
            lines = f.read().splitlines()
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue
            words = line.split()

            target = Target(task=task)
            target.name = words[0]
            target.smiles = words[1]
            target.T = int(words[2])
            target.P = int(words[3])
            target.density = float(words[4])
            target.wDens = float(words[5])
            target.hvap = float(words[6])
            target.wHvap = float(words[7])

            # TODO dielectric
            try:
                target.dielectric = float(words[8])
            except:
                pass

            target.calc_n_mol()
            self.db.session.add(target)

        try:
            os.makedirs(task.dir)
            self.db.session.commit()
        except:
            self.db.session.rollback()
            raise

    def list_task(self):
        for task in self.db.session.query(Task):
            print(task.name, task.dir)

    def reset_task(self, task_name, ppf_file=None):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        if ppf_file is not None:
            cycle = 0
            ppf = PPF(ppf_file)
            task.ppf = str(ppf)
        else:
            cycle = 1

        for target in task.targets:
            del_files = []
            try:
                files = os.listdir(target.dir_base_npt)
                del_files += list(map(lambda x: os.path.join(target.dir_base_npt, x), files))
            except:
                pass
            try:
                files = os.listdir(target.dir_base_vacuum)
                del_files += list(map(lambda x: os.path.join(target.dir_base_vacuum, x), files))
            except:
                pass

            for fpath in del_files:
                if os.path.isdir(fpath):
                    dircycle = int(fpath.split('-')[-1])
                    if dircycle > cycle:
                        try:
                            shutil.rmtree(fpath)
                        except Exception as e:
                            print(str(e))

        task.results.delete()
        task.iteration = 0
        self.db.session.commit()

    def remove_task(self, task_name):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        try:
            shutil.rmtree(task.dir)
        except Exception as e:
            print(str(e))

        task.targets.delete()
        task.results.delete()
        self.db.session.delete(task)
        self.db.session.commit()

    def optimize(self, task_name, max_iter=None, torsions=None, modify_torsions=None,
                 weight_expansivity=0,
                 penalty_sigma=0, penalty_epsilon=0, penalty_charge=0,
                 drde_atoms={}):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        if task.iteration != 0:
            print('Error: This task has been optimized before')
            sys.exit(0)

        task.check_need_hvap_files()

        if max_iter != None:
            self.max_iter = max_iter

        LOG = os.path.join(self.CWD, '%s.log' % task_name)

        self.R = []

        def get_penalty_for_para(key):
            if key.endswith('r0'):
                penalty = penalty_sigma
            elif key.endswith('e0'):
                penalty = penalty_epsilon
            elif key.endswith('bi'):
                penalty = penalty_charge
            else:
                penalty = 0

            return penalty

        def residual(params: Parameters):
            ### if result exist in database, ignore calculation
            result = self.db.session.query(Result).filter(Result.task == task) \
                .filter(Result.parameter == str(params)).first()
            if result is not None:
                R = result.residual
                if R is not None:
                    return json.loads(R)
            ###

            ### save ppf file and run NPT
            ppf = PPF(string=task.ppf)
            paras = OrderedDict()
            for k, v in params.items():
                print(v)
                paras[restore_para_name(k)] = v.value
            ppf.set_nb_paras(paras)

            # TODO Fit several torsions one by one
            if torsions is not None and len(torsions) > 0:
                from .config import Config
                print('Fit torsion based on new non-bonded parameters')
                for n, torsion in enumerate(torsions):
                    print(torsion)
                    ppf.fit_torsion(Config.DFF_ROOT, torsion[0], torsion[1], torsion[2], torsion[3],
                                    dfi_name='fit_torsion-%i-%i' % (task.iteration + 1, n))
            if modify_torsions is not None:
                for torsion in modify_torsions:
                    ppf.modify_torsion(torsion[0], torsion[1], torsion[2])

            ### new iteration
            task.iteration += 1
            self.db.session.commit()

            ppf_out = os.path.join(self.CWD, '%s-%i.ppf' % (task.name, task.iteration))
            ppf.write(ppf_out)

            if not task.npt_started():
                ### save gtx_dirs and gtx_cmds for running jobs on gtx queue
                gtx_dirs = []
                gtx_cmds = []
                ###
                for target in task.targets:
                    if not target.need_npt:
                        continue
                    if target.npt_started():
                        continue
                    cmds = target.run_npt(ppf_out, paras, drde_dict=self.drde_dict)
                    ### save gtx_dirs and gtx_cmds for running jobs on gtx queue
                    if cmds != []:
                        gtx_dirs.append(target.dir_npt)
                        gtx_cmds = cmds

                os.chdir(self.CWD)

                if gtx_dirs != []:
                    from .models import npt, jobmanager
                    commands_list = npt.gmx.generate_gpu_multidir_cmds(gtx_dirs, gtx_cmds,
                                                                       n_parallel=self.n_parallel,
                                                                       n_gpu=jobmanager.ngpu,
                                                                       n_procs=jobmanager.nprocs)
                    for i, commands in enumerate(commands_list):
                        sh = os.path.join(task.dir, '_job.npt-%i.sh' % i)
                        jobmanager.generate_sh(task.dir, commands,
                                               name='%s-%i-%i' % (task.name, task.iteration, i), sh=sh)
                        jobmanager.submit(sh)

            if not task.vacuum_started():
                gtx_dirs = []
                gtx_cmds = []
                for target in task.targets:
                    if not target.need_vacuum:
                        continue
                    if target.vacuum_started():
                        continue
                    cmds = target.run_vacuum(ppf_out, paras, drde_dict=self.drde_dict)
                    ### save gtx_dirs and gtx_cmds for running jobs on gtx queue
                    if cmds != []:
                        gtx_dirs.append(target.dir_vacuum)
                        gtx_cmds = cmds

                os.chdir(self.CWD)

                if gtx_dirs != []:
                    from .models import vacuum, jobmanager
                    commands_list = vacuum.gmx.generate_gpu_multidir_cmds(gtx_dirs, gtx_cmds,
                                                                          n_parallel=self.n_parallel,
                                                                          n_gpu=jobmanager.ngpu,
                                                                          n_procs=jobmanager.nprocs)
                    for i, commands in enumerate(commands_list):
                        sh = os.path.join(task.dir, '_job.vacuum-%i.sh' % i)
                        jobmanager.generate_sh(task.dir, commands,
                                               name='%s-%i-VAC%i' % (task.name, task.iteration, i), sh=sh)
                        jobmanager.submit(sh)

            while True:
                if task.npt_finished() and task.vacuum_finished():
                    break
                else:
                    current_time = time.strftime('%m-%d %H:%M')
                    print(current_time + ' Job still running. Wait ...')
                    time.sleep(60)

            Dens = []
            Hvap = []
            R_dens = []
            R_hvap = []
            targets = task.targets.all()
            for target in targets:
                if target.wDens > 1E-4:
                    dens = target.get_density()
                    R_dens.append((dens - target.density) / target.density * 100 * target.wDens)  # deviation  percent
                    Dens.append(dens)
                if target.wHvap > 1E-4:
                    hvap = target.get_hvap()
                    R_hvap.append((hvap - target.hvap) / target.hvap * 100 * target.wHvap)  # deviation percent
                    Hvap.append(hvap)
            R = R_dens + R_hvap
            os.chdir(self.CWD)

            ### expansivity
            if weight_expansivity != 0:
                R_expa = []
                for i_mol in range(len(targets) // 2):
                    target_T1 = targets[2 * i_mol]
                    target_T2 = targets[2 * i_mol + 1]
                    res_Kt = ((target_T1.sim_dens - target_T2.sim_dens) / (target_T1.density - target_T2.density) - 1) \
                             * 100 * weight_expansivity
                    R_expa.append(res_Kt)

                R += R_expa

            # parameter penalty
            R_pena = []
            for k, v in params.items():
                if k.endswith('r0') or k.endswith('e0'):
                    res = (v.value - adj_nb_paras[restore_para_name(k)]) / adj_nb_paras[restore_para_name(k)]
                elif k.endswith('bi'):
                    res = v.value - adj_nb_paras[restore_para_name(k)]
                else:
                    res = v.value
                penalty = get_penalty_for_para(k)
                R_pena.append(res * penalty * np.sqrt(len(R_dens)))
            R += R_pena

            ### save result to database
            result = Result(task=task)
            result.iteration = task.iteration
            result.ppf = str(ppf)
            result.parameter = str(params)
            result.residual = json.dumps(R)
            self.db.session.add(result)
            self.db.session.commit()
            ###

            ### write current parameters and residual to log
            txt = '\nITERATION %i, RSQ %.2f\n' % (task.iteration, np.sum(list(map(lambda x: x ** 2, R))))
            txt += '\nPARAMETERS:\n'
            for k, v in self.drde_dict.items():
                txt += '%10.5f  %-12s  Fixed\n' % (v, k)
            for k, v in params.items():
                txt += '%10.5f  %-12s  %10.5f\n' % (v.value, restore_para_name(k), init_params[k])
            txt += '\n%8s %8s %10s %8s %8s %8s %3s %3s %s %s\n' % (
                'RESIDUAL', 'Property', 'Deviation', 'Expt.', 'Simu.', 'Weight', 'T', 'P', 'Molecule', 'SMILES')

            targets_dens = task.targets.filter(Target.wDens > 1E-4).all()
            for i, r in enumerate(R_dens):
                target = targets_dens[i]
                prop = 'density'
                weight = target.wDens
                txt += '%8.2f %8s %8.2f %% %8.3f %8.3f %8.2f %3i %3i %s %s\n' % (
                    r, prop, r / weight, target.density, Dens[i], weight, target.T, target.P, target.name,
                    target.smiles)

            targets_hvap = task.targets.filter(Target.wHvap > 1E-4).all()
            for i, r in enumerate(R_hvap):
                target = targets_hvap[i]
                prop = 'hvap'
                weight = target.wHvap
                txt += '%8.2f %8s %8.2f %% %8.1f %8.1f %8.2f %3i %3i %s %s\n' % (
                    r, prop, r / weight, target.hvap, Hvap[i], weight, target.T, target.P, target.name,
                    target.smiles)

            if weight_expansivity != 0:
                for i, r in enumerate(R_expa):
                    target = targets[i * 2]
                    prop = 'expan'
                    weight = weight_expansivity
                    txt += '%8.2f %8s %8.2f %% %8s %8s %8.2f %3s %3s %s %s\n' % (
                        r, prop, r / weight, '', '', weight, '', '', target.name, target.smiles)

            for i, r in enumerate(R_pena):
                prop = 'penalty'
                k = list(params.keys())[i]
                txt += '%8.2f %8s %10s %8s %8s %8.2f\n' % (r, prop, k, '', '', get_penalty_for_para(k))

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return R

        def jacobian(params: Parameters):
            ### if result exist in database, ignore calculation
            result = self.db.session.query(Result).filter(Result.task == task) \
                .filter(Result.parameter == str(params)).first()
            if result is not None:
                J = result.jacobian
                if J is not None:
                    return json.loads(J)
            ###

            paras = OrderedDict()
            for k, v in params.items():
                paras[restore_para_name(k)] = v.value

            J_dens = []
            J_hvap = []
            targets = task.targets.all()
            for target in targets:
                if target.wDens > 1E-4:
                    dDdp_list = target.get_dDens_list_from_paras(paras)
                    J_dens.append([i / target.density * 100 * target.wDens for i in dDdp_list])  # deviation  percent
                if target.wHvap > 1E-4:
                    dHdp_list = target.get_dHvap_list_from_paras(paras)
                    J_hvap.append([i / target.hvap * 100 * target.wHvap for i in dHdp_list])  # deviation  percent
            J = J_dens + J_hvap
            os.chdir(self.CWD)

            ### expansivity
            if weight_expansivity != 0:
                J_expa = []
                for i_mol in range(len(targets) // 2):
                    target_T1 = targets[2 * i_mol]
                    target_T2 = targets[2 * i_mol + 1]
                    dExpa = (target_T1.dDdp_array - target_T2.dDdp_array) / (target_T1.density - target_T2.density) \
                            * 100 * weight_expansivity
                    J_expa.append(list(dExpa))

                J += J_expa

            ### parameter penalty
            J_pena = []
            for k, v in params.items():
                if k.endswith('r0') or k.endswith('e0'):
                    d = 1 / adj_nb_paras[restore_para_name(k)]
                else:
                    d = 1
                penalty = get_penalty_for_para(k)
                J_pena.append(d * penalty * np.sqrt(len(J_dens)))
            J_pena = [list(a) for a in np.diag(J_pena)]  # convert list to diagonal matrix
            J += J_pena

            ### save result to database
            result = self.db.session.query(Result).filter(Result.task == task) \
                .filter(Result.iteration == task.iteration).first()

            result.jacobian = json.dumps(J)
            self.db.session.commit()
            ###

            ### write Jacobian to log
            txt = '\nJACOBIAN MATRIX:\n'
            for k in params.keys():
                txt += '%10s' % restore_para_name(k)
            txt += '\n'

            targets_dens = task.targets.filter(Target.wDens > 1E-4).all()
            for i, row in enumerate(J_dens):
                name = targets_dens[i].name
                prop = 'density'
                for item in row:
                    txt += '%10.2f' % item
                txt += ' %8s %s\n' % (prop, name)

            targets_hvap = task.targets.filter(Target.wHvap > 1E-4).all()
            for i, row in enumerate(J_hvap):
                name = targets_hvap[i].name
                prop = 'hvap'
                for item in row:
                    txt += '%10.2f' % item
                txt += ' %8s %s\n' % (prop, name)

            if weight_expansivity != 0:
                for i, row in enumerate(J_expa):
                    name = targets[2 * i].name
                    prop = 'expan'
                    for item in row:
                        txt += '%10.2f' % item
                    txt += ' %8s %s\n' % (prop, name)

            for i, row in enumerate(J_pena):
                name = restore_para_name(list(params.keys())[i])
                prop = 'penalty'
                for item in row:
                    txt += '%10.2f' % item
                txt += ' %8s %s\n' % (prop, name)

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return J

        def callback(params: Parameters, iter: int, res: [float]):
            print(task.iteration, self.max_iter)
            if task.iteration >= self.max_iter:
                print('Max iter reached. Abort the optimization')
                sys.exit()
            print('Wait for 3 seconds ...')
            time.sleep(3)

        ppf = PPF(string=task.ppf)
        adj_nb_paras = ppf.get_adj_nb_paras()
        params = Parameters()
        for k, v in adj_nb_paras.items():
            bound = get_bound_for_para(k)
            params.add(replace_para_name(k), value=v, min=bound[0], max=bound[1])

        ### drde_atoms for temperature dependence
        for k, v in drde_atoms.items():
            bound = get_bound_for_para(k)
            params.add(k, value=v, min=bound[0], max=bound[1])

        init_params = copy.copy(params)
        minimize = Minimizer(residual, params, iter_cb=callback)
        result = minimize.leastsq(Dfun=jacobian, ftol=0.001)
        print(result.lmdif_message, '\n')

        return result.params
