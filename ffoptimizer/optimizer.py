import os
import sys
import time
import shutil
import json
from collections import OrderedDict

import numpy as np
from lmfit import Parameters, Minimizer

from .ppf import PPF
from .db import DB
from .models import Task, Target, Result

from sqlalchemy import and_


class Optimizer():
    def __init__(self, db_file):
        self.db = DB(db_file)
        self.db.conn()
        self.CWD = os.getcwd()

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

    def reset_task(self, task_name):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        for target in task.targets:
            files = os.listdir(target.dir_base_npt)
            for filename in files:
                filepath = os.path.join(target.dir_base_npt, filename)
                if os.path.isdir(filepath) and not filepath.endswith('-1'):
                    try:
                        shutil.rmtree(filepath)
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

    def optimize(self, task_name, power_residual=1, torsions=None,
                 weight_expansivity=0, penalty_sigma=0, penalty_epsilon=0, penalty_charge=0):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        if task.iteration != 0:
            print('Error: This task has been optimized before')
            sys.exit(0)

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
            result = self.db.session.query(Result).filter(
                and_(
                    Result.task == task,
                    Result.parameter == str(params)
                )).first()
            if result is not None:
                R = result.residual
                if R is not None:
                    return json.loads(R)
            ###

            ### save ppf file and run NPT
            ppf = PPF(string=task.ppf)
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value
            ppf.set_nb_paras(paras)

            # TODO fit several torsions iteratively. More torsion to fit, more cycles. This is inefficient
            if torsions is not None:
                for i in range(len(torsions)):
                    for torsion in torsions:
                        print('Fit torsion based on new non-bonded parameters. Cycle %i / %i ...' % (
                            i + 1, len(torsions)))
                        print(torsion)
                        ppf.fit_torsion(torsion[0], torsion[1], torsion[2], torsion[3])

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
                    cmds = target.run_npt(ppf_out, paras)
                    ### save gtx_dirs and gtx_cmds for running jobs on gtx queue
                    if cmds != []:
                        gtx_dirs.append(target.dir)
                        gtx_cmds = cmds

                os.chdir(self.CWD)

                if gtx_dirs != []:
                    from .models import npt
                    commands_list = npt.gmx.generate_gpu_multidir_cmds(gtx_dirs, gtx_cmds)
                    npt.jobmanager.queue = 'gtx'
                    npt.jobmanager.nprocs = 2
                    for i, commands in enumerate(commands_list):
                        sh = os.path.join(task.dir, '_job.multi-%i.sh' % i)
                        npt.jobmanager.generate_sh(task.dir, commands, name='NPT-GTX-%i-%i' % (task.iteration, i),
                                                   sh=sh)
                        npt.jobmanager.submit(sh)

            while True:
                if task.npt_finished():
                    break
                else:
                    current_time = time.strftime('%m-%d %H:%M')
                    print(current_time + ' Job still running. Wait ...')
                    time.sleep(60)

            R_dens = []
            R_hvap = []
            targets = task.targets.all()
            for target in targets:
                dens, hvap = target.get_npt_result()
                R_dens.append((dens - target.density) / target.density * 100 * target.wDens)  # deviation  percent
                R_hvap.append((hvap - target.hvap) / target.hvap * 100 * target.wHvap)  # deviation percent
            R = R_dens + R_hvap
            os.chdir(self.CWD)

            ### thermal expansivity
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
                if k.endswith('bi'):
                    res = v.value - adj_nb_paras[k]
                else:
                    res = (v.value - adj_nb_paras[k]) / adj_nb_paras[k]
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
            for k, v in params.items():
                txt += '%10.5f  %s\n' % (v.value, k)
            txt += '\n%8s %8s %10s %8s %8s %s %s\n' % (
                'RESIDUAL', 'Property', 'Deviation', 'Expt.', 'Weight', 'Molecule', 'SMILES')
            for i, r in enumerate(R_dens):
                target = targets[i]
                prop = 'density'
                weight = target.wDens
                txt += '%8.2f %8s %8.2f %% %8.3f %8.2f %s %s\n' \
                       % (r, prop, r / weight, target.density, weight, target.name, target.smiles)
            for i, r in enumerate(R_hvap):
                target = targets[i]
                prop = 'hvap'
                weight = target.wHvap
                txt += '%8.2f %8s %8.2f %% %8.1f %8.2f %s %s\n' \
                       % (r, prop, r / weight, target.hvap, weight, target.name, target.smiles)

            if weight_expansivity != 0:
                for i, r in enumerate(R_expa):
                    target = targets[i * 2]
                    prop = 'expan'
                    weight = weight_expansivity
                    txt += '%8.2f %8s %8.2f %% %8s %8.2f %s %s\n' \
                           % (r, prop, r / weight, '', weight, target.name, target.smiles)

            for i, r in enumerate(R_pena):
                prop = 'penalty'
                txt += '%8.2f %8s %10s\n' % (r, prop, list(params.keys())[i])

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return R

            # self.R = R
            # new_R = [r ** power_residual for r in R]
            # return new_R

        def jacobian(params: Parameters):
            ### if result exist in database, ignore calculation
            result = self.db.session.query(Result).filter(
                and_(
                    Result.task == task,
                    Result.parameter == str(params)
                )).first()
            if result is not None:
                J = result.jacobian
                if J is not None:
                    return json.loads(J)
            ###

            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value

            J_dens = []
            J_hvap = []
            targets = task.targets.all()
            for target in targets:
                dDdp_list, dHdp_list = target.get_dDens_dHvap_list_from_paras(paras)
                J_dens.append([i / target.density * 100 * target.wDens for i in dDdp_list])  # deviation  percent
                J_hvap.append([i / target.hvap * 100 * target.wHvap for i in dHdp_list])  # deviation  percent
            J = J_dens + J_hvap
            os.chdir(self.CWD)

            ### thermal expansivity
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
                if k.endswith('bi'):
                    d = 1
                else:
                    d = 1 / adj_nb_paras[k]
                penalty = get_penalty_for_para(k)
                J_pena.append(d * penalty * np.sqrt(len(J_dens)))
            J_pena = [list(a) for a in np.diag(J_pena)]  # convert list to diagonal matrix
            J += J_pena

            ### save result to database
            result = self.db.session.query(Result).filter(
                and_(
                    Result.task == task,
                    Result.iteration == task.iteration
                )).first()

            result.jacobian = json.dumps(J)
            self.db.session.commit()
            ###

            ### write Jacobian to log
            txt = '\nJACOBIAN MATRIX:\n'
            for k in params.keys():
                txt += '%10s' % k
            txt += '\n'
            for i, row in enumerate(J_dens):
                name = targets[i].name
                prop = 'density'
                for item in row:
                    txt += '%10.2f' % item
                txt += ' %8s %s\n' % (prop, name)
            for i, row in enumerate(J_hvap):
                name = targets[i].name
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
                name = list(params.keys())[i]
                prop = 'penalty'
                for item in row:
                    txt += '%10.2f' % item
                txt += ' %8s %s\n' % (prop, name)

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return J

            # new_J = []
            # for i, j_list in enumerate(J):
            #     new_J.append([power_residual * self.R[i] ** (power_residual - 1) * j for j in j_list])
            # return new_J

        def callback(params: Parameters, iter: int, res: [float]):
            print('Wait for 3 seconds ...')
            time.sleep(3)

        ppf = PPF(string=task.ppf)
        adj_nb_paras = ppf.get_adj_nb_paras()
        params = Parameters()
        for k, v in adj_nb_paras.items():
            bound = PPF.get_bound_for_para(k)
            params.add(k, value=v, min=bound[0], max=bound[1])

        minimize = Minimizer(residual, params, iter_cb=callback)
        result = minimize.leastsq(Dfun=jacobian, ftol=0.005)
        print(result.lmdif_message, '\n')

        return result.params

    def plot(self, task_name, iterations=None):
        try:
            import pylab
        except:
            print('matplotlib not found. will not plot.')
            return

        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        if iterations is None:
            iterations = (1, task.iteration)

        props = dict()
        for target in task.targets:
            mol = target.name
            if not mol in props.keys():
                props[mol] = {'smiles': target.smiles,
                              'T': [],
                              'dens': OrderedDict([('expt', [])]),
                              'hvap': OrderedDict([('expt', [])])
                              }
            props[mol]['T'].append(target.T)
            props[mol]['dens']['expt'].append(target.density)
            props[mol]['hvap']['expt'].append(target.hvap)

            for i in iterations:
                if i not in props[mol]['dens'].keys():
                    props[mol]['dens'][i] = []
                    props[mol]['hvap'][i] = []

                density, hvap = target.get_npt_result(i)

                props[mol]['dens'][i].append(density)
                props[mol]['hvap'][i].append(hvap)

        os.chdir(self.CWD)
        pylab.rcParams.update({'font.size': 12})
        for mol, prop in props.items():
            pylab.figure(figsize=(6, 8))
            pylab.subplot(211)
            for i, points in prop['dens'].items():
                if i == 'expt':
                    marker = '--'
                elif i == 1:
                    marker = 'x'
                else:
                    marker = 'o'
                pylab.plot(prop['T'], points, marker, label=i)
            y_mean = np.mean(prop['dens']['expt'])
            pylab.ylim(y_mean - 0.2, y_mean + 0.2)
            pylab.legend()
            pylab.title('Density %s %s (g/mL)' % (mol, prop['smiles']))

            pylab.subplot(212)
            for i, points in prop['hvap'].items():
                if i == 'expt':
                    marker = '--'
                elif i == 1:
                    marker = 'x'
                else:
                    marker = 'o'
                pylab.plot(prop['T'], points, marker, label=i)
            y_mean = np.mean(prop['hvap']['expt'])
            pylab.ylim(y_mean - 20, y_mean + 20)
            pylab.legend()
            pylab.title('HVap %s %s (kJ/mol)' % (mol, prop['smiles']))
            pylab.savefig('%s-%s.png' % (task.name, mol))
