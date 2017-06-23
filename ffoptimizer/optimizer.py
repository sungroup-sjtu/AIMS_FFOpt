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

    def init_task(self, task_name, data_file, ppf_file, cwd):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is not None:
            print('Error: Task %s already exist' % task_name)
            sys.exit(0)

        task = Task()
        task.name = task_name
        task.cwd = os.path.abspath(cwd)
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

    def remove_task(self, task_name):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        try:
            shutil.rmtree(task.dir)
        except:
            raise

        self.db.session.query(Target).filter(Target.task == task).delete()
        self.db.session.delete(task)
        self.db.session.commit()

    def list_task(self):
        for task in self.db.session.query(Task):
            print(task.name, task.dir)

    def optimize(self, task_name, wExpansivity, qmd=None, msd=None, torsion=None):
        task = self.db.session.query(Task).filter(Task.name == task_name).first()
        if task is None:
            print('Error: Task %s not exist' % task_name)
            sys.exit(0)

        if task.iteration != 0:
            print('Error: This task has been optimized before')
            sys.exit(0)

        LOG = os.path.join(self.CWD, '%s.log' % task_name)

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

            ### new iteration
            task.iteration += 1

            ### save ppf file and run NPT
            ppf = PPF(string=task.ppf)
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value
            ppf.set_nb_paras(paras)

            if torsion is not None:
                print('Fit torsion based on new non-bonded parameters...')
                ppf.fit_torsion(qmd, msd, torsion)

            ppf_out = os.path.join(self.CWD, '%s-%i.ppf' % (task.name, task.iteration))
            ppf.write(ppf_out)

            if not task.npt_started():
                for target in task.targets:
                    target.run_npt(ppf_out, paras)
            os.chdir(self.CWD)

            while True:
                if task.npt_finished():
                    break
                else:
                    current_time = time.strftime('%m-%d %H:%M')
                    print(current_time + ' Job still running. Wait ...')
                    time.sleep(60)

            R_dens = []
            R_hvap = []
            R_expa = []
            targets = task.targets.all()
            for target in targets:
                dens, hvap = target.get_npt_result()
                R_dens.append((dens - target.density) / target.density * 100 * target.wDens)  # deviation  percent
                R_hvap.append((hvap - target.hvap) / target.hvap * 100 * target.wHvap)  # deviation percent
            os.chdir(self.CWD)

            ### thermal expansivity
            for i_mol in range(len(targets) // 2):
                target_T1 = targets[2 * i_mol]
                target_T2 = targets[2 * i_mol + 1]
                res_Kt = ((target_T1.sim_dens - target_T2.sim_dens) / (target_T1.density - target_T2.density) - 1) \
                         * 100 * wExpansivity
                R_expa.append(res_Kt)

            R = R_dens + R_hvap + R_expa

            ### save result to database
            result = Result(task=task)
            result.iteration = task.iteration
            result.ppf = ppf
            result.parameter = str(params)
            result.residual = json.dumps(R)
            self.db.session.commit()
            ###

            ### write current parameters and residual to log
            txt = '\nITERATION: %i\n' % task.iteration
            txt += '\nRSQ: %.2f\n' % np.sum(list(map(lambda x: x ** 2, R)))
            txt += '\nPARAMETERS:\n'
            for k, v in params.items():
                txt += '%10.5f  %s\n' % (v.value, k)
            txt += '\nRESIDUE:\n'
            for i, r in enumerate(R_dens):
                target = targets[i]
                prop = 'density'
                weight = target.wDens
                txt += '%8.2f %10s %8.2f %8.2f %s\n' % (r, prop, weight, r / weight, target.name)
            for i, r in enumerate(R_hvap):
                target = targets[i]
                prop = 'hvap'
                weight = target.wHvap
                txt += '%8.2f %10s %8.2f %8.2f %s\n' % (r, prop, weight, r / weight, target.name)
            for i, r in enumerate(R_expa):
                target = targets[i * 2]
                prop = 'expan'
                weight = wExpansivity
                txt += '%8.2f %10s %8.2f %8.2f %s\n' % (r, prop, weight, r / weight, target.name)

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return R

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
            J_expa = []
            targets = task.targets.all()
            for target in targets:
                dDdp_list, dHdp_list = target.get_dDens_dHvap_list_from_paras(paras)
                J_dens.append([i / target.density * 100 * target.wDens for i in dDdp_list])  # deviation  percent
                J_hvap.append([i / target.hvap * 100 * target.wHvap for i in dHdp_list])  # deviation  percent
            os.chdir(self.CWD)

            ### thermal expansivity
            for i_mol in range(len(targets) // 2):
                target_T1 = targets[2 * i_mol]
                target_T2 = targets[2 * i_mol + 1]
                dExpa = (target_T1.dDdp_array - target_T2.dDdp_array) / (target_T1.density - target_T2.density) \
                        * 100 * wExpansivity
                J_expa.append(list(dExpa))

            J = J_dens + J_hvap + J_expa

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
                txt += ' %10s %s\n' % (prop, name)
            for i, row in enumerate(J_hvap):
                name = targets[i].name
                prop = 'hvap'
                for item in row:
                    txt += '%10.2f' % item
                txt += ' %10s %s\n' % (prop, name)
            for i, row in enumerate(J_expa):
                name = targets[2 * i].name
                prop = 'expan'
                for item in row:
                    txt += '%10.2f' % item
                txt += ' %10s %s\n' % (prop, name)

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return J

        def callback(params: Parameters, iter: int, res: [float]):
            print('Wait for 3 seconds ...')
            time.sleep(3)

        ppf = PPF(string=task.ppf)
        params = Parameters()
        for k, v in ppf.get_adj_nb_paras().items():
            bound = PPF.get_bound_for_para(k)
            params.add(k, value=v, min=bound[0], max=bound[1])

        minimize = Minimizer(residual, params, iter_cb=callback)
        result = minimize.leastsq(Dfun=jacobian, ftol=0.001)
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
        for target in self.db.session.query(Target).all():
            name = target.name
            if not name in props.keys():
                props[name] = {'smiles': target.smiles,
                               'T': [],
                               'dens': OrderedDict([('expt', [])]),
                               'hvap': OrderedDict([('expt', [])])
                               }
            props[name]['T'].append(target.T)
            props[name]['dens']['expt'].append(target.density)
            props[name]['hvap']['expt'].append(target.hvap)

            for i in iterations:
                if i not in props[name]['dens'].keys():
                    props[name]['dens'][i] = []
                    props[name]['hvap'][i] = []

                density, hvap = target.get_npt_result(i)

                props[name]['dens'][i].append(density)
                props[name]['hvap'][i].append(hvap)

        os.chdir(self.CWD)
        pylab.rcParams.update({'font.size': 12})
        for name, prop in props.items():
            pylab.figure(figsize=(6, 8))
            pylab.subplot(211)
            for i, points in prop['dens'].items():
                if i == 'expt':
                    marker = '--'
                elif i == 0:
                    marker = 'x'
                else:
                    marker = 'o'
                pylab.plot(prop['T'], points, marker, label=i)
            y_mean = np.mean(prop['dens']['expt'])
            pylab.ylim(y_mean - 0.2, y_mean + 0.2)
            pylab.legend()
            pylab.title('Density %s %s (g/mL)' % (name, prop['smiles']))

            pylab.subplot(212)
            for i, points in prop['hvap'].items():
                if i == 'expt':
                    marker = '--'
                elif i == 0:
                    marker = 'x'
                else:
                    marker = 'o'
                pylab.plot(prop['T'], points, marker, label=i)
            y_mean = np.mean(prop['hvap'])
            pylab.ylim(y_mean - 20, y_mean + 20)
            pylab.legend()
            pylab.title('HVap %s %s (kJ/mol)' % (name, prop['smiles']))
            pylab.savefig('property-%s.png' % name)
