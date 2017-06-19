import os
import time
import shutil
import json
from collections import OrderedDict

import numpy as np
from lmfit import Parameters, Minimizer

from .ppf import PPF
from .db import DB
from .models import Target, Result

from sqlalchemy import and_


class Optimizer():
    def __init__(self, db_file, cwd=os.getcwd()):
        self.db = DB(db_file)
        self.db.conn()
        self.CWD = cwd

    @property
    def iteration(self):
        target = self.db.session.query(Target).first()
        return target.iteration

    def optimize(self, ppf_file, wExpansivity, qmd=None, msd=None, torsion=None):
        LOG = os.path.join(self.CWD, 'Opt-%s.log' % os.path.basename(ppf_file)[:-4])

        def residual(params: Parameters):
            ### if result exist in database, ignore calculation
            result = self.db.session.query(Result).filter(
                and_(
                    Result.ppf == ppf_file,
                    Result.parameter == str(params)
                )).first()
            if result is not None:
                R = result.residual
                if R is not None:
                    return json.loads(R)
            ###

            ### save ppf file and run NPT
            ppf = PPF(ppf_file)
            paras = OrderedDict()
            for k, v in params.items():
                paras[k] = v.value
            ppf.set_nb_paras(paras)

            if torsion is not None:
                print('Fit torsion based on new non-bonded parameters...')
                ppf.fit_torsion(qmd, msd, torsion)

            shutil.copy(ppf_file, ppf_file + '.bak-%i' % self.iteration)
            ppf.write(ppf_file)

            for target in self.db.session.query(Target).all():
                if not target.npt_started(ppf_file):
                    ### next iteration
                    target.iteration += 1
                    target.run_npt(ppf_file, paras)
            self.db.session.commit()
            os.chdir(self.CWD)
            ###

            while True:
                FINISHED = True
                for target in self.db.session.query(Target).all():
                    if not target.npt_finished(ppf_file):
                        FINISHED = False
                        break

                if FINISHED:
                    break
                else:
                    current_time = time.strftime('%m-%d %H:%M')
                    print(current_time + ' Job still running. Wait ...')
                    time.sleep(60)

            R_dens = []
            R_hvap = []
            R_expa = []
            targets = self.db.session.query(Target).all()
            for target in targets:
                dens, hvap = target.get_npt_result(os.path.basename(ppf_file)[:-4])
                R_dens.append((dens - target.density) / target.density * 100 * target.wDensity)  # deviation  percent
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
            if result is None:
                result = Result(ppf=ppf_file, parameter=str(params))

            result.residual = json.dumps(R)
            self.db.session.add(result)
            self.db.session.commit()
            ###

            ### write current parameters and residual to log
            txt = '\nITERATION: %i\n' % self.iteration
            txt += '\nPARAMETERS:\n'
            for k, v in params.items():
                txt += '%10.5f  %s\n' % (v.value, k)
            txt += '\nRESIDUE:\n'
            for i, r in enumerate(R_dens):
                target = targets[i]
                prop = 'density'
                weight = target.wDensity
                txt += '%8.2f  %10s %10s %.2f %8.2f\n' % (r, prop, target.name, weight, r / weight)
            for i, r in enumerate(R_hvap):
                target = targets[i]
                prop = 'hvap'
                weight = target.wHvap
                txt += '%8.2f  %10s %10s %.2f %8.2f\n' % (r, prop, target.name, weight, r / weight)
            for i, r in enumerate(R_expa):
                target = targets[i * 2]
                prop = 'expan'
                weight = wExpansivity
                txt += '%8.2f  %10s %10s %.2f %8.2f\n' % (r, prop, target.name, weight, r / weight)

            txt += '\nRSQ: %.2f\n' % np.sum(list(map(lambda x: x ** 2, R)))

            print(txt)
            with open(LOG, 'a') as log:
                log.write(txt)
            ###

            return R

        def jacobian(params: Parameters):
            ### if result exist in database, ignore calculation
            result = self.db.session.query(Result).filter(
                and_(
                    Result.ppf == ppf_file,
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
            targets = self.db.session.query(Target).all()
            for target in targets:
                dDdp_list, dHdp_list = target.get_dDens_dHvap_list_from_paras(ppf_file, paras)
                J_dens.append([i / target.density * 100 * target.wDensity for i in dDdp_list])  # deviation  percent
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
            if result is None:
                result = Result(ppf=ppf_file, parameter=str(params))

            result.jacobian = json.dumps(J)
            self.db.session.add(result)
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
            ### clear _finished_ and job.sh for next iteration
            ### TODO this is not good
            for target in self.db.session.query(Target).all():
                target.clear_npt_result(ppf_file)

        ppf = PPF(ppf_file)
        params = Parameters()
        for k, v in ppf.get_adj_nb_paras().items():
            bound = PPF.get_bound_for_para(k)
            params.add(k, value=v, min=bound[0], max=bound[1])

        minimize = Minimizer(residual, params, iter_cb=callback)
        result = minimize.leastsq(Dfun=jacobian, ftol=0.0001)
        print(result.lmdif_message, '\n')

        return result.params

    def npt(self, ppf_file):
        for target in self.db.session.query(Target).all():
            target.run_npt(ppf_file)

    def plot(self, ppfs, iteration=None):
        try:
            import pylab
        except:
            print('matplotlib not found. will not plot.')
            return

        ppfs = [ppf[:-4] for ppf in ppfs]
        props = dict()
        for target in self.db.session.query(Target).all():
            if not target.name in props.keys():
                props[target.name] = {'smiles': target.smiles, 'T': [], 'd_exp': [], 'h_exp': [],
                                      'd_sim': OrderedDict(), 'h_sim': OrderedDict()}
            props[target.name]['T'].append(target.T)
            props[target.name]['d_exp'].append(target.density)
            props[target.name]['h_exp'].append(target.hvap)

            for ppf in ppfs:
                if ppf not in props[target.name]['d_sim'].keys():
                    props[target.name]['d_sim'][ppf] = []
                    props[target.name]['h_sim'][ppf] = []

                density, hvap = target.get_npt_result(ppf, iteration)

                props[target.name]['d_sim'][ppf].append(density)
                props[target.name]['h_sim'][ppf].append(hvap)

        os.chdir(self.CWD)
        pylab.rcParams.update({'font.size': 12})
        for tid, prop in props.items():
            pylab.figure(figsize=(6, 8))
            pylab.subplot(211)
            pylab.plot(prop['T'], prop['d_exp'], '--')
            for ppf, points in prop['d_sim'].items():
                marker = 'x' if ppf.endswith('init') else 'o'
                pylab.plot(prop['T'], points, marker, label=ppf)
            y_mean = np.mean(prop['d_exp'])
            pylab.ylim(y_mean - 0.2, y_mean + 0.2)
            pylab.legend()
            pylab.title('Density %s %s (g/mL)' % (tid, prop['smiles']))

            pylab.subplot(212)
            pylab.plot(prop['T'], prop['h_exp'], '--')
            for ppf, points in prop['h_sim'].items():
                marker = 'x' if ppf.endswith('init') else 'o'
                pylab.plot(prop['T'], points, marker, label=ppf)
            y_mean = np.mean(prop['h_exp'])
            pylab.ylim(y_mean - 20, y_mean + 20)
            pylab.legend()
            pylab.title('HVap %s %s (kJ/mol)' % (tid, prop['smiles']))
            pylab.savefig('property-%s.png' % tid)

    def init_db(self, filename):
        with open(filename) as f:
            lines = f.read().splitlines()
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue
            words = line.split()

            target = Target()
            target.name = words[0]
            target.smiles = words[1]
            target.T = int(words[2])
            target.P = int(words[3])

            if self.db.session.query(Target).filter(
                    and_(Target.smiles == target.smiles,
                         Target.T == target.T,
                         Target.P == target.P
                         )).count() > 0:
                continue

            target.density = float(words[4])
            target.wDensity = float(words[5])
            target.hvap = float(words[6])
            target.wHvap = float(words[7])
            target.calc_n_mol()
            self.db.session.add(target)

        try:
            self.db.session.commit()
        except:
            self.db.session.rollback()
