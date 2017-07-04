from collections import OrderedDict
from typing import Dict


class PPF():
    def __init__(self, ppf_file=None, string=None):
        lines = []
        if ppf_file is not None:
            with open(ppf_file) as f:
                lines = f.read().splitlines()
        if string is not None:
            lines = string.splitlines()
        self.terms = [l.strip() for l in lines]

    def __str__(self):
        return '\n'.join(self.terms)

    def get_adj_nb_paras(self) -> OrderedDict:
        terms = OrderedDict()
        for term in self.terms:
            if not (term.startswith('N12_6') or term.startswith('BINC')):
                continue

            words = term.split(':')
            words = [w.strip() for w in words]
            if term.startswith('N12_6'):
                a_type = words[1]
                paras = words[2]

                words = paras.split(',')
                words = [w.strip() for w in words]
                r0 = words[0]
                e0 = words[1]
                if not r0.endswith('*'):
                    terms['%s_r0' % a_type] = float(r0)
                if not e0.endswith('*'):
                    terms['%s_e0' % a_type] = float(e0)
            elif term.startswith('BINC'):
                a_types = words[1]
                para = words[2]

                words = a_types.split(',')
                words = [w.strip() for w in words]
                a1_type = words[0]
                a2_type = words[1]
                if not para.endswith('*'):
                    terms['%s_%s_bi' % (a1_type, a2_type)] = float(para)
        return terms

    def set_nb_paras(self, new_paras: Dict):
        terms = [''] * len(self.terms)
        replaced = False
        for i, term in enumerate(self.terms):
            if not (term.startswith('N12_6') or term.startswith('BINC')):
                terms[i] = term
                continue

            if term.startswith('N12_6'):
                words = term.split(':')
                words = [w.strip() for w in words]
                a_type = words[1]
                paras = words[2]
                words = paras.split(',')
                words = [w.strip() for w in words]
                r0 = words[0]
                e0 = words[1]

                r0_key = '%s_r0' % a_type
                if r0_key in new_paras.keys():
                    r0 = new_paras[r0_key]
                    replaced = True

                e0_key = '%s_e0' % a_type
                if e0_key in new_paras.keys():
                    e0 = new_paras[e0_key]
                    replaced = True

                new_term = 'N12_6: %s: %s, %s:' % (a_type, str(r0), str(e0))
                terms[i] = new_term

            elif term.startswith('BINC'):
                words = term.split(':')
                words = [w.strip() for w in words]
                a_types = words[1]
                bi = words[2]

                words = a_types.split(',')
                words = [w.strip() for w in words]
                a1_type = words[0]
                a2_type = words[1]

                bi_key = '%s_%s_bi' % (a1_type, a2_type)
                if bi_key in new_paras.keys():
                    bi = new_paras[bi_key]
                    replaced = True

                new_term = 'BINC: %s, %s: %s:' % (a1_type, a2_type, str(bi))
                terms[i] = new_term

        self.terms = terms
        return replaced

    def write(self, ppf_out):
        with open(ppf_out, 'w') as f:
            f.write(str(self))

    def freeze_torsions(self):
        terms = [''] * len(self.terms)
        for i, term in enumerate(self.terms):
            if not (term.startswith('TCOSP')):
                terms[i] = term
                continue
            words = term.split(':')
            words = [w.strip() for w in words]
            a_types = words[1]
            paras = words[2]
            para_words = paras.split(',')
            para_words = [w.strip() for w in para_words]
            for k in range(len(para_words)):
                para_word = para_words[k]
                if not para_word.endswith('*'):
                    para_words[k] += '*'
            new_paras = ', '.join(para_words)
            new_term = 'TCOSP: %s: %s:' % (a_types, new_paras)
            terms[i] = new_term
        self.terms = terms

    def relax_torsion(self, torsion_key):
        terms = [''] * len(self.terms)
        for i, term in enumerate(self.terms):
            if not (term.startswith('TCOSP')):
                terms[i] = term
                continue
            words = term.split(':')
            words = [w.strip() for w in words]
            a_types = words[1]

            a_type_words = a_types.split(',')
            a_type_words = [w.strip() for w in a_type_words]

            key_words = torsion_key.split(',')
            key_words = [w.strip() for w in key_words]

            if a_type_words == key_words or a_type_words == list(reversed(key_words)):
                new_paras = '0*, 0.0, 3*, 0*, 0.0, 1*, 180*, 0.0, 2*'
                new_term = 'TCOSP: %s: %s:' % (a_types, new_paras)
                terms[i] = new_term
            else:
                terms[i] = term
        self.terms = terms

    def fit_torsion(self, qmd=None, msd=None, restraint=None, torsion_key=None):
        import os
        import random
        from mstools.wrapper import DFF
        from .config import Config

        ### relax only one torsion
        self.freeze_torsions()
        self.relax_torsion(torsion_key)

        ### Backup adj_nb_paras
        adj_nb_paras = self.get_adj_nb_paras()

        dff = DFF(dff_root=Config.DFF_ROOT)
        ppf_tmp = 'tmp-%i.ppf' % random.randint(1E5, 1E6)
        ppf_tmp2 = 'tmp-%i.ppf' % random.randint(1E5, 1E6)
        self.write(ppf_tmp)
        try:
            dff.fit_torsion(qmd, msd, ppf_tmp, ppf_tmp2, restraint)
        except:
            raise
        else:
            self.__init__(ppf_tmp2)
            ### restore adj_nb_paras
            self.set_nb_paras(adj_nb_paras)
            os.remove(ppf_tmp)
            os.remove(ppf_tmp2)

    @staticmethod
    def get_delta_for_para(key):
        if key.endswith('r0'):
            delta = 0.05
        elif key.endswith('e0'):
            delta = 0.001
        elif key.endswith('bi'):
            delta = 0.005
        else:
            raise Exception('Unknown parameter: ' + key)

        return delta

    @staticmethod
    def get_bound_for_para(key):
        if key.endswith('r0'):
            bound = (3, 5)
            if key.startswith('h'):
                bound = (2, 3)
        elif key.endswith('e0'):
            bound = (0.01, 0.5)
            if key.startswith('h'):
                bound = (0.01, 0.1)
            if key.startswith('c'):
                bound = (0.01, 0.1)
        elif key.endswith('bi'):
            bound = (-0.6, 0.6)
        else:
            raise Exception('Unknown parameter: ' + key)

        return bound
