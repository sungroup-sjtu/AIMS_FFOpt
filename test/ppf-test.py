from mstools.wrapper.ppf import PPF

ppf = PPF('TEAM_LS.ppf')
ppf.write('new.ppf')
print(ppf.get_adj_nb_paras())