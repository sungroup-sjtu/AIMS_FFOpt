import socket


class BaseConfig():
    GMX_MDRUN = None
    DFF_TABLE = 'MGI'


class ClusterConfig(BaseConfig):
    MS_TOOLS_DIR = '/home/gongzheng/GitHub/ms-tools'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'

    JOB_MANAGER = 'slurm'
    PBS_ENV_CMD = '''
module purge
module load gcc gromacs/2016.5
'''
    # PBS_QUEUE_LIST = [('gtx', 32, 2, 16)]
    # PBS_QUEUE_LIST = [('cpu', 8, 0, 8)]
    # GMX_BIN = 'gmx_gpu'

    PBS_QUEUE_LIST = [('fast', 24, 0, 12)]
    GMX_BIN = 'gmx_fast'


class PiConfig(BaseConfig):
    pass


if socket.gethostname() == 'cluster.hpc.org':
    Config = ClusterConfig
elif socket.gethostname() == 'mu06.pi.sjtu.edu.cn':
    Config = PiConfig
else:
    raise Exception('ffoptimizer will not work on this machine')
