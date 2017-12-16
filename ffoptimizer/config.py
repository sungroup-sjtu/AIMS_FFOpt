import socket
from collections import OrderedDict


class BaseConfig():
    DFF_TABLE = 'TEAM_IL'


class MacConfig(BaseConfig):
    PACKMOL_BIN = '/Users/zheng/Projects/DFF/Developing/bin64m/Packmol/packmol.exe'
    GMX_BIN = '/opt/gromacs/2016.3/bin/gmx'
    DFF_ROOT = '/Users/zheng/Projects/DFF/Developing'

    MS_TOOLS_DIR = '/Users/zheng/Projects/msd-server'

    JOB_MANAGER = 'local'
    QUEUE_DICT = OrderedDict([(None, 2)])


class ClusterConfig(BaseConfig):
    PACKMOL_BIN = '/share/apps/tools/packmol'
    GMX_BIN = '/share/apps/gromacs/2016.3/bin/gmx_gpu'
    # GMX_BIN = '/share/apps/gromacs/2016.3/bin/gmx_fast'
    DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'

    MS_TOOLS_DIR = '/home/gongzheng/GitHub/ms-tools'

    JOB_MANAGER = 'torque'
    QUEUE_DICT = OrderedDict([('gtx', 2)])
    # QUEUE_DICT = OrderedDict([('cpu', 8)])
    # QUEUE_DICT = OrderedDict([('cpu', 8), ('fast', 12)])


class PiConfig(BaseConfig):
    PACKMOL_BIN = '/lustre/home/nishsun/software/tools/packmol'
    GMX_BIN = '/lustre/home/nishsun/software/gromacs/5.1.4-static/bin/gmx'
    DFF_ROOT = '/lustre/home/nishsun/gongzheng/apps/dff'

    MS_TOOLS_DIR = '/lustre/home/nishsun/gongzheng/workspace/msd-server'

    JOB_MANAGER = 'slurm'
    QUEUE_DICT = OrderedDict([('cpu', 8)])


if socket.gethostname() == 'cluster.hpc.org':
    Config = ClusterConfig
elif socket.gethostname() == 'mu06.pi.sjtu.edu.cn':
    Config = PiConfig
elif socket.gethostname() == 'z-Mac.local':
    Config = MacConfig
else:
    raise Exception('MSDServer will not work on this machine')
