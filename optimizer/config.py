import socket
from collections import OrderedDict


class MacConfig():
    DFF_ROOT = '/Users/zheng/Projects/DFF/Developing'
    PACKMOL_BIN = '/Users/zheng/Projects/DFF/Developing/bin32m/Packmol/packmol.exe'
    GMX_BIN = '/opt/gromacs/2016.3/bin/gmx'

    MS_TOOLS_DIR = '/Users/zheng/Projects/msd-server'
    WORK_DIR = '/tmp/MSDServer/'

    JOB_MANAGER = 'local'
    NPROC_PER_JOB = 2


class ClusterConfig():
    DFF_ROOT = '/share/apps/dff/msdserver'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    GMX_BIN = '/share/apps/gromacs/2016.3-static-compatible/bin/gmx'

    MS_TOOLS_DIR = '/home/msdserver/msd-server'
    WORK_DIR = '/share/workspace/msdserver/MSDServer/'

    JOB_MANAGER = 'torque'
    QUEUE_DICT = OrderedDict([('cpu', 8), ('fast', 12)])


class PiConfig():
    DFF_ROOT = '/lustre/home/nishsun/gongzheng/apps/dff'
    PACKMOL_BIN = '/lustre/home/nishsun/software/tools/packmol'
    GMX_BIN = '/lustre/home/nishsun/software/gromacs/5.1.4-static/bin/gmx'

    MS_TOOLS_DIR = '/lustre/home/nishsun/gongzheng/workspace/msd-server'
    WORK_DIR = '/lustre/home/nishsun/gongzheng/workspace/MSDServer'

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
