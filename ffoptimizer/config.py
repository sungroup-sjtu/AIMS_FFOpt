class Config:
    MS_TOOLS_DIR = '/home/gongzheng/GitHub/ms-tools'
    PACKMOL_BIN = '/share/apps/tools/packmol'
    DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
    DFF_TABLE = 'MGI'
    GMX_BIN = 'gmx_serial'

    JOB_MANAGER = 'slurm'
    PBS_ENV_CMD = '''
module purge
module load gromacs/2018.6
'''
    PBS_QUEUE = ('gtx', 32, 2, 16)
    # PBS_QUEUE= ('cpu', 8, 0, 8)
    GMX_MDRUN = 'gmx_gpu mdrun'

    # PBS_QUEUE= ('fast', 12, 0, 12)
    # GMX_MDRUN = 'gmx_fast mdrun'
