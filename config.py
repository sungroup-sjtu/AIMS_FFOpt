class Config:
    MS_TOOLS_DIR = '/path/of/ms-tools'
    PACKMOL_BIN = '/path/of/packmol'
    DFF_ROOT = '/path/of/dff_root'
    DFF_TABLE = 'MGI'
    GMX_BIN = 'gmx_serial'

    JOB_MANAGER = 'slurm'
    PBS_ENV_CMD = 'module purge; module load icc gromacs/2018.6'
    # PBS_QUEUE = ('gtx', 32, 2, 16) # For GROMACS 2018, -multidir is weird. Make sure you have even lines in expt data file.
    PBS_QUEUE = ('cpu', 8, 0, 8)
    GMX_MDRUN = 'gmx_gpu mdrun'

    # PBS_QUEUE = ('fast', 12, 0, 12)
    # GMX_MDRUN = 'gmx_fast mdrun'

    # TODO If you are trying to use LJ96 function form, set LJ96 = True
    # Note that LJ96 does not support GPU or OpenMP. So it will be very slow
    # Note that LJ96 and torsion optimization are not compatible. Use at your own risk
    LJ96 = False
