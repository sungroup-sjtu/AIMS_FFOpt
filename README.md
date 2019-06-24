# ffoptimizer
Automatically optimize the LJ parameters for TEAM Force Field against density and Hvap data.  
`ms-tools`, `DFF`, `Packmol` and `GROMACS` are required for this code.

## Steps

### 1. Specify paths for required packages in `config.py`  
For example
```
 MS_TOOLS_DIR = '/home/gongzheng/GitHub/ms-tools'
 PACKMOL_BIN = '/share/apps/tools/packmol'
 DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
 ...
```

### 2. Initialize optimization
* Prepare data file which contains SMILES and experimental density and Hvap data. An example is provided as `run-example/data.txt`
* Prepare ppf file which contains ff parameters. The parameters to be optimized should be **unfrozen**. An example is provided as `run-example/TEAM_LS.ppf`
* Prepare an empty directory `WORKDIR` for running simulation
* Init
```
  cd run-example
  ./run.py init task_name data.txt TEAM_LS.ppf WORKDIR
```

### 3. Start optimization
```
  ./run.py optimize task_name
```
