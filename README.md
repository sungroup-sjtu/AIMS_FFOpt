# ffoptimizer
Automatically optimize the LJ parameters for TEAM Force Field against experimental data of density and enthalpy of vaporization.  
`ms-tools`, `DFF`, `Packmol` and `GROMACS` are required.

## Steps

### 1. Configure environment in `config.py`
* Specify paths for required packages.
  The version of `mstools` should be `0.1`.
```
  MS_TOOLS_DIR = '/home/gongzheng/GitHub/ms-tools'
  PACKMOL_BIN = '/share/apps/tools/packmol'
  DFF_ROOT = '/home/gongzheng/apps/DFF/Developing'
```
* Specify the DFF table for assigning atom types
```
  DFF_TABLE = 'MGI'
```
* Specify the `Slurm` partitions and corresponding `GROMACS` binaries.
  The default in config.py is good for `gtx` partition.

### 2. Initialize optimization
* Prepare data file which contains SMILES and experimental density and Hvap data.
  An example is provided as `example_LJ/data.txt`.
  For better GPU performance, make sure that there are **even** lines in expt data file.
* Prepare PPF file which contains initial parameters.
  The parameters to be optimized should be **unfrozen**.
  An example is provided as `example_LJ/initial.ppf`.
* Prepare an empty directory `WORKDIR` for running simulation
* Init
```
  cd example_LJ
  ./run.py init task_name data.txt initial.ppf WORKDIR
```

### 3. Start optimization
```
  ./run.py optimize task_name
```

## Examples
Three examples are provided for optimized LJ, temperature-dependent LJ and dihedral parameters

### 1. Optimize LJ-12-6 parameters: `example_LJ`
Three files are required
* `data.txt` lists the molecules and experimental data and their weight for optimization.
  The names of molecules are arbitrary but should only contain alphabets and numbers.
* `initial.ppf` is the initial parameters exported from DFF.
  The `N12_6` lines for `c_4h3`, `c_4h2` and `h_1` are unfrozen
  so that the epsilon and sigma for these three atom types will be optimized.
* `run.py` is the controlling script. No modification is required for this script.

Start optimization by
```
  mkdir WORKDIR
  ./run.py init LJ data.txt initial.ppf WORKDIR
  ./run.py optimize LJ
```

Check the generated `log` file. Make sure `RSQ` is decreasing and new parameters are reasonable.
It usually converges in less than 6 iterations.

### 2. Optimize temperature-dependent LJ-12-6 parameters: `example_LJ_T`
`data.txt` and `initial.pff` are the same as previous example. Modifications should be made in `run.py` to optimize temperature parameters.
* Two variables `drde_dict` and `drde_atoms` should be specified in `run.py` to optimize the temperature parameter `\lambda` for different atom types.
* `optimizer.drde_dict` lists temperature parameters that **will be fixed** during the optimization.
* `optimizer.optimize(drde_atoms=...)` lists temperature parameters that are **subject to optimization**.
* the `\lambda` parameters are always named as `xxx_dl`. It will match all atom types starting with `xxx`.
  For example, `c_4_dl` is the `\lambda` parameter for `c_4h2`, `c_4h3`, `c_4o`...

Note that
* If a parameter appears in both `drde_dict` and `drde_atoms`, it will be subject to optimization.
* It is better to have at most one free `\lambda` parameter to get rid of coupling between parameters.

Start optimization and check the `log` file. Make sure `RSQ` is decreasing and new parameters are reasonable.
The `RSQ` should be smaller than previous example because of the introduction of temperature-dependence.

### 3. Optimize temperature-dependent LJ-12-6 parameters together with dihedral: `example_LJ_T_dihedral`
`data.txt` is the same as previous example. Modifications should be made in `initial.pff` and `run.py` to optimize dihedrals.
The dihedrals are fitted to QM energy surface by using DFF.
* Unfreeze the dihedral parameters in `initial.ppf` so that they can be optimized.
* Variable `torsions` lists the dihedrals subject to optimization.
* Corresponding `MSD` and `QMD` files should be provided.

Note that
* Dihedral fitting here only make sense for the backbone of linear molecules. For example, alkanes, diphenyl...

Start optimization and check the `log` file. Make sure `RSQ` is decreasing and new parameters are reasonable.
Check the `dft` file to make sure dihedrals are correctly fitted.
