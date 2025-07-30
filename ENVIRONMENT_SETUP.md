# Dislocate Environment Setup

This guide explains how to set up the conda environment for the dislocate project.

## Quick Setup

### 1. Create the conda environment

```bash
conda create -n dislocate python=3.9 -y
conda activate dislocate
```

### 2. Install core packages

```bash
# Core scientific packages
conda install numpy scipy matplotlib pandas -y

# Atomic simulation environment
conda install -c conda-forge ase -y

# YAML configuration
conda install pyyaml -y

# Materials science
pip install pymatgen

# Machine learning potentials
pip install mace-torch
```

### 3. Verify installation

```bash
python test_environment.py
```

## Detailed Installation

### Core Dependencies

The project requires the following core packages:

- **Python 3.9+**: Base Python environment
- **NumPy**: Numerical computing
- **SciPy**: Scientific computing
- **Matplotlib**: Plotting and visualization
- **Pandas**: Data manipulation
- **ASE**: Atomic Simulation Environment
- **PyYAML**: YAML configuration files
- **PyMatGen**: Materials science toolkit
- **MACE**: Machine learning atomic calculator
- **PyTorch**: Deep learning framework (required by MACE)

### Optional Dependencies

- **LAMMPS**: For LAMMPSlib calculator (if needed)
- **OVITO**: For structure analysis (external tool)

## Environment Management

### Activate the environment

```bash
conda activate dislocate
```

### Deactivate the environment

```bash
conda deactivate
```

### Remove the environment

```bash
conda env remove -n dislocate
```

## Testing the Environment

Run the test script to verify all packages are working:

```bash
python test_environment.py
```

This will test:
- All required package imports
- Project-specific imports
- Basic functionality of key packages

## Troubleshooting

### Import Errors

If you get import errors, make sure:
1. The conda environment is activated: `conda activate dislocate`
2. All packages are installed: Check with `conda list` or `pip list`
3. The project root is in Python path when running scripts

### MACE Installation Issues

If MACE installation fails:
1. Make sure PyTorch is installed first
2. Try installing from conda-forge: `conda install -c conda-forge mace-torch`
3. Check the MACE documentation for specific requirements

### LAMMPS Issues

LAMMPS is optional and only needed for LAMMPSlib calculator:
1. Install via conda: `conda install -c conda-forge lammps`
2. Or install via pip: `pip install lammps`
3. Or skip if not using LAMMPSlib

## Requirements File

A `requirements.txt` file is provided with all necessary packages:

```bash
pip install -r requirements.txt
```

Note: Some packages may need to be installed via conda for better compatibility. 