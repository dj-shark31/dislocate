# Material Property Scripts

This directory contains scripts for calculating fundamental material properties using atomistic simulations.

## Overview

The material property scripts provide automated workflows for computing:
- Elastic constants
- Lattice parameters
- Stacking fault energies

## Files

### `get_elastic_constants.py`
Calculates elastic constants from a reference cell and outputs results in a file. The reference cell is strained in multiple directions. Each strained image is relaxed and the stress tensor is computed. The elastic constant tensor is then computed form the strain and stress. 

**Output:**
- Elastic constants matrix (Voigt notation)

### `get_lattice_parameters.py`
Determines equilibrium lattice parameters given a list of elements and composition. SQS cells or arbitrary cells can be used.

**Output:**
- Equilibrium lattice constants (a and c)

### `get_stacking_fault_energy.py`
Calculates basal or prism stacking fault energies for an HCP crystal.

**Output:**
- Stacking fault energy (eV/Å²)

## Dependencies

- **ASE**: Atomic simulation environment
- **NumPy/SciPy**: Numerical computations and fitting
- **MACE-Torch**: Machine learning potentials
- **PyYAML**: Configuration management

## Output Formats

### Elastic Constants
c11 c11_std c33 c33_std c12 c12_std c13 c13_std c44 c44_std
180 0.2 190 0.5 84 0.2 74 0.6 55 0.1

### Lattice Parameters
a a_std c c_std
2.9 0.01 4.6 0.02

### Stacking Fault Energy
sf_energies sf_energies_std
50 0.4

## Tips

1. **Convergence**: Ensure energy and force convergence for accurate results
2. **System size**: Use appropriate supercell sizes for property calculations