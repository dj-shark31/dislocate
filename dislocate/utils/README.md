# Utils

This directory contains utility functions and tools for atomistic simulations and dislocation analysis.
## Overview

The utils module provides essential tools for:
- Atomic structure manipulation
- Dislocation translation and analysis
- Peierls barrier calculations
- Path retrieving for external tools
- Job submission utilities

## Files

### `atomistic_tools.py`
Core tools for atomistic simulations and structure manipulation.

**Key Functions:**
- `build_atoms()`: Create crystal structures (HCP, BCC, FCC)
- `alloy_sqs()`: alloy an SQS cells
- `alloy_randomly()`: randomly alloy an arbitrary cell
- `setup_calculator()`: Configure interatomic potentials
- `get_energy()` and `get_stress()`: Compute energy and stress
- `load_sqs()`: Generate an ASE Atom from a SQS structure in `sqs_cells`
- `get_lattice_parameters()`: compute lattice parameters of a cell
- `get_cij()`: compute elastic constants of a cell
- `get_lattice_parameters()`: compute the basal or prism stacking fault of an HCP cell.

### `translate_dislocation.py`
Tools for translating dislocation structures.

**Key Functions:**
- `translate_dislocation()`: Translate dislocation by specified vector
- `reorder_atoms()`: Reorder atoms for consistency to run a NEB or an dislocation analysis
- `main()`: outputs the translated structure in a POSCAR file

**Usage:**
```bash
# Translate dislocation by [1, 0, 0]
python -m dislocate.utils.translate_dislocation --input_cell core.poscar --output_file translated_core.poscar --translation 1 0 0 --n_cells 8 8 1
```

### `peierls_barrier_neb.py`
Nudged Elastic Band calculations for Peierls barriers.

**Key Functions:**
- `run_neb()`: Calculate Peierls barrier using NEB
- `relax_intermediate_images()`: Relax the intermediate images of the NEB
- `main()`: Computes Peierls barrier and can relax intermediate images

**Usage:**
```bash
# Calculate Peierls barrier
python -m dislocate.utils.peierls_barrier_neb --initial core1.poscar --final core2.poscar --potential_path /path/to/mace_model.pt --n_images 5 --output_dir neb_results --neb_fmax 0.005 --relax_fmax 0.001 --relax_intermediate_images true --output_preopt_images false --perform_neb true --device cpu
```

### `config_loader.py`
Functions to retrieve path for external tools

**Key Functions:**
- `load_config()`: Load YAML configuration files
- `get_tool_paths()`: Extract tool paths from config

**Usage:**
```python
from dislocate.utils.config_loader import get_tool_path

# Load configuration
babel_path = get_tool_path('babel')
```

### `convert_sqs.py`
Special Quasirandom Structure (SQS) conversion tools.

**Usage:**
```bash
# Convert SQS cell to POSCAR file
python -m dislocate.utils.convert_sqs --in_file bestsqs.out --out_file core.poscar --a 2.9 --c 4.7
```

### `highthroughput.py`
Job submission utilities for computing clusters.

**Key Functions:**
- `submit_savio_job()`: Submit jobs to queue systems
- `high_throughput_property()`: Submit hightrouhput computations of lattice parameters, elastic constants, or HCP stacking faults
- `find_compositions()`: Finds all compositions given a list of elements and a delta of concentration

## Dependencies

- **ASE**: Atomic simulation environment
- **NumPy/SciPy**: Numerical computations
- **PyYAML**: Configuration file parsing
- **MACE-Torch**: Machine learning potentials
- **LAMMPS**: Molecular dynamics (optional)

## Integration

These utilities are used throughout the dislocate package:
- **Core Analysis**: Structure manipulation and energy calculations
- **Core Search**: Configuration management and job submission
- **Property Scripts**: Structure setup and analysis
