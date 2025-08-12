# Core Search

This directory contains algorithms for searching and identifying dislocation core structures in crystalline materials.

## Overview

The core search module provides two main approaches for finding dislocation core configurations:

1. **Elastic Singularity Positioning Method (ESPM)**: A systematic approach to explore dislocation core configurations
2. **NEB Core Search Method**: Uses the Nudged Elastic Band method to find meta-stable core structure in the transition paths between stable core configurations.

## Files

### `ESPM.py`
Implements the Elastic Singularity Positioning Method for dislocation core search. Builds dislocation dipole from a reference cell using babel and then relaxes the cell.

**Usage:**
```bash
python -m dislocate.core_search.ESPM --reference_cell POSCAR --output_dir ./output --n_cells 8 8 --potential_path /path/to/mace_model.pt --xpos 0 0 10 -10 5 -5 --ypos 5 -5 0 0 5 -5 --babel_path /path/to/babel --fmax 0.0005 --meshing 20 30 --remove_original true --cij 177.1 84.8 82.9 193.8 54.8
"""
```

### `neb_core_search.py`
Implements NEB-based core search algorithms. Performs NEB between 2 core configurations. Relax intermediate images. Identify new core configurations. If new core configurations exist, we repeat the search.

**Usage:**
```bash
# Run NEB search between two core configurations
python -m dislocate.core_search.neb_core_search --initial_file initial.poscar --final_file final.poscar --potential_path /path/to/mace_model.pt --n_images 5 --output_dir neb_core_search --max_iterations 10 --energy_tolerance 0.001 --neb_fmax 0.005 --relax_fmax 0.001 --device cpu
```

### `highthroughput_ESPM.py`
Script for submitting ESPM calculations to computing clusters.

**Usage:**
```bash
# Run NEB search between two core configurations
python -m dislocate.core_search.submitruns_ESPM \
     --elements Hf Ti Zr --delta_x 0.1 --potential_path /path/to/mace_model.pt --reference_cell POSCAR \
    --n_cells 8 8 --xpos 0 0 10 -10 5 -5 --ypos 5 -5 0 0 5 -5 \
    --babel_path /path/to/babel --fmax 0.0005 --meshing 20 30 \
    --remove_binaries false --remove_unaries true \
    --cij 177.1 84.8 82.9 193.8 54.8 \
    --job_name ESPM --partition savio4_htc --cpus 10 --time 5:00:00 --qos savio_lowprio --account co_chrzangroup \
    --analyze_core true
```

### `input.babel_S`
Input file template to build dislocation dipole with Babel in ESPM.

## Dependencies

- **Babel**: Required for dislocation calculations
- **ASE**: Atomic simulation environment
- **NumPy/SciPy**: Numerical computations
- **PyYAML**: Configuration file parsing
