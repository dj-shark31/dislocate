# dislocate
Python tools for analyzing dislocation core structures.

## Dependencies

### Babel
This package requires Babel for dislocation calculations. You will need to:
1. Download and build Babel from source: http://emmanuel.clouet.free.fr/Programs/Babel
2. Provide path for executable in core_analysis/get_data.py and core_search/ESPM.py 

### OVITO (Modified Version)
A modified version of OVITO is required for visualization:
- The custom OVITO Singularity container (`ovito.sif`) is needed
- This container will be made available for download (location TBD)
- Once downloaded, place `ovito.sif` in the `bin/` directory

Note: The standard OVITO installation will not work - you must use the modified version.
