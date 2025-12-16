# dislocate

Python tools for analyzing dislocation core structures in crystalline materials and run hightrouput atomistic computations.

## Overview

The `dislocate` package provides a comprehensive toolkit for:
    - `core_analysis`: analyzing dislocation core structures using atomistic simulations
    - `core_search`: finding (meta)-stable dislocation core structure via the Elastic Singularity Placement Method and NEB Core Search method
    - `utils`: high-throughput computations of lattice parameters, elastic constants of HCP, BCC, and FCC crystal and stacking fault energies of an HCP crsytal and other utility functions.

It integrates with LAMMPS, OVITO, and Babel to perform automated dislocation analysis workflows, including core structure analysis (Nye tensor, elastic stability parameter, core angle, core eccentricty, core surface), and energy calculations.

## Installation

### Install
`pip install git+https://github.com/dj-shark31/dislocate.git`

### Update package
`pip install --upgrade git+https://github.com/dj-shark31/dislocate.git`

### Prerequisites

1. **Python 3.8+** with scientific Python stack
2. **Babel** for dislocation calculations:
   - Download and build from: http://emmanuel.clouet.free.fr/Programs/Babel
   - Configure paths in `config.yaml`

3. **OVITO** (Modified Version):
   - Modified OVITO by Eric Rothchild: https://github.com/Erothch/ovitoWallace
   - Docker container: `dj31/ovito-wallace:amd64`
   - Singularity: `singularity pull ovitowallace.sif docker://dj31/ovito-wallace:amd64`

## Project Structure

```
dislocate/
├── core_analysis/          # Main dislocation analysis workflow
│   ├── analyze_core.py     # Core analysis script
│   ├── get_data.py         # Orchestrate core analysis for one core
│   ├── assemble.py         # Compile results into one output file
│   ├── visualize_core.py   # Script to plot the dislocation core maps
│   ├── Babel/             # Nye tensor analysis
│   ├── EnergyStress/      # Energy and stress calculations
│   ├── OvitoFitting/      # OVITO integration (elastic stability paramater, DXA, CNA, core angle, eccentricity, area)
│   └── Pattern/           # Pattern detection. Required for the Nye tensor analysis
├── core_search/            # Dislocation core search algorithms
│   ├── ESPM.py            # Elastic Singularity Positioning Method
│   |── neb_core_search.py # NEB core search method
|   └── highthroughput_ESPM.py # Launch ESPM comoputation on supercomupter
├── utils/                  # Utility functions
│   ├── atomistic_tools.py  # Atomistic simulation tools
│   ├── translate_dislocation.py # Dislocation dipole translation
│   ├── peierls_barrier_neb.py # Peierls barrier calculation
│   └── config_loader.py    # Functions to retreive path from config.yaml for external tools
├── get_properties_scripts/ # Material property calculations
│   ├── get_elastic_constants.py
│   ├── get_lattice_parameters.py
│   └── get_stacking_fault_energy.py
├── examples/               # Example workflows
├── potentials/             # Interatomic potentials
└── sqs_cells/             # Special quasirandom structures
```

## Dependencies

### Core Dependencies
- **NumPy** (≥1.26.0): Numerical computing
- **SciPy** (≥1.13.0): Scientific computing
- **Matplotlib** (≥3.9.0): Plotting and visualization
- **Pandas** (≥2.3.0): Data manipulation
- **ASE** (≥3.23.0): Atomic simulation environment
- **PyMatGen** (≥2024.8.0): Materials analysis
- **PyYAML** (≥6.0): Configuration files

### Machine Learning
- **MACE-Torch** (≥0.3.0): Machine learning potentials
- **PyTorch** (≥2.2.0): Deep learning framework

### External Tools
- **LAMMPS** (≥2022.0): Molecular dynamics
- **OVITO**: Structure analysis (modified version)
- **Babel**: Dislocation calculations

## Configuration

Create a `config.yaml` file to specify paths to external tools:

```yaml
tools:
  ovitosif: /path/to/ovito.sif
  babel: /opt/babel/bin/babel
  displacement: /opt/babel/bin/displacement
  patternDetect: /opt/babel/bin/patternDetect
  patternInit: /opt/babel/bin/patternInit
```

## Usage Examples

See `*/README.md`

## Documentation

- **Core Analysis**: See `dislocate/core_analysis/README.md`
- **Core Search**: See `dislocate/core_search/README.md`
- **Utils**: See `dislocate/utils/README.md`
- **Environment Setup**: See `ENVIRONMENT_SETUP.md`

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{dislocate,
  title={dislocate: Python tools for analyzing dislocation core structures and hightroughput computations},
  author={Jany, David},
  year={2024},
  url={https://github.com/dj-shark31/dislocate}
}
```

## Acknowledgement

This repository was completed thanks to the funding provided by the U.S. Office of Naval Research under Grants No. N00014-17-1- 2283, No. N00014-19-1-2376. and the National Science Foundation under Grant No. 2324022.

## Contact

- **Author**: David Jany
- **Email**: djany31@gmail.com
- **GitHub**: https://github.com/dj-shark31/dislocate
