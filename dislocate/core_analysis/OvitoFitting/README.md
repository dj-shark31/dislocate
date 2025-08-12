# OVITO Integration and Fitting

This directory contains scripts for integrating OVITO analysis and computing core angle, eccentricity, area via fitting elastic stability parameter distributions.

## Files

### `get_ovito.py`
Orchestrates computations of the elastic stability parameter and CNA, DXA analysis, core angle, eccentricity, area.

**Usage:**
```bash
# Run OVITO analysis
python -m dislocate.core_analysis.OvitoFitting.get_ovito --dis_cell dis_cell.poscar --ref_cell ref_cell.poscar --b 2.9 \
              --thickness 1 --tmp_stab tmp_stab --tmp_dxa tmp_dxa --tmp_fitting tmp_fitting \
              --fitting true --oxygen 0 --pbc false --config S
```

### `fitting_core.py`
Computes core angle, eccentricity, area via fitting elastic stability parameter distributions.

**Usage:**
```bash
from dislocate.core_analysis.OvitoFitting.fitting_core import fit_core_properties

# Fit core properties
python -m dislocate.core_analysis.OvitoFitting.fitting_core --outStab tmp_stab --outDXA tmp_dxa \
                   --thickness 1 --outFitting tmp_fitting --pbc false
```

### `ovito_cna.py`
Common Neighbor Analysis implementation. Outputs the phase distribution (number of bcc atoms, fcc atoms, hcp atoms ico atoms and others). This script is not included in the `analyze_core.py` routine.

## Output Formats

```
bcc fcc hcp ico other
0 0 40 0 0
```

### `ovito_dxa.py`
DXA (Dislocation Extraction Algorithm) implementation.

## Output Formats

```
30.22 40.45 xDXA yDXA leftDis
10.50 13.24 xDXA yDXA rightDis
```

### `ovito_elastStab.py`
Elastic stability parameter and CNA calculation.

## Output Formats
#reference cell x, y, z - dislocation cell x, y, z - stab par - cna

```
0 0 0 0.12 0.02 0.08 0.85 1
```
