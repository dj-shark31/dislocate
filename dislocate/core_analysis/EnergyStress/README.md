# Energy and Stress Analysis

This directory contains scripts for calculating energy and stress fields in dislocation structures.

## Files

### `get_energy_stress.py`
Computes energy and stress tensor of cell.

**Usage:**
```bash
# Calculate energy field
python -m dislocate.core_analysis.EnergyStress.get_energy_stress --dis_cell dislocated_cell.poscar --tmp_energy_stress tmp_energy_stress --potential_type MACE --potential_path /path/to/mace_model.pt
```

## Output Formats

```
10 20 30 0 0 0 sxx syy szz syz sxz sxy (MPa)
1500.99 toteng (eV)
```