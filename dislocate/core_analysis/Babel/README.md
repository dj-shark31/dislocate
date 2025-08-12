# Babel Integration

This directory contains scripts for to compute the nye tensor with Babel 

## Files

### `get_babel.py`
Main interface for Babel calculations.

**Usage:**
```bash
# Run Babel analysis
python -m dislocate.core_analysis.Babel.get_babel --dis_cell dis_cell.poscar --ref_cell ref_cell.poscar \
        --thickness 1 --a0 2.9 --natom 192 --file_pattern tmp_pattern --tmp_babel tmp_babel --oxygen 0
```

### `input_nye.displacement`
Input file template for Babel displacement calculations.

## Output Formats

```
# x y z α_11 α_12 α_13 α_21 α_22 α_23 α_31 α_32 α_33 (header not in file)
0.0 0.0 0.0 0.123 0.045 0.067 0.045 0.234 0.089 0.067 0.089 0.345
```