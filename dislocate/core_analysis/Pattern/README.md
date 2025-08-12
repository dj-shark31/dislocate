# Pattern Detection and Analysis

This directory contains scripts for detecting patterns in dislocation structures. It generates a file required to run the nye tensor analysis.

## Files

### `get_pattern.py`
Doesn't include the detection of stacking faults when computing the Nye tensor.

**Usage:**
```bash
# Detect patterns in structure
python -m dislocate.core_analysis.Pattern.get_pattern --ref_cell ref_cell.poscar --thickness 1 --a0 2.9 --natom 192 --tmp_pattern tmp_pattern
```

### `get_patternInit.py`
Include the detection of stacking faults when computing the Nye tensor.

**Key Functions:**
- `initialize_patterns()`: Initialize pattern detection
- `setup_pattern_parameters()`: Setup pattern parameters
- `validate_pattern_setup()`: Validate pattern setup
- `optimize_pattern_detection()`: Optimize pattern detection

**Usage:**
```bash
# Detect patterns in structure
python -m dislocate.core_analysis.Pattern.get_patternInit --a0 2.9 --coa0 1.57 --tmp_pattern tmp_pattern
```

### `pattern.dat`
Pattern data file template.

### `patternInit.dat`
Pattern initialization data.
