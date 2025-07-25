# core_analysis: Dislocation Core Analysis Workflow

The `core_analysis` directory contains a modular Python workflow for analyzing dislocation core structures in crystalline materials. It automates the process of extracting, fitting, and assembling core properties from atomistic simulation data, using a combination of LAMMPS, OVITO, Babel and custom analysis scripts.

## Key Features

- **Automated Workflow:** Orchestrates the full analysis pipeline from input structure files to final core property output.
- **Flexible Analysis:** Supports LAMMPS calculations, OVITO-based analysis, Nye tensor computation, and custom fitting routines.
- **Parallel Execution:** Can process multiple dislocation configurations in parallel.
- **Extensible:** Modular scripts for each analysis step (LAMMPS, OVITO, Babel, fitting, assembly).

## Directory Structure

- `main.py`: Orchestrates the full workflow for multiple dislocation configurations.
- `get_data.py`: Runs the workflow for a single configuration (used internally).
- `assemble.py`: Assembles results from all analysis steps into a final output file.
- `Babel/`, `EnergyStress/`, `OvitoFitting/`, `Pattern/`: Submodules for specific analysis tasks.
- `input_file_template`: Example input file for running the workflow.

## How the Workflow Operates

1. **Input Preparation:**  
   - Prepare a reference POSCAR file and one or more dislocation POSCAR files.
   - Create an input file (see template below) specifying all required parameters and file paths.

2. **Running the Workflow:**  
   - From the `core_analysis` directory, run:
     ```bash
     python main.py input_file ncore
     ```
     - `input_file`: Path to your input file (see template below). Path provided in input_file can be relative to the working directory
     - `ncore`: Number of CPU cores to use for parallel processing.

3. **Analysis Steps (Automated):**
   - **Pattern Generation:** Generates stacking fault patterns if requested.
   - **Energy and Stress Calculation:** Computes energies and stresses if enabled.
   - **OVITO Analysis:** Extracts dislocation positions and atomic properties.
   - **Nye Tensor Calculation:** Computes Nye tensor fields if enabled.
   - **Fitting:** Fits core properties from the extracted data.
   - **Assembly:** Combines all results into output files, one per dislocation configuration.

4. **Output:**  
   - Results are written to the files specified in the `output_files` field of your input file.

## Minimal Example Command

```bash
python main.py input_file 4
```

- This will process all dislocation configurations listed in your input file using 4 CPU cores.

## Advanced Usage

- **Custom Analysis:** Modify or extend scripts in the `Babel/`, `EnergyStress/`, `OvitoFitting/`, or `Pattern/` subdirectories for specialized analysis.

## Requirements

- Python 3.x
- ASE, NumPy, Pandas, Matplotlib, and other scientific Python libraries
- OVITO (for structure analysis, via Singularity container if needed) 
- Babel (for Nye tensor calculation)