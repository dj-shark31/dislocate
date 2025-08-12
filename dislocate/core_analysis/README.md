# core_analysis: Dislocation Core Analysis Workflow

The `core_analysis` directory contains a modular Python workflow for analyzing dislocation core structures in crystalline materials. It automates the process of computing different core analysis parameter and compiling them into one output file.

## Directory Structure

- `analyze_core.py`: Orchestrates the full workflow for multiple dislocation configurations.
- `get_data.py`: Runs the workflow for a single configuration (used internally).
- `assemble.py`: Assembles results from all analysis steps into a final output file.
- `Babel/`, `EnergyStress/`, `OvitoFitting/`, `Pattern/`: Submodules for specific analysis tasks.
- `input_file_template`: Example input file for running the workflow.

## How the Workflow Operates

1. **Input Preparation:**  
   - Prepare a reference POSCAR file (`ref_cell`) and one or more dislocation POSCAR files (`dis_cells`). Note that atom ordering in the POSCAR file matters for the elastic stability paramter computation. If the dislocation poscar file is wrongly ordered, the user can input reference dislocation poscar file that is correctly ordered (under `ref_dis_cells`)  to correctly re-order atoms.
   - Create an input file (see template below) specifying all required parameters and file paths. The user can also pass arguments directly through the command line.

2. **Running the Workflow:**  
     ```bash
     python -m dislocate.core_analysis.analyze_core input_file ncore
     ```

     OR

     ```bash
     python -m dislocate.core_analysis.analyze_core --ncore 4 --thickness 4 --ref_cell POSCAR --output_files output1.dat output2.dat --dis_cells _dislocation1.poscar _dislocation2.poscar --energy_stress true --fitting true --ovito true --nye true --sf true --nrep 1 --oxygen 0 --pbc false --config S --nx 32 --potential_path /path/to/potential.meam --potential_type MEAM --ref_dis_cells ref_dislocation1.poscar ref2_dislocation1.poscar
     ```

     - `input_file`: Path to your input file (see template). Path provided in input_file can be relative to the working directory
     - `ncore`: Number of CPU cores to use for parallel processing.

3. **Analysis Steps (Automated):**
   - **Pattern Generation:** Uses Babel to identifies patterns for the Nye tensor computation.
   - **Energy and Stress Calculation:** Computes energies and stresses if enabled.
   - **OVITO Analysis:** Computes elastic stability parameter, performs a Common Neighbor Analysis, extract dislocation position via a DXA analysis, and extracts dislocation positions.
   - **Nye Tensor Calculation:** Computes Nye tensor if enabled.
   - **Fitting:** Computes core angle, core eccentricity, and core area.
   - **Assembly:** Combines all results into output files, one per dislocation configuration.

4. **Output:**  
   - Results are written to the files specified in the `output_files` field of your input file.

## Minimal Example Command

```bash
python main.py input_file 4
```

- This will process all dislocation configurations listed in your input file using 4 CPU cores.


## Requirements

- Python 3.x
- ASE, NumPy, Pandas, Matplotlib, and other scientific Python libraries
- Modified OVITO (for elastic stability parameter, CNA, DXA). The program automatically uses a remote Docker container (`dj31/ovito-wallace:amd64`) that contains the modified OVITO if possible (if Docker is installed). The other option is for the user to generate a Singularity container from the existing Docker container ()`singularity pull ovitowallace.sif docker://dj31/ovito-wallace:amd64`) and specify the sif file in `config.yaml`. 
- Babel (for Nye tensor calculation). Specify paths to the Babel executables (`babel`, `displacement`, `patternDetect`, and `patternInit`) in `config.yaml`. 