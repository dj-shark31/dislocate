#!/usr/bin/env python3
"""
Python version of analyze_core.sh: orchestrates the full workflow.
Assumes both ref_cell and dis_cells are POSCAR files. Handles nrep > 1 by replicating ref_cell along z using ToPoscar.py.

Example usage:
    python analyze_core.py input_file 4

This script performs dislocation analysis by:
1. Reading configuration from input file (input_file)
2. Processing reference and dislocation cell files (reordering if needed) in parallel if ncore > 1
3. Running displacement, energy/stress, and Nye tensor calculations
4. Generates output files for each configuration

Required input file parameters:
- thickness: Cell thickness (number of unit cells in z direction)
- output_files: Space-separated list of output files for each configuration
- ref_cell: Reference POSCAR file for lattice parameters
- dis_cells: Space-separated list of dislocation cell files
- ref_dis_cells: Reference dislocation cells with matching atom order
- energy_stress: Whether to run energy and stress calculations
- nrep: Number of times to replicate reference POSCAR along z
- ovito: Whether to run OVITO analysis
- nye: Whether to run Nye tensor analysis
- sf: Whether to include stacking fault patterns
- oxygen: Oxygen flag for ToPoscar.py (0=keep, 1=remove)
- potential: Interatomic potential type for LAMMPS
- pbc: Whether to use periodic boundary conditions
- disposcar: Whether to use a dislocated POSCAR
- config: Configuration string (e.g. S, O)
- nx: Number of unit cells in x direction
- potential_type: Potential for energy/stress calculation
- potential_path: Path to potential file

Example usage:
    python analyze_core.py --input_file input_file --ncore 4
    python analyze_core.py --input_file input_file --ncore 4 --thickness 4 --ref_cell POSCAR --output_files output1.dat output2.dat --dis_cells dis1.poscar dis2.poscar --energy_stress true --fitting true --ovito true --nye true --sf true --nrep 1 --oxygen 0 --pbc false --config S --nx 32 --potential_path /path/to/potential.meam --potential_type MEAM --ref_dis_cells ref1.poscar ref2.poscar
"""
import argparse
import os
import tempfile
import subprocess
import shlex
from concurrent.futures import ProcessPoolExecutor
from ase.io import read
from ase.atoms import Atoms
import numpy as np
from ase.geometry import get_distances

def parse_input_file(path):
    """Parse key=value pairs from the input file."""
    config = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                config[key.strip()] = value.strip().strip('"')
    return config

def get_lattice_params_and_natoms_from_poscar(poscar_file, thickness, nx):
    atoms = read(poscar_file, format='vasp')
    cell = atoms.cell
    a0 = cell.lengths()[2] / thickness  # z direction
    c0 = cell.lengths()[0] / nx         # x direction
    coa0 = c0 / a0
    natom = len(atoms)
    return a0, c0, coa0, natom

def run(cmd, **kwargs):
    print('Running:', ' '.join(str(x) for x in cmd))
    return subprocess.run(cmd, check=True, **kwargs)

def abspath_from_script(rel_path):
    """Return absolute path relative to the directory containing this script."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, rel_path))

def reorder_atoms_to_reference(atoms: Atoms, ref_atoms: Atoms) -> Atoms:
    """
    Reorder atoms to match reference structure with minimal distances.
    
    Args:
        atoms: Atoms object to reorder
        ref_atoms: Reference Atoms object
        
    Returns:
        Reordered Atoms object
    """
    if len(atoms) != len(ref_atoms):
        raise ValueError("Number of atoms must match between structures")
    
    # Calculate distance matrix
    pos1 = ref_atoms.get_positions()
    pos2 = atoms.get_positions()
    
    _, distances = get_distances(pos1, pos2, ref_atoms.cell, pbc=True)
    
    row_indices = []
    for i in range(len(atoms)):
        min_idx = np.argmin(distances[i])
        if min_idx in row_indices:
            # If duplicate found, get next closest atom not already used
            sorted_indices = np.argsort(distances[i])
            for idx in sorted_indices:
                if idx not in row_indices:
                    min_idx = idx
                    break
        row_indices.append(min_idx)
    # Create new atoms object with reordered positions
    reordered_atoms = atoms.copy()
    
    # Reorder atomic positions
    new_positions = atoms.get_positions()[row_indices]
    reordered_atoms.set_positions(new_positions)
    
    # Reorder atomic numbers and other properties if they exist
    if hasattr(atoms, 'get_atomic_numbers'):
        new_numbers = atoms.get_atomic_numbers()[row_indices]
        reordered_atoms.set_atomic_numbers(new_numbers)
    
    # Reorder tags if they exist
    if hasattr(atoms, 'get_tags'):
        new_tags = atoms.get_tags()[row_indices]
        reordered_atoms.set_tags(new_tags)
    
    # Reorder momenta if they exist
    if hasattr(atoms, 'get_momenta'):
        new_momenta = atoms.get_momenta()[row_indices]
        reordered_atoms.set_momenta(new_momenta)
    
    return reordered_atoms

def main():
    parser = argparse.ArgumentParser(description='Python version of main.sh (POSCAR-only version)')
    parser.add_argument('--input_file', help='Input file with variables')
    parser.add_argument('--ncore', type=int, help='Number of cores for parallel execution', default=1)
    parser.add_argument('--thickness', type=int, help='Cell thickness (number of unit cells in z direction)', default=1)
    parser.add_argument('--ref_cell', help='Reference POSCAR file')
    parser.add_argument('--output_files', help='Output files for each dislocation configuration (space-separated)')
    parser.add_argument('--dis_cells', help='Dislocation cell files (space-separated)')
    parser.add_argument('--energy_stress', help='Whether to run Energy and Stress calculations', default='true')
    parser.add_argument('--fitting', help='Whether to run fitting', default='true')
    parser.add_argument('--ovito', help='Whether to run OVITO analysis', default='true')
    parser.add_argument('--nye', help='Whether to run Nye tensor analysis', default='true')
    parser.add_argument('--sf', help='Whether to include stacking fault patterns', default='false')
    parser.add_argument('--nrep', type=int, help='Number of times to replicate reference POSCAR along z', default=1)
    parser.add_argument('--oxygen', type=int, help='Oxygen flag (0=keep oxygen, 1=remove oxygen)', default=0)
    parser.add_argument('--pbc', help='Whether to use periodic boundary conditions', default='false')
    parser.add_argument('--config', help='Configuration string', default='S')
    parser.add_argument('--nx', type=int, help='Number of unit cells in x direction', default=32)
    parser.add_argument('--potential_path', help='Path to potential file', default=abspath_from_script('../potentials/Ti.meam'))
    parser.add_argument('--potential_type', help='Potential type to use', default='MEAM')
    parser.add_argument('--ref_dis_cells', help='Reference dislocation cell files (space-separated)')
    args = parser.parse_args()

    if args.input_file:
        # Load from input file if provided
        config = parse_input_file(args.input_file)
        ncore = args.ncore
        thickness = int(config.get('thickness', 1))
        ref_cell = config['ref_cell']
        output_files = shlex.split(config['output_files'])
        dis_cells = shlex.split(config['dis_cells'])
        energy_stress = config.get('energy_stress', 'true').lower()
        fitting = config.get('fitting', 'true').lower()
        ovito = config.get('ovito', 'true').lower()
        nye = config.get('nye', 'true').lower()
        sf = config.get('sf', 'true').lower()
        nrep = int(config.get('nrep', 1))
        oxygen = int(config.get('oxygen', 0))
        pbc = config.get('pbc', 'false').lower()
        config_str = config.get('config', 'S')
        nx = int(config.get('nx', 32))
        potential_path = config.get('potential_path', abspath_from_script('../potentials/Ti.meam'))
        potential_type = config.get('potential_type', 'MEAM')
        ref_dis_cells = shlex.split(config.get('ref_dis_cells', ''))
    else:
        # Use command line arguments
        ncore = args.ncore
        thickness = args.thickness
        ref_cell = args.ref_cell
        output_files = shlex.split(args.output_files) if args.output_files else []
        dis_cells = shlex.split(args.dis_cells) if args.dis_cells else []
        energy_stress = args.energy_stress.lower()
        fitting = args.fitting.lower()
        ovito = args.ovito.lower()
        nye = args.nye.lower()
        sf = args.sf.lower()
        nrep = args.nrep
        oxygen = args.oxygen
        pbc = args.pbc.lower()
        config_str = args.config
        nx = args.nx
        potential_path = args.potential_path
        potential_type = args.potential_type
        ref_dis_cells = shlex.split(args.ref_dis_cells) if args.ref_dis_cells else []

    # Setup tmp dir
    tmp_dir = 'tmp'
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_pattern = tempfile.NamedTemporaryFile(prefix='pattern-', dir=tmp_dir, delete=False).name
    tmp_ref_poscar = tempfile.NamedTemporaryFile(prefix='ref-', dir=tmp_dir, delete=False).name

    if ref_dis_cells:
        # Create temporary files for reordered cells
        # Verify that ref_dis_cells and dis_cells have the same length
        if len(ref_dis_cells) != len(dis_cells):
            raise ValueError(f"Number of reference dislocation cells ({len(ref_dis_cells)}) must match number of dislocation cells ({len(dis_cells)})")

        reordered_dis_cells = []
        for dis_cell, ref_dis_cell in zip(dis_cells, ref_dis_cells):
            # Read the cells using ASE
            dis_atoms = read(dis_cell)
            ref_dis_atoms = read(ref_dis_cell)
            
            # Reorder atoms according to reference
            reordered_atoms = reorder_atoms_to_reference(dis_atoms, ref_dis_atoms)
            
            # Create temp file for reordered cell
            tmp_dis_cell = tempfile.NamedTemporaryFile(prefix='reordered-', dir=tmp_dir, delete=False).name
            reordered_atoms.write(tmp_dis_cell, format='vasp')
            reordered_dis_cells.append(tmp_dis_cell)
        
        # Update dis_cells to use reordered versions
        print("Dislocation cells reordered to match reference dislocation cells")
        dis_cells = reordered_dis_cells

    # Reference POSCAR replication (if nrep > 1)
    if nrep > 1:
        # Use ToPoscar.py to replicate along z-direction
        run(['python3', abspath_from_script('ToPoscar.py'), ref_cell, tmp_ref_poscar, str(oxygen), str(nrep)])
        ref_cell = tmp_ref_poscar
        thickness *= nrep

    # Lattice params and natom from POSCAR using ASE
    a0, c0, coa0, natom = get_lattice_params_and_natoms_from_poscar(ref_cell, thickness, nx)

    # Pattern detection
    if sf:
        run(['python3', abspath_from_script('Pattern/get_patternInit.py'), str(a0), str(coa0), tmp_pattern])
    else:
        run(['python3', abspath_from_script('Pattern/get_pattern.py'), ref_cell, str(thickness), str(a0), str(natom), tmp_pattern])

    # Prepare argument sets for get_data.py
    job_args = []
    for dis_cell, output_file in zip(dis_cells, output_files):
        job_args.append([
            '--ref_cell', ref_cell, '--thickness', thickness, '--a0', a0, 
            '--natom', natom, '--tmp_pattern', tmp_pattern,
            '--energy_stress', energy_stress, '--fitting', fitting, '--ovito', ovito, '--nye', nye,
            '--oxygen', oxygen, '--pbc', pbc,
            '--dis_cell', dis_cell, '--output_file', output_file, '--config', config_str,
            '--potential_type', potential_type, '--potential_path', potential_path
        ])

    # Parallel or serial execution
    if ncore > 1:
        print("Running in parallel with", ncore, "cores")
        with ProcessPoolExecutor(max_workers=ncore) as executor:
            futures = [executor.submit(run, [str(x) for x in ['python3', abspath_from_script('get_data.py')] + args]) for args in job_args]
            for f in futures:
                f.result()
    else:
        for args_ in job_args:
            print("Running in serial")
            run([str(x) for x in ['python3', abspath_from_script('get_data.py')] + args_])

    # Clean up tmp files
    for f in [tmp_pattern, tmp_ref_poscar]:
        try:
            os.remove(f)
            print(f"Removed {f}")
        except Exception:
            pass
    # Clean up tmp_dis_cells if defined
    if 'reordered_dis_cells' in locals():
        for f in reordered_dis_cells:
            try:
                os.remove(f)
                print(f"Removed {f}")
            except Exception:
                pass

if __name__ == '__main__':
    main() 