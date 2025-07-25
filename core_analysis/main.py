#!/usr/bin/env python3
"""
Python version of main.sh: orchestrates the full workflow.
Assumes both ref_cell and dis_cells are POSCAR files. Handles nrep > 1 by replicating ref_cell along z using ToPoscar.py.
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
    parser.add_argument('input_file', help='Input file with variables')
    parser.add_argument('ncore', type=int, help='Number of cores for parallel execution')
    args = parser.parse_args()

    config = parse_input_file(args.input_file)
    ncore = args.ncore
    thickness = int(config.get('thickness', 1))
    ref_cell = config['ref_cell']
    output_files = shlex.split(config['output_files'])
    dis_cells = shlex.split(config['dis_cells'])
    lammps = config.get('lammps', 'true').lower()
    fitting = config.get('fitting', 'true').lower()
    ovito = config.get('ovito', 'true').lower()
    nye = config.get('nye', 'true').lower()
    sf = config.get('sf', 'true').lower()
    nrep = int(config.get('nrep', 1))
    oxygen = int(config.get('oxygen', 0))
    pbc = config.get('pbc', 'false').lower()
    config_str = config.get('config', 'S')
    nx = int(config.get('nx', 32))
    potential_path  = config.get('potential_path', '/global/scratch/users/djany/lammps/potentials/Ti.meam')
    potential_type = config.get('potential_type', 'MEAM')
    ref_dis_cells = shlex.split(config.get('ref_dis_cells', ''))

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
            '--lammps', lammps, '--fitting', fitting, '--ovito', ovito, '--nye', nye,
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