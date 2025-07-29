#!/usr/bin/env python3
"""
Python version of get_data.sh: orchestrates the workflow for a single dislocation configuration.
Assumes both ref_cell and dis_cell are POSCAR files only.
"""
import argparse
import os
import tempfile
from analyze_core import run, abspath_from_script

def main():
    parser = argparse.ArgumentParser(description='Orchestrate workflow for a single dislocation configuration (POSCAR-only, keyword arguments)')
    parser.add_argument('--ref_cell', required=True, help='Reference POSCAR file')
    parser.add_argument('--thickness', required=True, help='Cell thickness (number of unit cells in z)')
    parser.add_argument('--a0', required=True, help='Lattice parameter a0 (z direction)')
    parser.add_argument('--natom', required=True, help='Number of atoms in the reference cell')
    parser.add_argument('--tmp_pattern', required=True, help='Pattern file (output from get_patternInit.py or get_pattern.py)')
    parser.add_argument('--energy_stress', required=True, help='Whether to run EnergyStress calculations (true/false)')
    parser.add_argument('--fitting', required=True, help='Whether to run fitting (true/false)')
    parser.add_argument('--ovito', required=True, help='Whether to run OVITO analysis (true/false)')
    parser.add_argument('--nye', required=True, help='Whether to run Nye tensor analysis (true/false)')
    parser.add_argument('--oxygen', required=True, help='Oxygen flag for ToPoscar.py (0=keep, 1=remove)')
    parser.add_argument('--pbc', required=True, help='Whether to use periodic boundary conditions (true/false)')
    parser.add_argument('--dis_cell', required=True, help='Dislocation POSCAR file')
    parser.add_argument('--output_file', required=True, help='Output file for results')
    parser.add_argument('--config', required=True, help='Configuration string (e.g., S, O, etc.)')
    parser.add_argument('--potential_type', required=False , choices=['MEAM', 'DMD', 'RANN', 'ACE', 'MACE'], help='Potential type to use')
    parser.add_argument('--potential_path', required=False, help='Path to potential file')
    args = parser.parse_args()

    # Create temp files
    tmp_dir = 'tmp'
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_babel = tempfile.NamedTemporaryFile(prefix='babel-', dir=tmp_dir, delete=False).name
    tmp_stab = tempfile.NamedTemporaryFile(prefix='stab-', dir=tmp_dir, delete=False).name
    tmp_dxa = tempfile.NamedTemporaryFile(prefix='dxa-', dir=tmp_dir, delete=False).name
    tmp_fitting = tempfile.NamedTemporaryFile(prefix='fitting-', dir=tmp_dir, delete=False).name
    tmp_energy_stress = tempfile.NamedTemporaryFile(prefix='energy_stress-', dir=tmp_dir, delete=False).name

    # Nye computation
    if args.nye == 'true':
        print("Nye computation started")
        run(['python3', abspath_from_script('Babel/get_babel.py'),
             args.dis_cell, args.ref_cell, args.thickness, args.a0, args.natom,
             args.tmp_pattern, tmp_babel, args.oxygen])
        print("Nye computation ended")

    # Ovito/fitting computation
    if args.ovito == 'true' or args.fitting == 'true':
        print("Ovito computation started")
        run(['python3', abspath_from_script('OvitoFitting/get_ovito.py'),
             args.dis_cell, args.ref_cell, args.a0, args.thickness,
             tmp_stab, tmp_dxa, tmp_fitting,
             '--fitting', args.fitting,
             '--oxygen', args.oxygen,
             '--pbc', args.pbc,
             '--config', args.config
            ])
        print("Ovito computation ended")

    # Lammps computation
    if args.energy_stress == 'true':
        print("EnergyStress computation started")
        run(['python3', abspath_from_script('EnergyStress/get_energy_stress.py'),
             args.dis_cell, tmp_energy_stress, args.potential_type, args.potential_path])
        print("EnergyStress computation ended")

    # Assemble data
    print("Data assemble started")
    run(['python3', abspath_from_script('assemble.py'),
            args.thickness, args.a0, args.natom, tmp_babel, tmp_stab, tmp_dxa, tmp_fitting, tmp_energy_stress,
            args.output_file, args.energy_stress, args.fitting, args.ovito, args.nye, args.pbc])
    print("Data assemble ended")

    # Clean up temp files
    for f in [tmp_babel, tmp_stab, tmp_dxa, tmp_fitting, tmp_energy_stress]:
        try:
            os.remove(f)
        except Exception:
            pass

if __name__ == '__main__':
    main() 