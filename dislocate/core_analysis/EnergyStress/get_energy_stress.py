#!/usr/bin/env python3
"""
Python equivalent of get_energy_stress.sh using ASE LAMMPSlib
Runs EnergyStress calculations using ASE's LAMMPSlib calculator or MACE calculator
"""

import argparse
from pathlib import Path
from ase.io import read
import os
import sys

# Add project root to Python path for subprocess compatibility
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))  # Go up two levels
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from dislocate.utils.atomistic_tools import set_calculator

def main():
    parser = argparse.ArgumentParser(description='Run EnergyStress calculations using ASE LAMMPSlib or MACE calculator')
    parser.add_argument('dis_cell', help='Dislocation cell POSCAR file')
    parser.add_argument('tmp_energy_stress', help='Output file for EnergyStress results')
    parser.add_argument('potential_type', choices=['MEAM', 'DMD', 'RANN', 'ACE', 'MACE'],
                       help='Potential type to use')
    parser.add_argument('potential_path', help='Path to potential file')
    args = parser.parse_args()

    Path('tmp').mkdir(exist_ok=True)

    try:
        atoms = read(args.dis_cell, format='vasp')
        set_calculator(atoms, potential_path=args.potential_path, potential_type=args.potential_type, device='cpu')
        energy, stress = atoms.get_potential_energy(), atoms.get_stress(voigt=True) * 160217.66208 
        # Energy in eV, stress in MPa [xx, yy, zz, yz, xz, xy]
        with open(args.tmp_energy_stress, 'w') as f:
            if stress is not None:
                stress_str = ' '.join([f"{s:.6f}" for s in stress])
                f.write(f"{stress_str} sxx syy szz syz sxz sxy (MPa)\n")
            if energy is not None:
                f.write(f"{energy:.6f} toteng (eV)\n")
    except Exception as e:
        print(f"Error in main execution: {e}")
        with open(args.tmp_energy_stress, 'w') as f:
            pass

if __name__ == "__main__":
    main() 