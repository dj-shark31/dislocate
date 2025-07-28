#!/usr/bin/env python3
"""
Python equivalent of get_energy_stress.sh using ASE LAMMPSlib
Runs EnergyStress calculations using ASE's LAMMPSlib calculator or MACE calculator
"""

import argparse
from pathlib import Path
from ase.io import read
from ase.calculators.lammpslib import LAMMPSlib
from mace.calculators import MACECalculator

def setup_lammpslib_calculator(potential_type, potential_path):
    """Setup LAMMPSlib calculator with the specified potential"""

    if potential_type == "MEAM":
        lmp_pair_style = 'meam/spline'
        lmp_pair_coeff = [f'* * {potential_path} Ti']
    elif potential_type == "DMD":
        lmp_pair_style = f'deepmd {potential_path}'
        lmp_pair_coeff = ['* *']
    elif potential_type == "RANN":
        lmp_pair_style = 'rann'
        lmp_pair_coeff = [f'* * {potential_path} Ti']
    elif potential_type == "ACE":
        lmp_pair_style = 'pace product'
        lmp_pair_coeff = [f'* * {potential_path} Ti']
    else:
        raise ValueError(f"Unsupported potential type: {potential_type}")
 
    lmp_cmds = [lmp_pair_style, lmp_pair_coeff]
    # LAMMPSlib calculator setup
    calc = LAMMPSlib(lmpcmds=lmp_cmds)
    return calc

def get_energy_stress(atoms, potential_type, potential_path):
    """Run calculation using ASE LAMMPSlib or MACE calculator"""
    if potential_type in ["MEAM", "DMD", "RANN", "ACE"]:
        calc = setup_lammpslib_calculator(potential_type, potential_path)
    elif potential_type == "MACE":
        calc = MACECalculator(model_paths=potential_path, device='cpu')
    else:
        raise ValueError(f"Unsupported potential type: {potential_type}")
    atoms.calc = calc
    try:
        energy = atoms.get_potential_energy()
        stress = atoms.get_stress(voigt=True)  # [xx, yy, zz, yz, xz, xy]
        # Convert stress from eV/Å³ to MPa
        stress_mpa = stress * 160217.66208
        return stress_mpa, energy
    except Exception as e:
        print(f"Error in LAMMPSlib calculation: {e}")
        return None, None
    finally:
        if hasattr(calc, 'clean'):
            calc.clean()

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
        stress, energy = get_energy_stress(atoms, args.potential_type, args.potential_path)
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