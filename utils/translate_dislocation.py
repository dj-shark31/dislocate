#!/usr/bin/env python3
"""
Script to translate dislocation dipole by a lattice vector. 
All atoms are translated by lattice vectors, wrapped in periodic boundaries, 
and reordered to match a reference structure with minimal distances.

Example command line: python translate_dislocation.py --input_cell POSCAR_reference --output_file POSCAR_translated --translation 0.5 0.5 0.5 --n_cells 8 8 1
"""

import numpy as np
import argparse
import os
from ase.io import read, write
from core_analysis.analyze_core import reorder_atoms_to_reference

def main():
    parser = argparse.ArgumentParser(description='Translate, wrap, and reorder atoms using ASE')
    parser.add_argument('--input_cell',required=True, help='Input cell structure file')
    parser.add_argument('--output_file',required=True, help='Output structure file')
    parser.add_argument('--translation', nargs=3, type=float, default=[0, 0, 0],
                       help='Translation vector in lattice units (default: 0 0 0)')
    parser.add_argument('--n_cells', nargs=3, type=int, default=[8, 8, 1],
                       help='Number of unit cells in x-, y-, and z-direction (default: 8 8 1)')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.input_cell):
        raise FileNotFoundError(f"Reference cell file not found: {args.input_cell}")
    
    print("=" * 60)
    print("Dislocation dipole translation and reordering")
    print("=" * 60)
    
    # Read structures
    print(f"2. Reading reference cell: {args.input_cell}")
    reference_atoms = read(args.input_cell)
    dislocation_atoms = reference_atoms.copy()
    
    vectors = np.array(reference_atoms.cell)
    a1 = vectors[0] / args.n_cells[0]   
    a2 = vectors[1] / args.n_cells[1]
    a3 = vectors[2] / args.n_cells[2]
    
    # Step 1: Translate atoms
    translation = a1 * args.translation[0] + a2 * args.translation[1] + args.translation[2] * a3
    if not np.allclose(translation, 0):
        print(f"4. Translating atoms by vector: {translation} (Cartesian coordinates)")
        dislocation_atoms.translate(translation)
    else:
        print(f"4. No translation applied")
    
    # Step 2: Wrap positions to unit cell
    print(f"5. Wrapping atoms to unit cell")
    dislocation_atoms.wrap()
    reference_atoms.wrap()
    
    # Step 3: Reorder atoms to match reference
    print(f"6. Reordering atoms to match reference structure")
    reordered_atoms = reorder_atoms_to_reference(dislocation_atoms, reference_atoms)
    
    # Step 4: Write output
    print(f"7. Writing output structure: {args.output_file}")
    write(args.output_file, reordered_atoms, format='vasp')
    
    print(f"=" * 60)
    print("Processing complete!")
    print(f"Output saved to: {args.output_file}")
    print("=" * 60)

if __name__ == "__main__":
    main() 
