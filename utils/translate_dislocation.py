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
from ase.geometry import get_distances
from ase.atoms import Atoms
from utils.atomistic_tools import get_unit_cell_vectors

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Translate, wrap, and reorder atoms using ASE')
    parser.add_argument('--input_cell',required=True, help='Input cell structure file')
    parser.add_argument('--output_file',required=True, help='Output structure file')
    parser.add_argument('--translation', nargs=3, type=float, default=[0, 0, 0],
                       help='Translation vector in lattice units (default: 0 0 0)')
    parser.add_argument('--n_cells', nargs=3, type=int, default=[8, 8, 1],
                       help='Number of unit cells in x-, y-, and z-direction (default: 8 8 1)')
    return parser.parse_args()

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

def translate_dislocation(atoms, translation, n_cells):
    """Translate dislocation by a lattice vector."""
    translated_atoms = atoms.copy()

    # Step 1: Translate atoms by a lattice vector
    a1, a2, a3 = get_unit_cell_vectors(atoms, n_cells)
    translation = a1 * translation[0] + a2 * translation[1] + translation[2] * a3
    print(f" Translating atoms by vector: {translation} (Cartesian coordinates)")
    translated_atoms.translate(translation)
    
    # Step 2: Wrap positions to unit cell
    translated_atoms.wrap()
    
    # Step 3: Reorder atoms to match reference
    print(f" Reordering atoms to match reference structure")
    reordered_atoms = reorder_atoms_to_reference(translated_atoms, atoms)
    return reordered_atoms

def main():
    args = parse_arguments()
    
    # Check if input files exist
    if not os.path.exists(args.input_cell):
        raise FileNotFoundError(f"Reference cell file not found: {args.input_cell}")
    
    print("=" * 60)
    print("Dislocation dipole translation and reordering")
    print("=" * 60)

    # Read structures
    print(f" Reading reference cell: {args.input_cell}")
    reference_atoms = read(args.input_cell)

    reordered_atoms = translate_dislocation(atoms=reference_atoms, translation=args.translation, n_cells=args.n_cells)
    
    print(f" Writing output structure: {args.output_file}")
    write(args.output_file, reordered_atoms, format='vasp')
    
    print(f"=" * 60)
    print("Processing complete!")
    print(f"Output saved to: {args.output_file}")
    print("=" * 60)

if __name__ == "__main__":
    main() 
