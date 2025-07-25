import sys
import numpy as np
import argparse
from ase.io import read, write
from ase.build import make_supercell

#Set up argument parser
parser = argparse.ArgumentParser(description='Convert LAMMPS data file to POSCAR format using ASE')
parser.add_argument('infile', help='Input LAMMPS data file')
parser.add_argument('outfile', help='Output POSCAR file')
parser.add_argument('oxygen', type=int, choices=[0, 1], 
                   help='Delete oxygen atoms (0= keep oxygen, 1= delete oxygen)')
parser.add_argument('nrep', type=int, default=1, help='Number of times to replicate the slab')

#Parse arguments
args = parser.parse_args()

#Import the structure using ASE
atoms = read(args.infile)

#Remove oxygen atoms if requested (type 2, which is typically oxygen)
if args.oxygen == 1:
    try:
        # Get atom types and remove oxygen (type 2)
        if hasattr(atoms, 'get_atomic_numbers'):
            # Alternative approach: remove atoms by atomic number (8 for oxygen)
            atomic_numbers = atoms.get_atomic_numbers()
            mask = atomic_numbers != 8  # Remove oxygen atoms
            atoms = atoms[mask]
        else:
            print("Warning: Cannot identify atom types for oxygen removal")
    except Exception as e:
        print(f"Warning: Could not remove oxygen atoms: {e}")

# Create a supercell matrix with the specified thickness
supercell_matrix = np.array([[1, 0, 0],
                            [0, 1, 0], 
                            [0, 0, args.nrep]])
atoms = make_supercell(atoms, supercell_matrix)

#Get particle identifiers if available and sort atoms
try:
    if hasattr(atoms, 'arrays') and 'id' in atoms.arrays:
        # Sort atoms by their ID
        ids = atoms.arrays['id']
        order = np.argsort(ids)
        atoms = atoms[order]
except Exception as e:
    print(f"Warning: Could not sort atoms by ID: {e}")

#Write outputfile using ASE's write function
write(args.outfile, atoms, format='vasp')
