import argparse
import os
import numpy as np
import sys

# Add project root to Python path for subprocess compatibility
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from utils.atomistic_tools import alloy_randomly, relax_cell, compute_stacking_fault_energy
from utils.submitruns import composition_dir

def build_cell(n_cells):
    a, c = 2.94, 4.64
    lattice_vectors = np.array([
    [a, 0.0, 0.0],  # a-axis
    [0.0, 3**0.5 * a, 0.0],  # b-axis
    [0.0, 0.0, c]   # c-axis
    ])

    # Define atomic positions (fractional coordinates)
    atomic_positions = [
        ("Ti", [0.0, 0.0, 0.0]),  
        ("Ti", [0.5, 0.5, 0.0]),  
        ("Ti", [0.5, 1/6, 0.5]),  
        ("Ti", [0.0, 4/6, 0.5])   
    ]

    # Create ASE Atoms object
    atoms = Atoms(
        symbols=[atom[0] for atom in atomic_positions],  # Atomic symbols
        positions=[np.dot(atom[1], lattice_vectors) for atom in atomic_positions],  # Convert fractional to Cartesian
        cell=lattice_vectors,  # Assign the cell
        pbc=True  # Periodic boundary conditions
    )
    return atoms.repeat((n_cells[0], n_cells[1], n_cells[2]))

def main():
    parser = argparse.ArgumentParser(description="Get lattice parameters for a given composition")
    parser.add_argument('--elements', nargs='+', type=str, required=True, help='Elements (Ti/V/Cr/Mn/Fe/Co/Ni)')
    parser.add_argument('--composition', nargs='+', type=float, required=True, help='Composition (0.5/0.75/1.0)')
    parser.add_argument('--base_dir', type=str, required=True, help='Base directory')
    parser.add_argument('--potential_path', type=str, required=True, help='Path to potential file')
    parser.add_argument('--potential_type', type=str, required=True, help='Potential type (MACE/MEAM/...)')
    parser.add_argument('--device', type=str, required=True, help='Device (cpu/gpu)')
    parser.add_argument('--fmax', type=float, required=True, help='Maximum force')
    parser.add_argument('--n_cell_plane', nargs='+', type=int, required=True, help='Number of cells (x/y)')
    parser.add_argument('--n_layers', type=int, required=True, help='Number of layers (z)')
    parser.add_argument('--plane', type=str, required=True, help='Plane (basal/prismatic)')
    parser.add_argument('--output_cell', type=str, default="true", help='Output relaxed cell')
    
    args = parser.parse_args()
    output_dir = os.path.join(args.base_dir, composition_dir(args.elements, args.composition))

    if args.plane == 'basal':
        n_cells = [args.n_cell_plane[0], args.n_cell_plane[1], args.n_layers]
    elif args.plane == 'prismatic':
        n_cells = [args.n_cell_plane[0], args.n_layers, args.n_cell_plane[1]]
    else:
        raise ValueError(f"Invalid plane: {args.plane}")

    sf_energies = []    # list of stacking fault energies for each layer
    for i in range(args.n_layers):
        output_cell = os.path.join(output_dir, f"{args.plane}_nx{args.n_cell_plane[0]}_ny{args.n_cell_plane[1]}_nz{args.n_layers}_sf_{i}.poscar")
        output_properties = os.path.join(output_dir, f"{args.plane}_nx{args.n_cell_plane[0]}_ny{args.n_cell_plane[1]}_nz{args.n_layers}_sf_energy_{i}.dat")

        atoms = build_cell(n_cells)
        atoms = alloy_randomly(atoms, args.elements, args.composition)
        atoms = relax_cell(atoms, potential_path = args.potential_path, potential_type = args.potential_type, device = args.device, dof = args.dof, fmax = args.fmax)

        energy_reference = atoms.get_potential_energy()
        if args.output_cell == "true":
            atoms.write(output_cell, format='vasp')

        energy = compute_stacking_fault_energy(atoms, args.plane, i, energy_reference, n_cells = args.n_cells, output=output_cell, hard_plane = args.hard_plane, potential_path = args.potential_path, potential_type = args.potential_type, device = args.device, fmax = args.fmax)
        sf_energies.append(energy)

    # Write results to file
    output_properties = os.path.join(output_dir, f"{args.plane}_nx{args.n_cell_plane[0]}_ny{args.n_cell_plane[1]}_nz{args.n_layers}_sf_energy.dat")
    with open(output_properties, "w") as f:    
        # Write header with parameter names and their std
        f.write(f"sf_energies sf_energies_std \n")
        f.write(f"{np.mean(sf_energies):.6f} {np.std(sf_energies) / np.sqrt(len(sf_energies)):.6f}\n")
        f.write(f"layer sf_energy \n")
        for i, energy in enumerate(sf_energies):
            f.write(f"{i} {energy:.6f} \n")

if __name__ == "__main__":
    main()
