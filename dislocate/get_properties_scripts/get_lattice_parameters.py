import argparse
import os
import numpy as np
import sys

# Add project root to Python path for subprocess compatibility
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from dislocate.utils.atomistic_tools import alloy_randomly, load_sqs, build_atoms, get_lattice_parameters
from dislocate.utils.submitruns import composition_dir

def main():
    parser = argparse.ArgumentParser(description="Get lattice parameters for a given composition")
    parser.add_argument('--structure', type=str, required=True, help='Crystal structure (hcp/bcc/fcc/omega)')
    parser.add_argument('--elements', nargs='+', type=str, required=True, help='Elements (Ti/V/Cr/Mn/Fe/Co/Ni)')
    parser.add_argument('--composition', nargs='+', type=float, required=True, help='Composition (0.5/0.75/1.0)')
    parser.add_argument('--repetitions', type=int, required=True, help='Number of repetitions')
    parser.add_argument('--base_dir', type=str, required=True, help='Base directory')
    parser.add_argument('--potential_path', type=str, required=True, help='Path to potential file')
    parser.add_argument('--potential_type', type=str, required=True, help='Potential type (MACE/MEAM/...)')
    parser.add_argument('--device', type=str, required=True, help='Device (cpu/gpu)')
    parser.add_argument('--fmax', type=float, required=True, help='Maximum force')
    parser.add_argument('--dof', nargs='+', type=bool, required=True, help='Degrees of freedom (x/y/z/rx/ry/rz)')
    parser.add_argument('--n_cells', nargs='+', type=int, required=True, help='Number of cells (x/y/z)')
    parser.add_argument('--sqs', type=str, default="false", help='SQS order parameter')
    parser.add_argument('--output_cell', type=str, default="true", help='Output relaxed cell')
    
    args = parser.parse_args()
    output_dir = os.path.join(args.base_dir, composition_dir(args.elements, args.composition))

    lattice_params = []
    for i in range(args.repetitions):
        # Load SQS structure if it exists or build randomly alloyed structure
        if args.sqs == "true":
            atoms = load_sqs(args.elements, args.composition, args.structure, i)
            output_cell = os.path.join(output_dir, f"{args.structure}_nx{args.n_cells[0]}_ny{args.n_cells[1]}_nz{args.n_cells[2]}_sqs_{i}.poscar")
        else:
            atoms = build_atoms(args.structure, args.n_cells)
            atoms = alloy_randomly(atoms, args.elements, args.composition)
            output_cell = os.path.join(output_dir, f"{args.structure}_nx{args.n_cells[0]}_ny{args.n_cells[1]}_nz{args.n_cells[2]}_random_{i}.poscar")

        # Get lattice parameters
        lattice_params, atoms = get_lattice_parameters(atoms, n_cells = args.n_cells, potential_path = args.potential_path, potential_type = args.potential_type, device = args.device, dof = args.dof, fmax = args.fmax, phase = args.structure)
        if args.output_cell == "true":
            atoms.write(output_cell, format='vasp')
        lattice_params.append(lattice_params)

    # Write results to file
    if args.sqs == "true":
        output_properties = os.path.join(output_dir, f"{args.structure}_nx{args.n_cells[0]}_ny{args.n_cells[1]}_nz{args.n_cells[2]}_lattice_sqs.dat")
    else:
        output_properties = os.path.join(output_dir, f"{args.structure}_nx{args.n_cells[0]}_ny{args.n_cells[1]}_nz{args.n_cells[2]}_lattice_random.dat")

    with open(output_properties, "w") as f:
        # Get all parameter names from first dictionary
        param_names = list(lattice_params[0].keys())
        
        # Write header with parameter names and their std
        header = " ".join([f"{param} {param}_std" for param in param_names])
        f.write(f"{header}\n")
        
        # Write values and standard deviations
        values_str = []
        for param in param_names:
            param_values = [lp[param] for lp in lattice_params]
            values_str.extend([f"{np.mean(param_values):.6f}", f"{np.std(param_values) / np.sqrt(len(param_values)):.6f}"])
        f.write(" ".join(values_str) + "\n")

if __name__ == "__main__":
    main()
