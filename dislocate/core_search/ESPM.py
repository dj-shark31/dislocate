import os
import tempfile
import subprocess
import sys
from ase.io import read, write
from dislocate.utils.atomistic_tools import set_calculator
from ase.optimize import BFGS
import argparse
import numpy as np
from dislocate.utils.config_loader import get_tool_path

"""
ESPM (Elastic Singularity Positioning Method) script for dislocation dipole generation and relaxation.

This script performs the following tasks:
1. Takes a reference cell structure (perfect crystal) and generates dislocation configurations using Babel
2. Allows specification of dislocation positions in the unit cell
3. Handles elastic constants
4. Supports MACE potential for structure relaxation
5. Outputs relaxed dislocation structures and stress/energy data

The generated structures follow the S configuration pattern typical for dislocation studies.

Example usage: python ESPM.py --reference_cell POSCAR --output_dir ./output --n_cells 8 8 --potential_path /path/to/mace_model.pt --xpos 0 0 10 -10 5 -5 --ypos 5 -5 0 0 5 -5 --babel_path /path/to/babel --fmax 0.0005 --meshing 20 30 --remove_original false --cij 177.1 84.8 82.9 193.8 54.8
"""

script_dir = os.path.dirname(os.path.abspath(__file__))
babel_path = get_tool_path('babel')

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="ESPM dislocation structure and relaxation script.")
    parser.add_argument("--reference_cell", required=True, type=str, help="Reference cell structure file.")
    parser.add_argument("--output_dir", required=True, type=str, help="Output directory.")
    parser.add_argument("--n_cells", nargs=2, type=int, default=[8, 8],
                       help='Number of unit cells in x- and y-direction (default: 8 8)')
    parser.add_argument("--potential_path", required=True, type=str, help="Path to the MACE model.")
    parser.add_argument("--xpos", nargs="+", type=int, default=[0, 0, 10, -10, 5, -5],
        help="List of x positions, e.g. 0 0 10 -10 5 -5 (default: 0 0 10 -10 5 -5)"
    )
    parser.add_argument("--ypos", nargs="+", type=int, default=[5, -5, 0, 0, 5, -5],
        help="List of y positions, e.g. 5 -5 0 0 5 -5 (default: 5 -5 0 0 5 -5)"
    )
    parser.add_argument('--fmax', type=float, default=0.0005, help="Force convergence criterion (default: 0.0005)")
    parser.add_argument('--meshing', nargs=2, type=int, default=[20, 30], help="Meshing for x and y directions (default: 20 30)")
    parser.add_argument('--remove_original', type=str, default="false", help="Remove original structure (default: false)")
    parser.add_argument('--cij', nargs=5, type=float, default=[177.1, 84.8, 82.9, 193.8, 54.8], help="Cij matrix [c11 c12 c13 c33 c44] (default: 177.1 84.8 82.9 193.8 54.8)")
    parser.add_argument('--analyze_core', type=str, default="false", help="Run analyze_core.py (default: false)")
    parser.add_argument('--device', type=str, default='cpu', help="Device for MACE calculation (default: cpu)")
    args = parser.parse_args()

    # Parse arguments
    if len(args.xpos) != len(args.ypos):
        raise ValueError("--xpos and --ypos must have the same number of elements")

    # Read the reference cell structure
    atoms = read(args.reference_cell)
    if isinstance(atoms, list):
        atoms = atoms[0]

    # Extract cell dimensions and lattice parameter
    c_tot, y_tot, a = atoms.cell.lengths()
    c, y = c_tot / args.n_cells[0], y_tot / args.n_cells[1]

    # Find the atom closest to the bottom left (min x + y) using atoms.positions
    positions = atoms.positions  # shape (N, 3)
    min_idx = np.argmin(positions[:, 0] + positions[:, 1])
    r01, r02 = positions[min_idx, 0], positions[min_idx, 1]
    print(f"r01 : {r01}, r02 : {r02}")

    # Calculate reference positions for dislocation insertion (S configuration)
    r11, r12 = 0.25 * c_tot + r01, 0.25 * y_tot + r02
    r21, r22 = 0.75 * c_tot + r01, 0.75 * y_tot + r02

    # Generate new POSCAR files for each dislocation configuration
    for i in range(len(args.xpos)):
        r11_i = r11 + c / 2 * args.xpos[i] / args.meshing[0]
        r12_i = r12 + y / 2 * args.ypos[i] / args.meshing[1]
        r21_i = r21 + c / 2 * args.xpos[i] / args.meshing[0]
        r22_i = r22 + y / 2 * args.ypos[i] / args.meshing[1]
        outfile = os.path.join(args.output_dir, f"x{args.xpos[i]}_y{args.ypos[i]}_nx{args.n_cells[0]}ny{args.n_cells[1]}.poscar")
        # Read and modify the input template for Babel
        with open(os.path.join(script_dir, "input.babel_S"), 'r') as f:
            content = f.read()
        # Replace tokens in the template with calculated values
        content = content.replace("A0", str(a))
        content = content.replace("R11", str(r11_i))
        content = content.replace("R21", str(r21_i))
        content = content.replace("R12", str(r12_i))
        content = content.replace("R22", str(r22_i))
        content = content.replace("REFFILE", args.reference_cell)
        content = content.replace("OUTFILE", outfile)
        content = content.replace("C11", str(args.cij[0]))
        content = content.replace("C12", str(args.cij[1]))
        content = content.replace("C13", str(args.cij[2]))
        content = content.replace("C33", str(args.cij[3]))
        content = content.replace("C44", str(args.cij[4]))
        content = content.replace("C66", str((args.cij[0] - args.cij[1]) / 2))
        # Write the modified template to a temporary file and run Babel
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp:
            tmp.write(content)
            tmpbabel = tmp.name
        subprocess.run([args.babel_path, tmpbabel], check=True)
        os.remove(tmpbabel)

    # Set up the MACE calculator for relaxation
    for i in range(len(args.xpos)):
        cell = os.path.join(args.output_dir, f"x{args.xpos[i]}_y{args.ypos[i]}_nx{args.n_cells[0]}ny{args.n_cells[1]}.poscar")
        # Read the generated structure
        atoms = read(cell, format="vasp")
        set_calculator(atoms, args.potential_path, device=args.device)
        # Relax the structure using BFGS
        dyn = BFGS(atoms)
        dyn.run(fmax=args.fmax)
        # Write the relaxed structure and results
        write(os.path.join(args.output_dir, f"x{args.xpos[i]}_y{args.ypos[i]}_nx{args.n_cells[0]}ny{args.n_cells[1]}_relax.poscar"), atoms, format="vasp", sort=True)

        if args.remove_original == "true" and args.analyze_core == "false":
            os.remove(cell)
    # Run analyze_core.py on all relaxed structures if requested
    if args.analyze_core == "true":
        # Build list of relaxed structures and output files
        relaxed_cells = [os.path.join(args.output_dir, f"x{x}_y{y}_nx{args.n_cells[0]}ny{args.n_cells[1]}_relax.poscar") 
                        for x,y in zip(args.xpos, args.ypos)]
        ref_dis_cells = [os.path.join(args.output_dir, f"x{x}_y{y}_nx{args.n_cells[0]}ny{args.n_cells[1]}.poscar") 
                        for x,y in zip(args.xpos, args.ypos)]
        output_files = [os.path.join(args.output_dir, f"x{x}_y{y}_nx{args.n_cells[0]}ny{args.n_cells[1]}.dat")
                       for x,y in zip(args.xpos, args.ypos)]
        
        # Build command to run analyze_core.py
        analyze_cmd = [
            sys.executable,
            "-m",
            "dislocate.core_analysis.analyze_core",
            "--ref_cell", args.reference_cell,
            "--dis_cells", " ".join(relaxed_cells),
            "--ref_dis_cells", " ".join(ref_dis_cells),
            "--output_files", " ".join(output_files), 
            "--potential_path", args.potential_path,
            "--potential_type", "MACE",
            "--nx", str(args.n_cells[0]),
            "--ncore", str(os.cpu_count()),
            "--thickness", "1",
            "--fitting", "true",
            "--ovito", "true",
            "--nye", "true",
        ]
        
        # Run the analysis
        subprocess.run(analyze_cmd, check=True)

        if args.remove_original == "true":
            for cell in ref_dis_cells:
                os.remove(cell)



# Run the main function if this script is executed
if __name__ == "__main__":
    main()
