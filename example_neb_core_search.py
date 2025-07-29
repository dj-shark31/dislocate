#!/usr/bin/env python3
"""
Example script demonstrating how to use neb_core_search.py

This script shows how to run iterative NEB calculations between two dislocation dipole cells
to find the minimum energy path for dislocation motion.

Prerequisites:
- Two POSCAR files representing initial and final dislocation configurations
- A trained MACE potential file
- The peierls_barrier_neb.py script in the utils/ directory
"""

import subprocess
import os
import sys

def run_neb_core_search_example():
    """Run an example NEB core search calculation."""
    
    # Example parameters - modify these for your specific case
    initial_poscar = "POSCAR_initial"  # Initial dislocation configuration
    final_poscar = "POSCAR_final"      # Final dislocation configuration
    potential_path = "potentials/mace_model.pt"  # Path to MACE potential
    output_dir = "neb_example_results"
    
    # Check if input files exist
    if not os.path.exists(initial_poscar):
        print(f"Error: Initial POSCAR file '{initial_poscar}' not found.")
        print("Please create or specify the correct path to your initial dislocation configuration.")
        return False
    
    if not os.path.exists(final_poscar):
        print(f"Error: Final POSCAR file '{final_poscar}' not found.")
        print("Please create or specify the correct path to your final dislocation configuration.")
        return False
    
    if not os.path.exists(potential_path):
        print(f"Error: MACE potential file '{potential_path}' not found.")
        print("Please specify the correct path to your trained MACE potential.")
        return False
    
    # Build the command
    cmd = [
        "python3", "neb_core_search.py",
        "--initial", initial_poscar,
        "--final", final_poscar,
        "--potential_path", potential_path,
        "--output_dir", output_dir,
        "--n_images", "5",           # Number of intermediate images
        "--max_iterations", "5",     # Maximum iterations
        "--energy_tolerance", "0.001", # Energy convergence tolerance (eV)
        "--neb_fmax", "0.005",       # NEB force convergence
        "--relax_fmax", "0.001",     # Individual relaxation force convergence
        "--device", "cpu",           # Use CPU (change to 'cuda' for GPU)
        "--refinement_strategy", "interpolated"  # Use interpolated refinement
    ]
    
    print("Running NEB core search with the following parameters:")
    print(f"  Initial configuration: {initial_poscar}")
    print(f"  Final configuration: {final_poscar}")
    print(f"  MACE potential: {potential_path}")
    print(f"  Output directory: {output_dir}")
    print(f"  Number of images: 5")
    print(f"  Max iterations: 5")
    print(f"  Energy tolerance: 0.001 eV")
    print(f"  Device: cpu")
    print(f"  Refinement strategy: interpolated")
    print()
    
    # Run the command
    try:
        result = subprocess.run(cmd, check=True, text=True)
        print("NEB core search completed successfully!")
        print(f"Results saved in: {output_dir}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running NEB core search: {e}")
        return False
    except FileNotFoundError:
        print("Error: neb_core_search.py not found in current directory.")
        print("Make sure you're running this script from the correct directory.")
        return False

def print_usage_instructions():
    """Print usage instructions for the NEB core search script."""
    print("""
NEB Core Search Usage Instructions
=================================

The neb_core_search.py script performs iterative NEB calculations to find the minimum 
energy path between two dislocation configurations.

Basic Usage:
    python neb_core_search.py --initial POSCAR_initial --final POSCAR_final --potential_path /path/to/mace_model.pt

Required Arguments:
    --initial          Initial dislocation configuration (POSCAR format)
    --final            Final dislocation configuration (POSCAR format)
    --potential_path   Path to trained MACE potential file

Optional Arguments:
    --n_images          Number of intermediate images (default: 5)
    --output_dir        Output directory (default: neb_core_search)
    --max_iterations    Maximum iterations (default: 10)
    --energy_tolerance  Energy convergence tolerance in eV (default: 0.001)
    --neb_fmax          NEB force convergence (default: 0.005)
    --relax_fmax        Individual relaxation force convergence (default: 0.001)
    --device            Device for calculations: cpu or cuda (default: cpu)
    --refinement_strategy  Refinement strategy: simple or interpolated (default: interpolated)

Example Commands:
    1. Basic run with default parameters:
       python neb_core_search.py --initial POSCAR_initial --final POSCAR_final --potential_path model.pt

    2. Custom parameters:
       python neb_core_search.py --initial POSCAR_initial --final POSCAR_final --potential_path model.pt --n_images 7 --max_iterations 15 --energy_tolerance 0.0005 --device cuda

    3. Simple refinement strategy:
       python neb_core_search.py --initial POSCAR_initial --final POSCAR_final --potential_path model.pt --refinement_strategy simple

Output:
    The script creates an output directory with subdirectories for each iteration.
    Final results are saved in 'final_results' subdirectory.
    Check 'energies_relaxed_cell.dat' files for energy profiles.
""")

if __name__ == "__main__":
    print("NEB Core Search Example")
    print("=======================")
    print()
    
    # Check if user wants to see usage instructions
    if len(sys.argv) > 1 and sys.argv[1] in ['--help', '-h', 'help']:
        print_usage_instructions()
    else:
        # Run the example
        success = run_neb_core_search_example()
        if not success:
            print("\nFor usage instructions, run: python example_neb_core_search.py --help") 