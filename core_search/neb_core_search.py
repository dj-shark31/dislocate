#!/usr/bin/env python3
"""
NEB Core Search Script

This script performs iterative NEB (Nudged Elastic Band) calculations between two dislocation 
dipole cells to find the minimum energy path. It iteratively refines the path until all 
intermediate images have the same energy as the initial and final images.

The script:
1. Takes two dislocation dipole cells as input
2. Performs NEB calculation using peierls_barrier_neb.py
3. Relaxes intermediate images
4. Identifies images with different energy than endpoints
5. Repeats the cycle until all images have the same energy

Usage:
    python neb_core_search.py --initial_file POSCAR_initial --final_file POSCAR_final --potential_path /path/to/mace_model.pt --n_images 5 --output_dir neb_core_search --max_iterations 10 --energy_tolerance 0.001 --neb_fmax 0.005 --relax_fmax 0.001 --device cpu
"""

import argparse
import os
from ase.io import read, write
from utils.peierls_barrier_neb import run_neb, relax_intermediate_images
from ase.calculators.mace import MACECalculator

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Iterative NEB core search between dislocation dipole cells')
    parser.add_argument('--initial_file', type=str, required=True,
                       help='Initial dislocation dipole cell file (e.g., POSCAR_initial)')
    parser.add_argument('--final_file', type=str, required=True,
                       help='Final dislocation dipole cell file (e.g., POSCAR_final)')
    parser.add_argument('--potential_path', type=str, required=True,
                       help='Path to MACE potential file')
    parser.add_argument('--n_images', type=int, default=5,
                       help='Number of intermediate images for NEB (default: 5)')
    parser.add_argument('--output_dir', type=str, default='neb_core_search',
                       help='Output directory for results (default: neb_core_search)')
    parser.add_argument('--max_iterations', type=int, default=10,
                       help='Maximum number of iterations (default: 10)')
    parser.add_argument('--energy_tolerance', type=float, default=0.001,
                       help='Energy tolerance for convergence in eV (default: 0.001)')
    parser.add_argument('--neb_fmax', type=float, default=0.005,
                       help='Force convergence criterion for NEB optimization (default: 0.005)')
    parser.add_argument('--relax_fmax', type=float, default=0.001,
                       help='Force convergence criterion for individual image relaxation (default: 0.001)')
    parser.add_argument('--device', type=str, default='cpu',
                       choices=['cpu', 'cuda'], help='Device for MACE calculation (default: cpu)')
    return parser.parse_args()

def get_energy(image, potential_path, device):
    """
    Get the energy of an image.
    """
    image.calc = MACECalculator(model_paths=potential_path, device=device)
    return image.get_potential_energy()

def get_unique_state_energy(initial_energy, final_energy, intermediate_energies, energy_tolerance):
    unique_states = []
    for i, energy in enumerate(intermediate_energies):
        # Check if energy is significantly different from both endpoints
        if (abs(energy - initial_energy) > energy_tolerance and 
            abs(energy - final_energy) > energy_tolerance):
                
            # Check if this state is different from previously found states
            is_unique = True
            for existing_state in unique_states:
                if abs(energy - existing_state['energy']) < energy_tolerance:
                    is_unique = False
                    break
                        
            if is_unique:
                unique_states.append({'index': i, 'energy': energy})
    return unique_states

def get_new_cell(initial_image, final_image, n_images, potential_path, neb_fmax=0.005, relax_fmax=0.001, energy_tolerance=0.001, device='cpu'):
    """
    Get a new cell by running NEB and relaxing intermediate images.
    """
    initial_energy = get_energy(initial_image, potential_path, device)
    final_energy = get_energy(final_image, potential_path, device)
    if abs(initial_energy - final_energy) < energy_tolerance:
        print("Initial and final energies are too close to each other. No need to run NEB.")
        return None

    images = run_neb(initial_image, final_image, n_images, potential_path, neb_fmax=neb_fmax, device=device)
    intermediate_energies, intermediate_images = relax_intermediate_images(images[1:-1], potential_path, relax_fmax=relax_fmax, device=device)
    
    # Find intermediate states with different energies than endpoints
    different_energy_states = get_unique_state_energy(initial_energy, final_energy, intermediate_energies, energy_tolerance)
    print(f"Found {len(different_energy_states)} intermediate states with distinct energies:")
        
    if len(different_energy_states) == 0:
        print("No intermediate states with distinct energies found. Returning None.")
        return None
    else:
        for state in different_energy_states:
            state['image'] = intermediate_images[state['index']]
            print(f"Image {state['index']}: {state['energy']:.6f} eV")
        return different_energy_states

def iteration(current_initial, current_final, args, iteration):
    if iteration > args.max_iterations:
        return 0
    else:
        new_cells = get_new_cell(current_initial, current_final, args.n_images, args.potential_path, neb_fmax=args.neb_fmax, relax_fmax=args.relax_fmax, energy_tolerance=args.energy_tolerance, device=args.device)
        if new_cells is None:
            return 0
        else:
            for i, cell in enumerate(new_cells):
                write(os.path.join(args.output_dir, f'iteration_{iteration}_{i}.poscar'), cell['image'], format='vasp')
            
            for i in range(len(new_cells)):
                if i == 0:
                    cell_i = read(os.path.join(args.output_dir, f'iteration_{iteration}_{i}.poscar'))
                    iteration(current_initial, cell_i, args, iteration + 1)
                if i == len(new_cells) - 1:
                    cell_i = read(os.path.join(args.output_dir, f'iteration_{iteration}_{i}.poscar'))
                    iteration(cell_i, current_final, args, iteration + 1)
                if i != 0 and i != len(new_cells) - 1:
                    cell_i_minus_1 = read(os.path.join(args.output_dir, f'iteration_{iteration}_{i-1}.poscar'))
                    cell_i = read(os.path.join(args.output_dir, f'iteration_{iteration}_{i}.poscar'))
                    iteration(cell_i_minus_1, cell_i, args, iteration + 1)
    

def main():
    """Main function for iterative NEB core search."""
    args = parse_arguments()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize files
    initial = read(args.initial_file)
    final = read(args.final_file)
    
    print(f"Starting iterative NEB core search")
    print(f"Initial cell: {initial}")
    print(f"Final cell: {final}")
    print(f"Potential: {args.potential_path}")
    print(f"Energy tolerance: {args.energy_tolerance} eV")
    print(f"Max iterations: {args.max_iterations}")
    print(f"Number of images: {args.n_images}")
    print(f"Device: {args.device}")
    
    iteration(initial, final, args, 1)
    

if __name__ == "__main__":
    main() 