from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE, BFGS
from mace.calculators import MACECalculator
from ase import Atoms
import os
import argparse
from typing import List

"""
Script to calculate the Peierls barrier using the Nudged Elastic Band (NEB) method.

This script performs NEB calculations between initial and final dislocation configurations
to find the minimum energy path for dislocation motion. It uses the MACE potential for 
energy/force calculations and supports:

- NEB optimization with customizable number of images
- Individual image relaxation after NEB optimization
- Output of energies and atomic configurations at each step
- Flexible device selection (CPU/GPU) for MACE calculations
- Customizable convergence criteria for both NEB and individual relaxations

Example command line: python peierls_barrier_neb.py --initial POSCAR_initial --final POSCAR_final --potential_path /path/to/mace_model.pt --n_images 5 --output_dir neb_results --neb_fmax 0.005 --relax_fmax 0.001 --relax_intermediate_images true --output_preopt_images false --perform_neb true --device cpu
"""

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run Peierls barrier NEB calculation with MACE potential')
    parser.add_argument('--initial', type=str, required=True,
                       help='Name of the initial image file (e.g., POSCAR_initial)')
    parser.add_argument('--final', type=str, required=True,
                       help='Name of the final image file (e.g., POSCAR_final)')
    parser.add_argument('--potential_path', type=str, required=True,
                       help='Path to MACE potential file')
    parser.add_argument('--n_images', type=int, default=5,
                       help='Number of intermediate images (default: 5)')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory for energies and images ({output_dir}/energies.dat) and {output_dir}/{number}.poscar')
    parser.add_argument('--neb_fmax', type=float, default=0.005,
                       help='Force convergence criterion for NEB optimization (default: 0.005)')
    parser.add_argument('--relax_fmax', type=float, default=0.001,
                       help='Force convergence criterion for individual image relaxation (default: 0.001)')
    parser.add_argument('--relax_intermediate_images', default="true",
                       help='Relax all intermediate images individually after NEB optimization (default: true)')
    parser.add_argument('--output_preopt_images', default="false",
                       help='Output pre-optimized images (default: false)')
    parser.add_argument('--perform_neb', default="true",
                       help='Perform NEB optimization (default: true)')
    parser.add_argument('--device', type=str, default='cpu',
                       choices=['cpu', 'cuda'], help='Device for MACE calculation (default: cpu)')
    return parser.parse_args()

def write_energies(images: List[Atoms], output_file: str) -> None:
    """Write energies of all images to output file."""
    with open(output_file, 'w') as f:
        f.write("# NEB Image Energies\n")
        f.write("# Image_Index  Energy_eV\n")
        
        for i, image in enumerate(images):
            energy = image.get_potential_energy()
            f.write(f"{i:10d}  {energy:12.6f}\n")
    
    print(f"Energies written to {output_file}")

def relax_intermediate_images(args, output_dir: str) -> None:
    """Relax all intermediate images individually and output their energies."""
    
    print(f"Relaxing {args.n_images} intermediate images individually...")
    
    # Create output file for relaxed energies
    relaxed_energies_file = os.path.join(output_dir, 'energies_relaxed_cell.dat')
    
    with open(relaxed_energies_file, 'w') as f:
        f.write("# Relaxed Intermediate Image Energies\n")
        f.write("# Image_Index  Energy_eV \n")
        
        # Relax intermediate images (skip first and last)
        for i in range(1, args.n_images + 1):
            print(f"  Relaxing image {i}...")

            intermediate_image = read(os.path.join(output_dir, f'{i}.poscar'))

            intermediate_image.calc = MACECalculator(model_paths=args.potential_path, device=args.device)
            
            # Create optimizer for this image
            optimizer = BFGS(intermediate_image, logfile=None)
            
            # Relax the image
            optimizer.run(fmax=args.relax_fmax, steps=1000)
            
            # Get final energy and steps
            energy = intermediate_image.get_potential_energy()
            
            # Write to file
            f.write(f"{i:10d}  {energy:12.6f} \n")
            
            # Save relaxed structure
            write(os.path.join(output_dir, f'{i}_relaxed.poscar'), intermediate_image, format='vasp')
            
            print(f"Image {i} relaxed: Energy = {energy:.6f} eV")
    
    print(f"Relaxed intermediate energies written to {relaxed_energies_file}")


def main():
    # Parse command line arguments
    args = parse_arguments()

    # Create directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.perform_neb == "true":
        # --- LOAD INITIAL AND FINAL IMAGES ---
        
        print(f"Loading initial image from: {args.initial}")
        print(f"Loading final image from: {args.final}")
        print(f"Using MACE potential: {args.potential_path}")
        
        initial = read(args.initial)
        final = read(args.final)
        
        # --- SETUP IMAGES ---
        
        # Create images list: initial + intermediate + final
        images = [initial]
        for _ in range(args.n_images):
            images.append(initial.copy())
        images.append(final)
        
        # Interpolate linearly between initial and final
        neb = NEB(images)
        neb.interpolate(mic=True)
        
        # --- ASSIGN MACE CALCULATOR ---
        
        print(f"Loading MACE potential from: {args.potential_path}")
        
        for image in images:
            image.calc = MACECalculator(model_paths=args.potential_path, device=args.device)

        # --- OUTPUT IMAGE INTERPOLATION ---

        if args.output_preopt_images == "true":
            print(f"Outputting pre-optimized images to {args.output_dir}")
            for i, image in enumerate(images):
                write(os.path.join(args.output_dir, f'{i}_preopt.poscar'), image, format='vasp')

        # --- OPTIMIZATION ---
    
        print(f"Starting NEB optimization with {args.n_images} intermediate images...")
        optimizer = FIRE(neb, trajectory=os.path.join(args.output_dir,'neb.traj'), logfile=os.path.join(args.output_dir,'neb.log'))
        optimizer.run(fmax=args.neb_fmax, steps=1000)
    
        # --- WRITE ENERGIES ---
        
        write_energies(images, os.path.join(args.output_dir, 'energies.dat'))

        # --- SAVE OUTPUT ---
        
        for i, image in enumerate(images):
            write(os.path.join(args.output_dir, f'{i}.poscar'), image, format='vasp')

    else:
        print("Skipping NEB optimization")
        
       # --- RELAX INTERMEDIATE IMAGES (OPTIONAL) ---
    if args.relax_intermediate_images == "true":
        print(f"Relaxing intermediate images to {args.output_dir}")
        relax_intermediate_images(args, args.output_dir)
        
    print("NEB simulation complete. Results saved in 'neb_images'.")

if __name__ == "__main__":
    main()
