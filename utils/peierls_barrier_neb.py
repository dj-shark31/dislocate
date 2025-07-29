from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE, BFGS
from utils.atomistic_tools import set_calculator
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

def relax_intermediate_images(images: List[Atoms], potential_path: str, relax_fmax: float = 0.001, output_dir: str = None, device: str = "cpu") -> None:
    """Relax all intermediate images individually and output their energies."""
    
    print(f"Relaxing {len(images)} intermediate images individually...")
    
    relaxed_energies = []
    relaxed_images = []
    for i, image in enumerate(images):
        print(f"  Relaxing image {i}...")

        set_calculator(image, potential_path, device=device)
        
        # Create optimizer for this image
        optimizer = BFGS(image, logfile=None)
            
        # Relax the image
        optimizer.run(fmax=relax_fmax, steps=1000)
            
        # Get final energy and steps
        relaxed_energies.append(image.get_potential_energy())
        relaxed_images.append(image)
        print(f"Image {i} relaxed: Energy = {relaxed_energies[-1]:.6f} eV")


    if output_dir is not None:
        relaxed_energies_file = os.path.join(output_dir, 'energies_relaxed_cell.dat')
        with open(relaxed_energies_file, 'w') as f:
            f.write("# Relaxed Intermediate Image Energies\n")
            f.write("# Image_Index  Energy_eV \n")
            for i, energy in enumerate(relaxed_energies):
                f.write(f"{i:10d}  {energy:12.6f} \n")
                write(os.path.join(output_dir, f'{i}_relaxed.poscar'), relaxed_images[i], format='vasp')
        print(f"Relaxed intermediate energies written to {relaxed_energies_file}")

    return relaxed_energies, relaxed_images
    

def run_neb(initial: Atoms, final: Atoms, n_images: int, potential_path: str, 
            output_dir: str = None, neb_fmax: float = 0.005, output_preopt: str = "false", device: str = "cpu") -> List[Atoms]:
    """
    Run NEB calculation between initial and final images and save outputs.
    
    Args:
        initial: Initial atomic configuration
        final: Final atomic configuration 
        n_images: Number of intermediate images
        potential_path: Path to MACE potential file
        device: Device for MACE calculation (cpu/cuda)
        output_dir: Directory to save output files
        neb_fmax: Force convergence criterion for NEB
        output_preopt: Whether to output pre-optimized images
        
    Returns:
        List of optimized images including initial and final
    """
    # Create images list: initial + intermediate + final
    images = [initial]
    for _ in range(n_images):
        images.append(initial.copy())
    images.append(final)
    
    # Interpolate linearly between initial and final
    neb = NEB(images)
    neb.interpolate(mic=True)
    
    # Assign MACE calculator to all images
    print(f"Loading MACE potential from: {potential_path}")
    for image in images:
        set_calculator(image, potential_path, device=device)

    # Output pre-optimized images if requested
    if output_preopt == "true":
        print(f"Outputting pre-optimized images to {output_dir}")
        for i, image in enumerate(images):
            write(os.path.join(output_dir, f'{i}_preopt.poscar'), image, format='vasp')

    # Run NEB optimization
    print(f"Starting NEB optimization with {n_images} intermediate images...")
    optimizer = FIRE(neb, trajectory=os.path.join(output_dir,'neb.traj'), 
                    logfile=os.path.join(output_dir,'neb.log'))
    optimizer.run(fmax=neb_fmax, steps=1000)

    if output_dir is not None:
        # Write final energies
        write_energies(images, os.path.join(output_dir, 'energies.dat'))

        # Save optimized images
        for i, image in enumerate(images):
            write(os.path.join(output_dir, f'{i}.poscar'), image, format='vasp')

    print("NEB optimization completed")
    return images


def main():
    # Parse command line arguments
    args = parse_arguments()

    # Create directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.perform_neb == "true":
        # Load initial and final images
        print(f"Loading initial image from: {args.initial}")
        print(f"Loading final image from: {args.final}")
        print(f"Using MACE potential: {args.potential_path}")
        
        initial = read(args.initial)
        final = read(args.final)

        # Run NEB optimization
        images = run_neb(initial, final, args.n_images, args.potential_path, args.output_dir, 
                                        args.neb_fmax, args.output_preopt_images, args.device)

    else:
        print("Skipping NEB optimization")
        
    # Relax intermediate images if requested
    if args.relax_intermediate_images == "true":
        if args.perform_neb == "false":
            print(f"Loading images from {args.output_dir}")
            images = []
            for i in range(1, args.n_images + 1):
                images.append(read(os.path.join(args.output_dir, f'{i}.poscar')))
        else:
            images = images[1:-1] # Remove initial and final images
        print(f"Relaxing intermediate images to {args.output_dir}")
        relax_intermediate_images(images, args.potential_path, args.relax_fmax, args.output_dir, args.device)
        
    print("NEB simulation complete. Results saved in 'neb_images'.")

if __name__ == "__main__":
    main()
