import argparse
from dislocate.utils.submitruns import find_compositions, submit_savio_job, composition_dir

"""
This script creates the batch files that are then submitted to the queue for processing

Example command line:
python submitruns_ESPM.py --elements Hf Ti Zr --delta_x 0.1 --potential_path /path/to/mace_model.pt --reference_cell POSCAR \
    --n_cells 8 8 --xpos 0 0 10 -10 5 -5 --ypos 5 -5 0 0 5 -5 \
    --babel_path /path/to/babel --fmax 0.0005 --meshing 20 30 \
    --remove_binaries false --remove_unaries true \
    --cij 177.1 84.8 82.9 193.8 54.8 \
    --job_name ESPM --partition savio4_htc --cpus 10 --time 5:00:00 --qos savio_lowprio --account co_chrzangroup \
    --analyze_core true

"""

def main():
    parser = argparse.ArgumentParser(description="Submit batch jobs for ESPM dislocation calculations.")
    parser.add_argument('--elements', required=True, nargs='+', help='Element symbols in the alloy (e.g. Hf Ti Zr)')
    parser.add_argument('--n_cells', nargs=2, type=int, default=[8, 8],
                       help='Number of unit cells in x- and y-direction (default: 8 8)')
    parser.add_argument('--delta_x', type=float, required=True, help='Composition increment (e.g. 0.1)')
    parser.add_argument('--potential_path', type=str, required=True, help='Path to the model to test')
    parser.add_argument('--reference_cell', type=str, required=True, help='Path to reference cell structure file')
    parser.add_argument('--xpos', nargs="+", type=int, default=[0, 0, 10, -10, 5, -5],
        help="List of x positions, e.g. 0 0 10 -10 5 -5 (default: 0 0 10 -10 5 -5)"
    )
    parser.add_argument('--ypos', nargs="+", type=int, default=[5, -5, 0, 0, 5, -5],
        help="List of y positions, e.g. 5 -5 0 0 5 -5 (default: 5 -5 0 0 5 -5)"
    )
    parser.add_argument('--remove_binaries', type=str, default="false", help='Remove binary compositions (true/false)')
    parser.add_argument('--remove_unaries', type=str, default="true", help='Remove unary compositions (true/false)')
    parser.add_argument('--babel_path', type=str, help='Path to Babel executable')
    parser.add_argument('--fmax', type=float, default=0.0005, help='Force convergence criterion')
    parser.add_argument('--meshing', nargs=2, type=int, default=[20, 30], help='Meshing for x and y directions')
    parser.add_argument('--remove_original', type=str, default="false", help='Remove original structure')
    parser.add_argument('--cij', nargs=5, type=float, default=[177.1, 84.8, 82.9, 193.8, 54.8], 
                       help='Cij matrix [c11 c12 c13 c33 c44]')
    parser.add_argument('--analyze_core', type=str, default="false", help="Run analyze_core.py (default: false)")
    parser.add_argument('--job_name', type=str, default="ESPM", help='SLURM job name')
    parser.add_argument('--partition', type=str, default="savio4_htc", help='SLURM partition')
    parser.add_argument('--cpus', type=int, default=10, help='Number of CPUs per task')
    parser.add_argument('--time', type=str, default="5:00:00", help='Wall time limit')
    parser.add_argument('--qos', type=str, default="savio_lowprio", help='Quality of service')
    parser.add_argument('--account', type=str, default="co_chrzangroup", help='SLURM account')
    args = parser.parse_args()

    # Find all compositions using the function from submitruns.py
    comps_to_run = find_compositions(args.elements, args.delta_x, 
                                   remove_binaries=args.remove_binaries, 
                                   remove_unaries=args.remove_unaries)

    print(f"There are {len(comps_to_run)} composition configurations.")

    jobs_submitted = 0
    for composition in comps_to_run:
        output_dir = composition_dir(args.elements, composition)
        
        # Prepare arguments for ESPM.py
        script_args = {
            "reference_cell": args.reference_cell,
            "output_dir": output_dir,
            "n_cells": f"{args.n_cells[0]} {args.n_cells[1]}",
            "potential_path": args.potential_path,
            "xpos": " ".join(map(str, args.xpos)),
            "ypos": " ".join(map(str, args.ypos)),
            "fmax": str(args.fmax),
            "meshing": f"{args.meshing[0]} {args.meshing[1]}",
            "remove_original": args.remove_original,
            "cij": " ".join(map(str, args.cij)),
            "analyze_core": args.analyze_core
        }
        
        # Add optional arguments if provided
        if args.babel_path:
            script_args["babel_path"] = args.babel_path
        
        # Submit job using the function from submitruns.py
        job_name = f"{args.job_name}_{output_dir.replace('/', '_').rstrip('_')}"
        submit_savio_job(
            script_name="ESPM.py",
            script_args_dict=script_args,
            job_name=job_name,
            partition=args.partition,
            cpus=args.cpus,
            time=args.time,
            qos=args.qos,
            account=args.account
        )
        jobs_submitted += 1

    print(f"The number of jobs submitted: {jobs_submitted}")

if __name__ == "__main__":
    main()
