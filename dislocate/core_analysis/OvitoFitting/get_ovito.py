#!/usr/bin/env python3
"""
Python version of get_ovito.sh
Runs the OVITO and fitting pipeline using Python subprocess calls.
Accounts for the use of fitting_core.py and updated argument names.
"""
import argparse
import sys
import os
import subprocess

# Add project root to Python path for subprocess compatibility
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from dislocate.core_analysis.analyze_core import run

def main():
    parser = argparse.ArgumentParser(description='Python version of get_ovito.sh')
    parser.add_argument('--dis_cell', type=str, required=True, help='Dislocation cell file')
    parser.add_argument('--ref_cell', type=str, required=True, help='Reference cell file')
    parser.add_argument('--b', type=float, required=True, help='Burgers vector magnitude (b)')
    parser.add_argument('--thickness', type=float, required=True, help='Cell thickness (thickness)')
    parser.add_argument('--tmp_stab', type=str, required=True, help='Temporary stability file')
    parser.add_argument('--tmp_dxa', type=str, required=True, help='Temporary dxa file')
    parser.add_argument('--tmp_fitting', type=str, required=True, help='Temporary fitting file')
    parser.add_argument('--fitting', type=str, default='true', help='Fitting flag')
    parser.add_argument('--oxygen', type=int, default=0, help='Oxygen flag')
    parser.add_argument('--pbc', type=str, default='false', help='PBC flag')
    parser.add_argument('--config', default='S', help='Configuration file')
    args = parser.parse_args()

    # Try to check if docker is available by running a simple docker command
    try:
        result = run(['docker', 'info'], capture_output=True)
        # If docker command succeeds, use docker
        docker_container = 'dj31/ovito-wallace:amd64'
        # Mount current directory and use absolute paths
        container_cmd = ['docker', 'run', '--rm',
                        '-v', f'{os.getcwd()}:/mnt',
                        docker_container]
        print("Using docker container")
    except (FileNotFoundError, subprocess.CalledProcessError):
        # If docker fails, fall back to singularity
        print("Using singularity container") 
        from utils.config_loader import get_tool_path
        singularity_container = get_tool_path('ovitosif')
        container_cmd = ['singularity', 'exec', '-B', f'{os.getcwd()}:/mnt', singularity_container]

    # Run ovito_elastStab.py
    run(container_cmd + ['ovitos',
         "/mnt/" + os.path.relpath(os.path.join(script_dir, "ovito_elastStab.py"), os.getcwd()),
         "--reffile", "/mnt/" + os.path.relpath(args.ref_cell, os.getcwd()), 
         "--infile", "/mnt/" + os.path.relpath(args.dis_cell, os.getcwd()), 
         "--outfile", "/mnt/" + os.path.relpath(args.tmp_stab, os.getcwd()),
         "--b", str(args.b), "--oxygen", str(args.oxygen)])

    # Run ovito_dxa.py
    run(container_cmd + ['ovitos',
     "/mnt/" + os.path.relpath(os.path.join(script_dir, "ovito_dxa.py"), os.getcwd()),
     "--infile", "/mnt/" + os.path.relpath(args.dis_cell, os.getcwd()), 
     "--thickness", str(int(args.thickness)), 
     "--outfile", "/mnt/" + os.path.relpath(args.tmp_dxa, os.getcwd()), 
     "--oxygen", str(args.oxygen), 
     "--config", args.config])

    # Run fitting_core.py if requested
    if args.fitting == 'true':
        fit_cmd = [sys.executable, os.path.join(script_dir,"fitting_core.py"),
                   "--outStab", os.path.relpath(args.tmp_stab, os.getcwd()),
                   "--outDXA", os.path.relpath(args.tmp_dxa, os.getcwd()),
                   "--thickness", str(int(args.thickness)),
                   "--outFitting", os.path.relpath(args.tmp_fitting, os.getcwd()),
                   '--pbc', args.pbc]
        run(fit_cmd)

if __name__ == '__main__':
    main() 