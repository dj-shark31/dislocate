#!/usr/bin/env python3
"""
Python version of get_ovito.sh
Runs the OVITO and fitting pipeline using Python subprocess calls.
Accounts for the use of fitting_core.py and updated argument names.
"""
import argparse
import subprocess
import os
import sys

def run(cmd, check=True):
    print('Running:', ' '.join(str(x) for x in cmd))
    result = subprocess.run(cmd, check=check)
    return result

def abspath_from_script(rel_path):
    """Return absolute path relative to the directory containing this script."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, rel_path))

def main():
    parser = argparse.ArgumentParser(description='Python version of get_ovito.sh')
    parser.add_argument('dis_cell')
    parser.add_argument('ref_cell')
    parser.add_argument('b', type=float, help='Burgers vector magnitude (b)')
    parser.add_argument('thickness', type=float, help='Cell thickness (thickness)')
    parser.add_argument('tmp_stab')
    parser.add_argument('tmp_dxa')
    parser.add_argument('tmp_fitting')
    parser.add_argument('--fitting', type=str, default='true')
    parser.add_argument('--oxygen', type=int, default=0)
    parser.add_argument('--pbc', type=str, default='false')
    parser.add_argument('--config', default='S')
    parser.add_argument('--ovito_bin', default='ovitos')
    parser.add_argument('--ovito_sif', default='/global/scratch/users/djany/ovito.sif')
    args = parser.parse_args()

    # Run ovito_elastStab.py
    run(['singularity', 'exec', '-B', '.:/mnt/', args.ovito_sif, args.ovito_bin,
         abspath_from_script("ovito_elastStab.py"),
         args.ref_cell, args.dis_cell, args.tmp_stab,
         str(args.b), str(args.oxygen)])

    # Run ovito_dxa.py
    run(['singularity', 'exec', '-B', '.:/mnt/', args.ovito_sif, args.ovito_bin,
         abspath_from_script("ovito_dxa.py"),
         args.dis_cell, str(int(args.thickness)), args.tmp_dxa, str(args.oxygen), args.config])

    # Run fitting_core.py if requested
    if args.fitting == 'true':
        fit_cmd = [sys.executable, abspath_from_script("fitting_core.py"),
                   args.tmp_stab, args.tmp_dxa, str(int(args.thickness)), args.tmp_fitting,
                   '--pbc', args.pbc]
        run(fit_cmd)

if __name__ == '__main__':
    main() 