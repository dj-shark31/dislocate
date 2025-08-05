#!/usr/bin/env python3
"""
Python version of get_pattern.sh: pattern detection for POSCAR files only.
"""
import argparse
import os
import sys
import tempfile
import subprocess

# Add project root to Python path for subprocess compatibility
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from dislocate.utils.config_loader import get_tool_path

def escape_path(path):
    """Escape forward slashes in path for sed replacement"""
    return path.replace('/', '\\/')

def main():
    parser = argparse.ArgumentParser(description='Python version of get_pattern.sh')
    parser.add_argument('ref_cell')
    parser.add_argument('thickness', type=int)
    parser.add_argument('a0')
    parser.add_argument('natom', type=int)
    parser.add_argument('tmp_pattern')
    args = parser.parse_args()

    patternDetect_path = get_tool_path('patternDetect')

    # Determine nrep, duplicate, imm
    if args.thickness == 1:
        nrep = 3
        duplicate = 'true'
    elif args.thickness == 2:
        nrep = 2
        duplicate = 'true'
    else:
        nrep = 1
        duplicate = 'false'
    imm = args.natom * nrep

    # Read template
    with open(os.path.join(script_dir, 'pattern.dat'), 'r') as f:
        template = f.read()

    # Substitute values
    filled = template.replace('A0', str(args.a0)) \
        .replace('REFFILE', os.path.relpath(args.ref_cell, os.getcwd())) \
        .replace('OUTPATTERN', os.path.relpath(args.tmp_pattern, os.getcwd())) \
        .replace('DUPLICATE', duplicate) \
        .replace('IMM', str(imm)) \
        .replace('NREP', str(nrep))

    # Write to temp file
    with tempfile.NamedTemporaryFile('w', delete=False, prefix='pattern-', dir='tmp') as tmp:
        tmp_input = tmp.name
        tmp.write(filled)

    # Run patternDetect
    subprocess.run([patternDetect_path, tmp_input], check=True)

    # Clean up
    os.remove(tmp_input)

if __name__ == '__main__':
    main() 