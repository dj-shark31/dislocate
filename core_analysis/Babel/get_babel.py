#!/usr/bin/env python3
"""
Python version of get_babel.sh
Handles Babel displacement calculations for dislocation analysis
"""

import argparse
import os
import subprocess
import tempfile
from analyze_core import abspath_from_script

def create_temp_file(prefix='babel-', suffix='', directory='tmp'):
    """Create a temporary file with the specified prefix and suffix"""
    if not os.path.exists(directory):
        os.makedirs(directory)
    return tempfile.NamedTemporaryFile(prefix=prefix, suffix=suffix, 
                                      dir=directory, delete=False)

def create_babel_input(template_file, output_file, params):
    """Create Babel input file by substituting parameters in template"""
    with open(template_file, 'r') as f:
        content = f.read()
    
    # Substitute parameters
    for key, value in params.items():
        content = content.replace(key, str(value))
    
    with open(output_file, 'w') as f:
        f.write(content)

def run_babel_displacement(input_file, babel_path):
    """Run Babel displacement calculation"""
    subprocess.run([babel_path, input_file], check=True)

def main():
    parser = argparse.ArgumentParser(description='Run Babel displacement calculations')
    parser.add_argument('dis_cell', help='Dislocation cell file')
    parser.add_argument('ref_cell', help='Reference cell file')
    parser.add_argument('thickness', type=int, help='Thickness of the slab')
    parser.add_argument('a0', type=float, help='Lattice parameter')
    parser.add_argument('natom', type=int, help='Number of atoms')
    parser.add_argument('file_pattern', help='Pattern file')
    parser.add_argument('tmp_babel', help='Temporary Babel output file')
    parser.add_argument('oxygen', type=int, choices=[0, 1], 
                       help='Oxygen flag (0=keep, 1=delete)')
    parser.add_argument('--babel_path', default=abspath_from_script('../../bin/Babel_V10.8/bin/displacement'))
    
    args = parser.parse_args()
    
    # Create temporary files
    tmp_input = create_temp_file(prefix='babel-', directory='tmp')
    
    try:
        # Determine replication parameters based on b value
        if args.thickness == 1:
            nrep = 3
            duplicate = True
        elif args.thickness == 2:
            nrep = 2
            duplicate = True
        else:
            nrep = 1
            duplicate = False
        
        imm = args.natom * nrep
        
        # Create parameter dictionary for template substitution
        params = {
            'A0': args.a0,
            'PATTERN': os.path.relpath(args.file_pattern, os.getcwd()),
            'DUPLICATE': str(duplicate).lower(),
            'IMM': imm,
            'NREP': nrep,
            'REFFILE': os.path.relpath(args.ref_cell, os.getcwd()),
            'INFILE': os.path.relpath(args.dis_cell, os.getcwd()),
            'OUTFILE': os.path.relpath(args.tmp_babel, os.getcwd())
        }
        
        # Create Babel input file
        template_file = abspath_from_script('input_nye.displacement')
        create_babel_input(template_file, tmp_input.name, params)
        
        # Run Babel displacement calculation
        run_babel_displacement(tmp_input.name, args.babel_path)
        
    finally:
        # Clean up temporary files
        if os.path.exists(tmp_input.name):
            os.unlink(tmp_input.name)

if __name__ == '__main__':
    main() 