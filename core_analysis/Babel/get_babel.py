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
        
        # Create Babel input file
        with open(abspath_from_script('input_nye.displacement'), 'r') as f:
            content = f.read()

        content = content.replace("A0", str(args.a0))
        content = content.replace("PATTERN", os.path.relpath(args.file_pattern, os.getcwd()))
        content = content.replace("DUPLICATE", str(duplicate).lower())
        content = content.replace("IMM", str(imm))
        content = content.replace("NREP", str(nrep))
        content = content.replace("REFFILE", os.path.relpath(args.ref_cell, os.getcwd()))
        content = content.replace("INFILE", os.path.relpath(args.dis_cell, os.getcwd()))
        content = content.replace("OUTFILE", os.path.relpath(args.tmp_babel, os.getcwd()))

        # Write the modified template to a temporary file and run Babel
        with tempfile.NamedTemporaryFile(prefix='babel-', directory='tmp', mode='w+', delete=False) as tmp:
            tmp.write(content)
            tmpbabel = tmp.name
        
        # Run Babel displacement calculation
        subprocess.run([args.babel_path, tmpbabel], check=True)
        
    finally:
        # Clean up temporary files
        if os.path.exists(tmpbabel):
            os.unlink(tmpbabel)

if __name__ == '__main__':
    main() 