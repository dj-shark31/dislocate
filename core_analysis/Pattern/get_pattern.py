#!/usr/bin/env python3
"""
Python version of get_pattern.sh: pattern detection for POSCAR files only.
"""
import argparse
import subprocess
import os
import tempfile

def escape_path(path):
    """Escape forward slashes in path for sed replacement"""
    return path.replace('/', '\\/')

def run(cmd, **kwargs):
    print('Running:', ' '.join(str(x) for x in cmd))
    return subprocess.run(cmd, check=True, **kwargs)

def abspath_from_script(rel_path):
    """Return absolute path relative to the directory containing this script."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, rel_path))

def main():
    parser = argparse.ArgumentParser(description='Python version of get_pattern.sh')
    parser.add_argument('ref_cell')
    parser.add_argument('thickness', type=int)
    parser.add_argument('a0')
    parser.add_argument('natom', type=int)
    parser.add_argument('tmp_pattern')
    parser.add_argument('--pattern_dat', default=abspath_from_script('pattern.dat'))
    parser.add_argument('--patternDetect_bin', default='/global/home/users/djany/Babel_V10.8/bin/patternDetect')
    args = parser.parse_args()

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
    with open(args.pattern_dat, 'r') as f:
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
    subprocess.run([args.patternDetect_bin, tmp_input], check=True)

    # Clean up
    os.remove(tmp_input)

if __name__ == '__main__':
    main() 