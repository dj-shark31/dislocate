#!/usr/bin/env python3
"""
Python version of get_patternInit.sh
Fills in patternInit.dat template, runs patternInit, and cleans up.
"""
import argparse
import tempfile
import subprocess
import os
from analyze_core import abspath_from_script
import subprocess

def escape_path(path):
    """Escape forward slashes in path for sed replacement"""
    return path.replace('/', '\\/')

def main():
    parser = argparse.ArgumentParser(description='Python version of get_patternInit.sh')
    parser.add_argument('a0')
    parser.add_argument('coa0')
    parser.add_argument('tmp_pattern')
    parser.add_argument('--patternInit_dat', default=abspath_from_script('patternInit.dat'))
    parser.add_argument('--patternInit_bin', default=abspath_from_script('../../bin/Babel_V10.8/bin/patternInit'))
    args = parser.parse_args()

    # Read template
    with open(args.patternInit_dat, 'r') as f:
        template = f.read()

    # Substitute values
    filled = template.replace('A0', str(args.a0)) \
        .replace('COA1', str(args.coa0)) \
        .replace('OUTPATTERN', os.path.relpath(args.tmp_pattern, os.getcwd()))

    # Write to temp file
    with tempfile.NamedTemporaryFile('w', delete=False, prefix='patternInit-', dir='tmp') as tmp:
        tmp_input = tmp.name
        tmp.write(filled)

    # Run patternInit
    subprocess.run([args.patternInit_bin, tmp_input], check=True)

    # Clean up
    os.remove(tmp_input)

if __name__ == '__main__':
    main() 