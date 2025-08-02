#!/usr/bin/env python3
"""
Python version of get_patternInit.sh
Fills in patternInit.dat template, runs patternInit, and cleans up.
"""
import argparse
import tempfile
import subprocess
import os
import sys

# Add project root to Python path for subprocess compatibility
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from utils.config_loader import get_tool_path

def escape_path(path):
    """Escape forward slashes in path for sed replacement"""
    return path.replace('/', '\\/')

def main():
    parser = argparse.ArgumentParser(description='Python version of get_patternInit.sh')
    parser.add_argument('a0')
    parser.add_argument('coa0')
    parser.add_argument('tmp_pattern')
    args = parser.parse_args()

    patternInit_path = get_tool_path('patternInit')

    # Read template
    with open(os.path.join(script_dir, 'patternInit.dat'), 'r') as f:
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
    subprocess.run([patternInit_path, tmp_input], check=True)

    # Clean up
    os.remove(tmp_input)

if __name__ == '__main__':
    main() 