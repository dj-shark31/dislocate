#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 17:05:18 2021
@author: davidjany
"""

import argparse
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.data import DislocationNetwork
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='OVITO script for Dislocation Analysis (DXA)')
    parser.add_argument('infile', help='Input filename')
    parser.add_argument('thickness', type=int, help='Cell thickness')
    parser.add_argument('outfile', help='Output filename')
    parser.add_argument('oxygen', type=int, choices=[0, 1], 
                       help='Delete oxygen or not (0= no oxygen, 1= delete oxygen)')
    parser.add_argument('config', choices=['S', 'O'], 
                       help='Dislocation line configuration: S or O')
    
    args = parser.parse_args()
    
    node = import_file(args.infile)
    
    if args.oxygen == 1:
        # Select type:
        node.modifiers.append(SelectTypeModifier(types = {2}))
        # Delete selected:
        node.modifiers.append(DeleteSelectedModifier())
    
    # Replicate:
    if args.thickness < 2:
        node.modifiers.append(ReplicateModifier(num_z = 3))
    elif args.thickness == 2:
        node.modifiers.append(ReplicateModifier(num_z = 2))
    else:
        pass
    
    # Dislocation analysis (DXA):
    node.modifiers.append(DislocationAnalysisModifier(
        trial_circuit_length = 25, 
        input_crystal_structure = DislocationAnalysisModifier.Lattice.HCP, 
        only_perfect_dislocations = True))
    
    data = node.compute()
    j = 0
    for segment in data.dislocations.segments:
        if j==0:
            left=segment.points
            j+=1
        else:
            right=segment.points
                    
    leftPos=np.mean(left,axis=0)
    rightPos=np.mean(right,axis=0)
    
    with open(args.outfile, 'w') as f:
        if args.config == "O": # For O-configuration
            if int(leftPos[0]) < int(rightPos[0]):
                f.write(f"{leftPos[0]:.3f} {leftPos[1]:.3f} xDXA yDXA leftDis \n")
                f.write(f"{rightPos[0]:.3f} {rightPos[1]:.3f} xDXA yDXA rightDis \n")
            else:
                f.write(f"{rightPos[0]:.3f} {rightPos[1]:.3f} xDXA yDXA leftDis \n")
                f.write(f"{leftPos[0]:.3f} {leftPos[1]:.3f} xDXA yDXA rightDis \n")
        else:  #For S-configuration (left ->bottom dis, right -> top dis)
            if int(leftPos[1]) < int(rightPos[1]):
                f.write(f"{leftPos[0]:.3f} {leftPos[1]:.3f} xDXA yDXA leftDis \n")
                f.write(f"{rightPos[0]:.3f} {rightPos[1]:.3f} xDXA yDXA rightDis \n")
            else:
                f.write(f"{rightPos[0]:.3f} {rightPos[1]:.3f} xDXA yDXA leftDis \n")
                f.write(f"{leftPos[0]:.3f} {leftPos[1]:.3f} xDXA yDXA rightDis \n")

if __name__ == "__main__":
    main()
