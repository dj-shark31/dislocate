#an ovito script to be run as follows 
#singularity exec --bind 
# .:/mnt ovitos /mnt/ovito_elastStab.py [infile] [reffile] [outfile]
#
# all files must be in folder this is called from
# Make sure to change burgers content and elastic stability first as appropriate
import sys
import argparse
from ovito.io import import_file, export_file
from ovito.modifiers import *
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='OVITO script for Common Neighbor Analysis (CNA)')
    parser.add_argument('infile', help='Input filename')
    parser.add_argument('outfile', help='Output filename')
    parser.add_argument('oxygen', type=int, choices=[0, 1], 
                       help='Delete oxygen or not (0= no oxygen, 1= delete oxygen)')
    parser.add_argument('cell_poscar', choices=['true', 'false'], 
                       help='Whether using cell POSCAR format')
    
    args = parser.parse_args()
    
    #Import the file
    node = import_file(args.infile)
    
    if args.oxygen == 1:
        # Select type:
        node.modifiers.append(SelectTypeModifier(types = {2}))
        # Delete selected:
        node.modifiers.append(DeleteSelectedModifier())
    
    # Common neighbor analysis:
    node.modifiers.append(CommonNeighborAnalysisModifier(cutoff = 3.4))
    
    #Calculate
    data = node.compute()
    
    #Get atoms id
    if args.cell_poscar == "false":
        id = data.particles['Particle Identifier'][:]
    else:
        id = np.arange(data.number_of_particles)
    
    order = id.argsort()
    
    #Structure type
    cna = data.particles['Structure Type'][order]
    
    #Print out cna
    with open(args.outfile, 'w') as f:
        for i in range(0, data.number_of_particles): 
            f.write(f"{cna[i]:d}\n")

if __name__ == "__main__":
    main()

