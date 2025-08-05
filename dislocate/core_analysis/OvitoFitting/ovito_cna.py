#an ovito script to be run as follows 
#singularity exec --bind 
# .:/mnt ovitos /mnt/ovito_elastStab.py [infile] [reffile] [outfile]
#
# all files must be in folder this is called from
# Make sure to change burgers content and elastic stability first as appropriate
import argparse
from ovito.io import import_file    
from ovito.modifiers import CommonNeighborAnalysisModifier, SelectTypeModifier, DeleteSelectedModifier
import numpy as np

def get_cna(infile, oxygen):
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
    
    structure_types = np.array(list(data.particles['Structure Type']))
                            
    total_atoms = len(structure_types)
    bcc_count = np.count_nonzero(structure_types == CommonNeighborAnalysisModifier.Type.BCC)
    fcc_count = np.count_nonzero(structure_types == CommonNeighborAnalysisModifier.Type.FCC)
    hcp_count = np.count_nonzero(structure_types == CommonNeighborAnalysisModifier.Type.HCP)
    ico_count = np.count_nonzero(structure_types == CommonNeighborAnalysisModifier.Type.ICO)
    other_count = total_atoms - (bcc_count + fcc_count + hcp_count + ico_count)
    
    return {'bcc': bcc_count, 'fcc': fcc_count, 'hcp': hcp_count, 'ico': ico_count, 'other': other_count}

def main():
    parser = argparse.ArgumentParser(description='OVITO script for Common Neighbor Analysis (CNA)')
    parser.add_argument('infile', help='Input filename')
    parser.add_argument('outfile', help='Output filename')
    parser.add_argument('oxygen', type=int, choices=[0, 1], 
                       help='Delete oxygen or not (0= no oxygen, 1= delete oxygen)')
    
    args = parser.parse_args()
    
    cna = get_cna(args.infile, args.oxygen)

    #Print out cna
    with open(args.outfile, 'w') as f:
        f.write(f"bcc fcc hcp ico other\n")
        f.write(f"{cna['bcc']} {cna['fcc']} {cna['hcp']} {cna['ico']} {cna['other']}\n")

if __name__ == "__main__":
    main()

