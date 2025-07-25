#an ovito script to be run as follows 
#singularity exec --bind 
# .:/mnt ovitos /mnt/ovito_elastStab.py [infile] [reffile] [outfile]
#
# all files must be in folder this is called from
# Make sure to change burgers content and elastic stability first as appropriate
import argparse
from ovito.io import import_file, export_file
from ovito.modifiers import *

def main():
    parser = argparse.ArgumentParser(description='OVITO script for Elastic Stability Analysis')
    parser.add_argument('reffile', help='Reference filename')
    parser.add_argument('infile', help='Input filename')
    parser.add_argument('outfile', help='Output filename')
    parser.add_argument('b', type=float, help='Burgers vector magnitude')
    parser.add_argument('oxygen', type=int, choices=[0, 1], 
                       help='Delete oxygen or not (0= no oxygen, 1= delete oxygen)')
    
    args = parser.parse_args()
    
    #Import the file
    node = import_file(args.infile)
    nodeRef = import_file(args.reffile)
    
    #Get z-periodicity vector
    #nodeCell = import_file(infile)
    #dataCell = nodeCell.compute()
    #c = dataCell.cell[:,2]/2
    
    if args.oxygen == 1:
        # Select type:
        node.modifiers.append(SelectTypeModifier(types = {2}))
        # Delete selected:
        node.modifiers.append(DeleteSelectedModifier())
        
    #Add atomic strain modifier - will need to change burgers content appropriately
    modifier = AtomicStrainModBurgers(#burgersContent = [c[0], c[1], c[2]],
        burgersContent = [0, 0, args.b], #burgers content is with respect to the cartesian frame, our burgers vector is in z direction
        output_deformation_gradients = True,
            output_strain_tensors = True)
    
    modifier.reference.load(args.reffile)
    node.modifiers.append(modifier)
    
    #Add stability modifier - will need to change elastic constants appropriately
    #For Hex
    #input order c11, c12, c13, c33, c44
    #c111, c112, c113, c123, c133, c144, c155, c222, c333, c344
    stabmod = CalculateElasticStabilityModifier(
            set_structure = CalculateElasticStabilityModifier.Lattice.Hexagonal_High,
            set_soec = [174, 95, 72, 188, 58],
            set_toec = [-1850, -965, -443, -628, -292, -129, 170, -1597.25, -93.5, 117],
            set_transformation_matrix = [[0.0, 0.0, 1.0 , 0],
                                         [0.0, 1.0, 0, 0],
                                         [-1.0, 0, 0, 0]] #3x4 matrix, let the last column be 0 and the first three columns act as 3x3 rotation matrix from typical directions to actual directions
            )
    node.modifiers.append(stabmod)
    
    #Add displacement modifier, this is not mod anything
    #displacement = CalculateDisplacementsModifier()
    #displacement.reference.load(reffile)
    #node.modifiers.append(displacement)
    
    # Common neighbor analysis:
    node.modifiers.append(CommonNeighborAnalysisModifier(cutoff = 3.4))
    
    #Calculate
    data = node.compute()
    dataRef = nodeRef.compute()
    
    #Structure type
    cna = data.particles['Structure Type']
    
    #Origin of cell                                                                                                                                                                                                                      
    o = dataRef.cell[:,3]
    
    #Get positions in reference cell
    x = dataRef.particle_properties.position.array[:, 0]
    y = dataRef.particle_properties.position.array[:, 1]
    z = dataRef.particle_properties.position.array[:, 2]
    
    #Get positions in dislocation cell                                                         
    xd = data.particle_properties.position.array[:, 0]
    yd = data.particle_properties.position.array[:, 1]
    zd = data.particle_properties.position.array[:, 2]
    
    #Get stability parameter
    stab = data.particles['Stability Parameter']
    
    #Print out reference cell x, y, z - dislocation cell x, y, z - stab par - cna
    with open(args.outfile, 'w') as f:
        for i in range(0, data.number_of_particles): 
            f.write(f"{x[i] - o[0]:.4f} {y[i] - o[1]:.4f} {z[i] - o[2]:.4f} {xd[i] - o[0]:.4f} {yd[i] - o[1]:.4f} {zd[i] - o[2]:.4f} {stab[i]:.5f} {cna[i]:d}\n")

if __name__ == "__main__":
    main()

