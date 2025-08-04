from re import A
from ase.calculators.lammpslib import LAMMPSlib
from mace.calculators import MACECalculator
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from pymatgen.analysis.elasticity import diff_fit, Strain, DeformedStructureSet
import os
import sys
from subprocess import run
from ase.io import read

def set_lammpslib_calculator(potential_type, potential_path):
    """Setup LAMMPSlib calculator with the specified potential"""

    if potential_type == "MEAM":
        lmp_pair_style = 'meam/spline'
        lmp_pair_coeff = [f'* * {potential_path} Ti']
    elif potential_type == "DMD":
        lmp_pair_style = f'deepmd {potential_path}'
        lmp_pair_coeff = ['* *']
    elif potential_type == "RANN":
        lmp_pair_style = 'rann'
        lmp_pair_coeff = [f'* * {potential_path} Ti']
    elif potential_type == "ACE":
        lmp_pair_style = 'pace product'
        lmp_pair_coeff = [f'* * {potential_path} Ti']
    else:
        raise ValueError(f"Unsupported potential type: {potential_type}")
 
    lmp_cmds = [lmp_pair_style, lmp_pair_coeff]
    # LAMMPSlib calculator setup
    calc = LAMMPSlib(lmpcmds=lmp_cmds)
    return calc

def set_calculator(atoms, potential_path, potential_type='MACE', device='cpu'):
    if potential_type in ["MEAM", "DMD", "RANN", "ACE"]:
        calc = set_lammpslib_calculator(potential_type, potential_path)
    elif potential_type == 'MACE':
        calc = MACECalculator(model_paths=potential_path, device=device)
    else:
        raise ValueError(f"Unsupported potential type: {potential_type}")
    atoms.calc = calc

def get_energy(atoms, potential_path, potential_type='MACE', device='cpu'):
    """
    Get the energy of an image.
    """
    set_calculator(atoms, potential_path, potential_type, device)
    return atoms.get_potential_energy()

def get_stress(atoms, potential_path, potential_type='MACE', device='cpu'):
    """
    Get the energy of an image.
    """
    set_calculator(atoms, potential_path, potential_type, device)
    return atoms.get_stress(voigt=True) * 160217.66208 # Energy in eV, stress in MPa [xx, yy, zz, yz, xz, xy]

def get_unit_cell_vectors(atoms, n_cells):
    """Get lattice vectors from atoms object and number of unit cells."""
    vectors = np.array(atoms.cell)
    a1 = vectors[0] / n_cells[0]
    a2 = vectors[1] / n_cells[1]
    a3 = vectors[2] / n_cells[2]
    return a1, a2, a3

def alloy_sqs(atoms, elements, composition):
    """Alloy the atoms to the composition."""
    symbols = atoms.get_chemical_symbols()
    for k, symbol in enumerate(symbols):
        if len(elements) == 1 or composition[0] == 1:
            symbols[k] = elements[0]
        elif len(composition) > 1 and composition[1] == 1:
            symbols[k] = elements[1]
        elif len(composition) > 2 and composition[2] == 1:
            symbols[k] = elements[2]
        elif symbol == "A":
            symbols[k] = elements[0]
        elif symbol == "B":
            symbols[k] = elements[1]
        else:
            symbols[k] = elements[2]
    atoms.set_chemical_symbols(symbols)
    return atoms

def alloy_randomly(atoms, elements, composition):
    """Alloy atoms randomly to the composition."""
    n_atoms = len(atoms)
    atomic_distribution = [round(f * n_atoms) for f in composition[:-1]]
    atomic_distribution.append(n_atoms - sum(atomic_distribution))
    list_of_elements = [[elements[j]] * atomic_distribution[j] for j in range(len(atomic_distribution))]
    list_of_elements = [item for sublist in list_of_elements for item in sublist]
    atoms.set_chemical_symbols(list_of_elements)
    return atoms

def build_atoms(structure, n_cells):
    if structure == 'hcp':
        lat = Lattice.from_parameters(2.94, 2.94, 4.64, 90, 90, 120)
        cell = Structure(lat, ["Ti", "Ti"], [[0.0, 0.0, 0.0], [0.333333, 0.666667, 0.5]])
    elif structure == 'omega':
        lat = Lattice.from_parameters(4.57, 4.57, 2.83, 90, 90, 120)
        cell = Structure(lat, ["Ti", "Ti"], [[0.0, 0.0, 0.0], [0.333333, 0.666667, 0.5]], [0.666667, 0.333333, 0.5])
    elif structure == 'bcc':
        lat = Lattice.from_parameters(3.4, 3.4, 3.4, 90, 90, 90)
        cell = Structure(lat, ["Ti", "Ti"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    elif structure == 'fcc':
        lat = Lattice.from_parameters(4.3, 4.3, 4.3, 90, 90, 90)
        cell = Structure(lat, ["Ti"] * 4, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]])
    else:
        raise ValueError('Unknown structure type')
    atoms = AseAtomsAdaptor.get_atoms(cell.make_supercell(n_cells))
    atoms.pbc = (True, True, True)
    return atoms

def relax_cell(atoms, potential_path = None, potential_type = 'MACE', device = 'cpu', dof = [True, True, True, False, False, False], fmax = 0.004):
    """Relax the cell."""
    set_calculator(atoms, potential_path, potential_type, device)
    filter = FrechetCellFilter(atoms, dof)
    dyn = BFGS(filter)
    dyn.run(fmax=fmax)
    return atoms

def relax_atoms(atoms, potential_path = None, potential_type = 'MACE', device = 'cpu', fmax = 0.004):
    """Relax the atoms."""
    set_calculator(atoms, potential_path, potential_type, device)
    dyn = BFGS(atoms)
    dyn.run(fmax=fmax)
    return atoms

def get_lattice_parameters(atoms, n_cells = [1, 1, 1], potential_path = None, potential_type = 'MACE', device = 'cpu', dof = [True, True, True, False, False, False], fmax = 0.004, structure = 'hcp'):
    """Get the lattice parameters."""
    atoms = relax_cell(atoms, potential_path = potential_path, potential_type = potential_type, device = device, dof = dof, fmax = fmax)
    a1, a2, a3 = get_unit_cell_vectors(atoms, n_cells)
    if structure == 'hcp' or structure == 'omega':
        a = (a1 + a2) / 2
        c = a3
        return {'a': a, 'c': c}, atoms
    else:
        a = (a1 + a2 + a3) / 3
        return {'a': a}, atoms

def get_cij(atoms, potential_path = None, potential_type = 'MACE', device = 'cpu', fmax = 0.01, structure = None):
    """Get the elastic constants."""
    stresses, strains = [], []
    eq_stress = get_stress(atoms, potential_path = potential_path, potential_type = potential_type, device = device) * 1e-3 # Convert to GPa
    structure = AseAtomsAdaptor.get_structure(atoms)
    deform_set_sym = DeformedStructureSet(structure, symmetry=True)
    for i, structure in enumerate(deform_set_sym.deformed_structures):
        atoms_i = structure.to_ase_atoms()
        atoms_i = relax_atoms(atoms_i, potential_path = potential_path, potential_type = potential_type, device = device, fmax = fmax)
        stresses.append(get_stress(atoms_i, potential_path = potential_path, potential_type = potential_type, device = device) * 1e-3)
        strains.append(Strain.from_deformation(deform_set_sym.deformations[i]))
    
    cij = diff_fit(strains, stresses, eq_stress=eq_stress, order=2)[0].voigt
    if structure == 'hcp' or structure == 'omega':
        c11, c33, c12, c13, c44 = cij[0, 0], cij[2, 2], (cij[0, 1] + cij[1, 0]) / 2, (cij[0, 2] + cij[2, 0]) / 2, cij[3, 3]
        return {'c11': c11, 'c33': c33, 'c12': c12, 'c13': c13, 'c44': c44}
    elif structure == 'bcc' or structure == 'fcc':
        c11, c12, c44 = cij[0, 0], (cij[0, 1] + cij[1, 0]) / 2, cij[3, 3]
        return {'c11': c11, 'c12': c12, 'c44': c44}
    else:
        return cij
    
def compute_stacking_fault_energy(atoms, plane, layer, energy_reference, n_cells = [1, 1, 1], output=None, hard_plane = False, potential_path = None, potential_type = 'MACE', device = 'cpu', fmax = 0.005):
    """
    Compute the stacking fault energy for a given configuration.
    
    Returns:
        - Stacking fault energy in eV/Å².
    """
    if plane == "basal":
        index = [0,1,2]
        slip_vector = 0.5 / n_cells[0] * atoms.cell[:][0,:] - 1/6 / n_cells[1] * atoms.cell[:][1,:]
    elif plane == "prism":
        index =	[0,2,1]
        slip_vector = 0.5 / n_cells[0] * atoms.cell[:][0,:]
    else:
        raise ValueError(f"Invalid plane: {plane}")
    
    X_vector = atoms.cell[:][index[0],:]
    Y_vector = atoms.cell[:][index[1],:]
    Z_vector = atoms.cell[:][index[2],:]
    
    Z_vector_unit = Z_vector / np.linalg.norm(Z_vector)

    if hard_plane:
        dz = np.linalg.norm(atoms.cell[:][index[2],:]) / (n_cells[2] * 2) / 4
    else:
        dz = 0

    # Apply the shift to atoms above the slip plane
    z_cut = np.mean(np.dot(atoms.positions[:,:], Z_vector_unit)) / (n_cells[2] / 2) * layer + dz
    atoms.positions[np.dot(atoms.positions[:, :], Z_vector_unit) > z_cut] += slip_vector

    #Change periodic vector to have only one SF
    Z_vector_new = Z_vector + slip_vector
    new_cell = np.zeros((3,3))
    new_cell[index[0],:], new_cell[index[1],:], new_cell[index[2],:] = X_vector, Y_vector, Z_vector_new
    atoms.set_cell(new_cell , scale_atoms=False)

    # Apply a constraint to fix x and y positions (only allow z-movement)
    if plane != "basal":
        constraint = FixedLine(indices=[atom.index for atom in atoms], direction=Z_vector_new) #direction=atom_copy.cell[:][index[2],:])
        atoms.set_constraint(constraint)

    # Assign calculators (use a realistic potential)
    set_calculator(atoms, potential_path, potential_type, device)
    
    # Relax the faulted structure
    dyn = BFGS(atoms)
    dyn.run(fmax=fmax, steps=150)  # Converge forces to 0.005 eV/Å

    if output != None:
        atoms.write(output, format="vasp")
        
    # Compute the energies
    energy_faulted = atoms.get_potential_energy()

    # Compute stacking fault energy per area
    fault_energy_per_area = (energy_faulted - energy_reference) / atoms.cell.area(index[2])
    return fault_energy_per_area

def load_sqs(elements, composition, structure, i):
    """Load SQS structure."""

    a_lattice = {"hcp":{"Ti": 2.94, "Zr": 3.24, "Hf":3.19, "Sc":3.29, "Y":3.64}, "bcc":{"Ti" : 3.25, "Zr":3.58, "Hf": 3.53, "Sc":3.69, "Y":3.69}, "fcc" : {"Ti" : 4.11, "Zr":4.53, "Hf": 4.47, "Sc":4.65, "Y":5.10}, "omega":{"Ti": 2.83, "Zr": 3.09, "Hf":3.09, "Sc":3.20, "Y":3.20}}
    c_lattice = {"hcp":{"Ti":4.64, "Zr":5.17, "Hf":5.05, "Sc":5.33, "Y":5.88}, "omega":{"Ti":4.57, "Zr":4.96, "Hf":4.96, "Sc":5.06, "Y":5.06}} 

    a, c = 0, 0
    for i in range(len(elements)):
        a += a_lattice[structure][elements[i]] * composition[i]
        if structure == "hcp" or structure == "omega":
            c += c_lattice[structure][elements[i]] * composition[i]

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tmp_dir = os.path.join(project_root, 'tmp_sqs')
    os.makedirs(tmp_dir, exist_ok=True)
    # Sort elements and composition by composition values (ascending order)
    sorted_indices = np.argsort(composition)
    elements = [elements[i] for i in sorted_indices]
    composition = [composition[i] for i in sorted_indices]
    if len(elements) == 2:
        sqs_file = os.path.join(project_root, f"sqs_cells/binary/{structure}_A{composition[0]:.2f}_B{composition[1]:.2f}/bestsqs{i}.out")
        sqs_out = os.path.join(tmp_dir, f"{structure}_A{composition[0]:.2f}_B{composition[1]:.2f}_sqs_{i}.poscar")
    elif len(elements) == 3:
        sqs_file = os.path.join(project_root, f"sqs_cells/ternary/{structure}_A{composition[0]:.3f}_B{composition[1]:.3f}_C{composition[2]:.3f}/bestsqs{i}.out")
        sqs_out = os.path.join(tmp_dir, f"{structure}_A{composition[0]:.3f}_B{composition[1]:.3f}_C{composition[2]:.3f}_sqs_{i}.poscar")
    else:
        raise ValueError(f"Invalid number of elements: {len(elements)}")

    if structure == "hcp" or structure == "omega":
        run([sys.executable, os.path.join(project_root, 'utils', 'convert_sqs.py'),
            '--in_file', sqs_file,
            '--out_file', sqs_out,
            '--a', a,
            '--c', c,
            '--coordinate_type', 2])
    else:
        run([sys.executable, os.path.join(project_root, 'utils', 'convert_sqs.py'),
             '--in_file', sqs_file,
             '--out_file', sqs_out,
             '--a', a])
             
    atoms = read(sqs_out, format='vasp')
    atoms = alloy_sqs(atoms, elements, composition)
    os.remove(sqs_out)
    return atoms