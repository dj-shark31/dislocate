from ase.calculators.lammpslib import LAMMPSlib
from mace.calculators import MACECalculator
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from pymatgen.analysis.elasticity import diff_fit, Strain, DeformedStructureSet

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
        elif symbol == "Hf":
            symbols[k] = elements[0]
        elif symbol == "Ti":
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
        lat = Lattice.from_parameters(2.86, 2.86, 4.77, 90, 90, 120)
        cell = Structure(lat, ["Ti", "Ti"], [[0.0, 0.0, 0.0], [0.333333, 0.666667, 0.5]])
    elif structure == 'omega':
        lat = Lattice.from_parameters(2.86, 2.86, 4.77, 90, 90, 120)
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

def get_lattice_parameters(atoms, n_cells = [1, 1, 1], potential_path = None, potential_type = 'MACE', device = 'cpu', dof = [True, True, True, False, False, False], fmax = 0.004, phase = 'hcp'):
    """Get the lattice parameters."""
    atoms = relax_cell(atoms, potential_path = potential_path, potential_type = potential_type, device = device, dof = dof, fmax = fmax)
    a1, a2, a3 = get_unit_cell_vectors(atoms, n_cells)
    if phase == 'hcp' or phase == 'omega':
        a = (a1 + a2) / 2
        c = a3
        return {'a': a, 'c': c}
    else:
        a = (a1 + a2 + a3) / 3
        return {'a': a}

def get_cij(atoms, potential_path = None, potential_type = 'MACE', device = 'cpu', fmax = 0.01, phase = None):
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
    if phase == 'hcp' or phase == 'omega':
        c11, c33, c12, c13, c44 = cij[0, 0], cij[2, 2], (cij[0, 1] + cij[1, 0]) / 2, (cij[0, 2] + cij[2, 0]) / 2, cij[3, 3]
        return {'c11': c11, 'c33': c33, 'c12': c12, 'c13': c13, 'c44': c44}
    elif phase == 'bcc' or phase == 'fcc':
        c11, c12, c44 = cij[0, 0], (cij[0, 1] + cij[1, 0]) / 2, cij[3, 3]
        return {'c11': c11, 'c12': c12, 'c44': c44}
    else:
        return cij
    
