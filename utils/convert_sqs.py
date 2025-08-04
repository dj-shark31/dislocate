'''
This program reads bestsqs file from user and generates structure file
either in quantum espresso or VASP format
To run use "python3 convert_sqs.py bestsqs.out"
'''

import sys
import numpy as np
import os

def read_bestsqs(file_name):
    with open(file_name, 'r') as reader:
        data = reader.read().rstrip().split("\n")
    return data

def process_bestsqs(data):
    def extract_vals(s_i, e_i):
        vec = []
        for i in range(s_i, e_i):
            vec.append(data[i].split())
        return vec

    unit_cell = np.array(extract_vals(0, 3)).astype(float)
    lattice_vectors = np.array(extract_vals(3, 6)).astype(float)
    atomic_positions = np.array(extract_vals(6, len(data)))

    elem_arr = np.unique(atomic_positions[:, 3])
    a_pos = np.array([['', '', '', '']])
    num_atoms_arr = []
    for elem in elem_arr:
        tmp_arr = atomic_positions[atomic_positions[:, 3] == elem]
        num_atoms_arr.append(len(tmp_arr))
        a_pos = np.concatenate((a_pos, tmp_arr), axis=0)
    a_pos = a_pos[1:]

    l_vec = np.matmul(lattice_vectors, unit_cell)
    a_pos = np.matmul(a_pos[:, 0:3].astype(float), unit_cell)

    return l_vec, a_pos, elem_arr, num_atoms_arr

def create_qe_inp(l_vec, a_pos, elem_arr, num_atoms_arr, coordinate_type, unit_cell_inv):
    str_arr = [
        '&CONTROL',
        '/',
        '&SYSTEM',
        '  ibrav = 0',
        f'  nat = {len(a_pos)}',
        f'  ntyp = {len(elem_arr)}',
        '/',
        '&ELECTRONS',
        '/',
        'ATOMIC_SPECIES'
    ]
    for elem in elem_arr:
        str_arr.append(f"  {elem} xx xxxx")

    str_arr.append("CELL_PARAMETERS angstrom")
    for i in range(0, 3):
        str_arr.append(f'  {l_vec[i,0]:.10f}  {l_vec[i,1]:.10f}  {l_vec[i,2]:.10f}')

    str_arr.append(f'ATOMIC_POSITIONS {coordinate_type}')
    ind = 0
    arr_ind = 0
    for a in a_pos:
        if coordinate_type == 'crystal':
            a = np.matmul(a, unit_cell_inv)
        if(ind == num_atoms_arr[arr_ind]):
            arr_ind += 1
            ind = 0
        str_arr.append(f'{elem_arr[arr_ind]}  {a[0]:.10f}  {a[1]:.10f}  {a[2]:.10f}')
        ind += 1

    str_arr.append("K_POINTS automatic")
    str_arr.append("  1 1 1 0 0 0")

    return str_arr

def create_poscar(l_vec, a_pos, elem_arr, num_atoms_arr, coordinate_type, unit_cell_inv):
    str_arr = [
        'POSCAR',
        '1.0'
    ]

    for i in range(0, 3):
        str_arr.append(f'  {l_vec[i,0]:.10f}  {l_vec[i,1]:.10f}  {l_vec[i,2]:.10f}')

    str_arr.append('  ' + ' '.join(elem_arr))
    str_arr.append('  ' + ' '.join(map(str, num_atoms_arr)))
    str_arr.append(coordinate_type.capitalize())

    for a in a_pos:
        if coordinate_type == 'fractional' or coordinate_type == 'Direct':
            a = np.matmul(a, unit_cell_inv)
        str_arr.append(f'  {a[0]:.10f}  {a[1]:.10f}  {a[2]:.10f}')

    return str_arr

def save_to_file(out_file, str_arr):
    with open(out_file, 'w') as f:
        f.write('\n'.join(str_arr))

def main():
    parser = argparse.ArgumentParser(description='Convert SQS output to POSCAR/QE format')
    parser.add_argument('--in_file', type=str, required=True, help='Input bestsqs.out file')
    parser.add_argument('--out_file', type=str, required=True, help='Output file')
    parser.add_argument('--a', type=float, required=True, help='a lattice parameter')
    parser.add_argument('--c', type=float, required=False, help='c lattice parameter')
    parser.add_argument('--coordinate_type', type=int, help='Coordinate type (1: cartesian, 2: fractional)')
    args = parser.parse_args()

    data = read_bestsqs(args.in_file)
    l_vec, a_pos, elem_arr, num_atoms_arr = process_bestsqs(data)

    str_out = 2  # Creates POSCAR for VASP
    coordinate_type = 2  # Fractional coordinates
    a_lat = args.a
    if args.c:
        c_lat = args.c
    else:
        c_lat = a_lat

    #Re-scale lattice vectors and atomic positions
    if c_lat != a_lat:
        l_vec = np.matmul(l_vec, np.diagflat([a_lat, a_lat, c_lat/1.632990]))
        a_pos = np.matmul(a_pos, np.diagflat([a_lat, a_lat, c_lat/1.632990]))
    else:
        l_vec = np.matmul(l_vec, np.diagflat([a_lat, a_lat, a_lat]))
        a_pos = np.matmul(a_pos, np.diagflat([a_lat, a_lat, a_lat]))

    if coordinate_type == 1:
        coordinate_type = 'cartesian'
    elif coordinate_type == 2:
        coordinate_type = 'fractional'
    else:
        print("Invalid coordinate type option. Exiting.")
        sys.exit(1)

    unit_cell_inv = np.linalg.inv(l_vec)

    if str_out == 1:
        coordinate_type = 'crystal' if coordinate_type == 'fractional' else coordinate_type
        str_arr = create_qe_inp(l_vec, a_pos, elem_arr, num_atoms_arr, coordinate_type, unit_cell_inv)
    elif str_out == 2:
        coordinate_type = 'Direct' if coordinate_type == 'fractional' else coordinate_type
        str_arr = create_poscar(l_vec, a_pos, elem_arr, num_atoms_arr, coordinate_type, unit_cell_inv)

    save_to_file(args.out_file, str_arr)

if __name__ == "__main__":
    main()
