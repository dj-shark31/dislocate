import os
import argparse

def read_first_lines(filename, n):
    with open(filename, "r") as f:
        return [f.readline().strip() for _ in range(n)]

def write_lammps_data(fout, tmp_lammps):
    with open(tmp_lammps, "r") as f:
        stress = f.readline().strip()
        energy = f.readline().split()[0]
    fout.write('# Cell energy and stress \n')
    fout.write(f"{energy} cell total energy (eV) \n")
    fout.write(f"{stress}\n")

def write_dxa_data(fout, tmp_dxa):
    with open(tmp_dxa, "r") as f:
        left_dis = f.readline().split()
        right_dis = f.readline().split()
    fout.write('# Dislocation position dxa \n')
    fout.write(f"{left_dis[0]} {left_dis[1]} x_left y_left \n")
    fout.write(f"{right_dis[0]} {right_dis[1]} x_right y_right \n")
    return [left_dis[:2], right_dis[:2]]

def write_fitting_data(fout, tmp_fitting, thickness, pbc):
    left_core, right_core = [], []
    with open(tmp_fitting, "r") as f:
        if pbc == 'true' and thickness > 1:
            ncore = 2 * thickness
        else:
            ncore = thickness
        for i in range(2 * ncore):
            data = f.readline().split()
            (left_core if i < ncore else right_core).append(data)
    fout.write('# Fitting data: core angle, eccentricity, area, x_pos, y_pos, r_squared fitting, core angle uncertainty, eccentricity uncertainty, area uncertainty \n')
    fout.write('# Left dislocation: \n')
    for row in left_core:
        fout.write(' '.join(row[:9]) + '\n')
    fout.write('# Right dislocation: \n')
    for row in right_core:
        fout.write(' '.join(row[:9]) + '\n')

def read_atom_data(tmp_stab, natom):
    atom_data = []
    with open(tmp_stab, "r") as f:
        for _ in range(natom):
            atom_data.append(f.readline().split())
    return atom_data

def add_nye_tensor(atom_data, tmp_babel, natom):
    with open(tmp_babel, "r") as f:
        lines = f.readlines()[1:]
    order = [4,5,6,7,8,9,10,11,12,3]
    for line in lines:
        data = line.split()
        for j in range(natom):
            if (int(float(data[0])*10) == int(float(atom_data[j][0])*10) and
                int(float(data[1])*10) == int(float(atom_data[j][1])*10) and
                int(float(data[2])*10) == int(float(atom_data[j][2])*10)):
                atom_data[j] += [data[i] for i in order]
                break

def write_atomic_data(fout, atom_data, data_out):
    fout.write('# Data from ovito and Babel (atomic positions in reference cell, atomic positions in dislocated cell, elastic stability parameter, pattern cna, Nye tensor, pattern babel): \n')
    fout.write('# ref_x ref_y ref_z dis_x dis_y dis_z elastic_stability pattern_cna nye_xx nye_xy nye_xz nye_yx nye_yy nye_yz nye_zx nye_zy nye_zz pattern_babel \n')
    for row in atom_data:
        fout.write(' '.join([row[k] for k in data_out]) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Assemble core characteristics and atomic data.')
    parser.add_argument('thickness', type=int, help='Cell thickness (number of unit cells in z)')
    parser.add_argument('a0', help='Lattice parameter a0 (z direction)')
    parser.add_argument('natom', type=int, help='Number of atoms in cell')
    parser.add_argument('tmp_babel', help='Temporary Babel file')
    parser.add_argument('tmp_stab', help='Temporary stability file')
    parser.add_argument('tmp_dxa', help='Temporary DXA file')
    parser.add_argument('tmp_fitting', help='Temporary fitting file')
    parser.add_argument('tmp_lammps', help='Temporary LAMMPS file')
    parser.add_argument('output_file', help='Output file')
    parser.add_argument('lammps', help='Whether to include LAMMPS data (true/false)')
    parser.add_argument('fitting', help='Whether to include fitting data (true/false)')
    parser.add_argument('ovito', help='Whether to include OVITO data (true/false)')
    parser.add_argument('nye', help='Whether to include Nye tensor data (true/false)')
    parser.add_argument('pbc', help='Whether periodic boundary conditions are used (true/false)')
    args = parser.parse_args()

    with open(args.output_file, "w") as fout:
        fout.write('# Cell characteristics \n')
        fout.write(f"{args.natom} atoms \n")
        fout.write(f"{args.thickness} thickness in z-direction in # of Burgers vectors \n")
        fout.write(f"{float(args.a0):.5f} Burgers vector \n")

        if args.lammps == 'true':
            write_lammps_data(fout, args.tmp_lammps)

        if args.ovito == 'true':
            write_dxa_data(fout, args.tmp_dxa)

        if args.fitting == 'true':
            write_fitting_data(fout, args.tmp_fitting, args.thickness, args.pbc)

        if args.ovito == 'true' or args.nye == 'true':
            atom_data = read_atom_data(args.tmp_stab, args.natom)
            data_out = [0, 1, 2, 3, 4, 5]
            if args.ovito == 'true':
                data_out += [6, 7]
            if args.nye == 'true':
                data_out += [8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
                add_nye_tensor(atom_data, args.tmp_babel, args.natom)
            write_atomic_data(fout, atom_data, data_out)

if __name__ == '__main__':
    main()

 
