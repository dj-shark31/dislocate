#!/usr/bin/env python3
"""
Python translation of fitting_core.wls and fitting_core_pbc.wls
Usage:
    python fitting_core.py outStab outDXA thickness outFitting [--pbc]
"""
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# 2D elliptical Gaussian with rotation
# (x, y, v) data, fit params: theta, s1, s2, dis1, dis2

def elliptical_gaussian(XY, max_v, theta, s1, s2, dis1, dis2):
    x, y = XY
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    x_shifted = x - dis1
    y_shifted = y - dis2
    a = (cos_theta**2)/(2*s1**2) + (sin_theta**2)/(2*s2**2)
    b = (1/s1**2 - 1/s2**2)*sin_theta*cos_theta/2
    c = (sin_theta**2)/(2*s1**2) + (cos_theta**2)/(2*s2**2)
    return max_v * np.exp(-(a*x_shifted**2 + 2*b*x_shifted*y_shifted + c*y_shifted**2))

def fit_gaussian(core_atoms, dislocation_pos, vmin=0, vmax=2):
    # core_atoms: Nx3 array (x, y, v)
    core_data = np.copy(core_atoms)
    core_data[:,2] = np.where(core_data[:,2] < vmin, 0, core_data[:,2] - vmin)
    core_data[:,2] = np.where(core_data[:,2] < vmax, core_data[:,2], vmax)
    max_v = np.max(core_data[:,2])
    x = core_data[:,0]
    y = core_data[:,1]
    v = core_data[:,2]
    # Only fit theta, s1, s2, dis1, dis2; max_v is fixed
    def elliptical_gaussian_fixed_max(XY, theta, s1, s2, dis1, dis2):
        return elliptical_gaussian(XY, max_v, theta, s1, s2, dis1, dis2)
    initial_guess = [1.57, 1.5, 3.0, dislocation_pos[0], dislocation_pos[1]]
    bounds = ([0, 0.1, 0.1, dislocation_pos[0]-10, dislocation_pos[1]-10],
              [np.pi, 20, 20, dislocation_pos[0]+10, dislocation_pos[1]+10])
    try:
        popt, pcov = curve_fit(
            elliptical_gaussian_fixed_max,
            (x, y), v, p0=initial_guess, bounds=bounds, maxfev=1000000
        )
        # Calculate R^2
        residuals = v - elliptical_gaussian_fixed_max((x, y), *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((v - np.mean(v))**2)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
        # Calculate uncertainties
        perr = np.sqrt(np.diag(pcov)) if pcov is not None else [0]*5
    except Exception:
        popt = [np.nan]*5
        r2 = 0
        perr = [np.nan]*5
    theta, s1, s2, dis1, dis2 = popt
    theta_deg = np.mod(theta, np.pi)*180/np.pi
    aspect_ratio = abs(s2/s1) if s1 != 0 else np.nan
    area = abs(np.pi*s1*s2)

    core_properties = [theta_deg, aspect_ratio, area, dis1, dis2, r2, perr[0]*180/np.pi, perr[1], perr[2]]
    if aspect_ratio > 1:
        return core_properties
    else:
        return [np.mod(theta_deg+90,180), 1/aspect_ratio, area] + core_properties[3:]

def main():
    parser = argparse.ArgumentParser(description='Fit core properties from stability and DXA data.')
    parser.add_argument('outStab', help='Input stability data file (tabular)')
    parser.add_argument('outDXA', help='Input DXA data file (tabular)')
    parser.add_argument('thickness', type=int, help='Cell thickness')
    parser.add_argument('outFitting', help='Output file for fit parameters')
    parser.add_argument('--pbc', type=str, default='false', help='Enable periodic boundary condition logic')
    args = parser.parse_args()

    # Parameters
    core_radius = 20
    vmin = 0
    vmax = 2

    # Read and sort stability data
    stability_df = pd.read_csv(args.outStab, sep='\s+', header=None, comment='#')
    # Only keep columns 0,1,2,6 (x, y, z, v)
    stability_data = stability_df[[0,1,2,6]].values
    stability_data = stability_data[np.argsort(stability_data[:,2])]
    zmin = stability_data[0,2]
    num_atoms = stability_data.shape[0]

    if args.pbc == 'true':
        atoms_per_slice = int(round(num_atoms/args.thickness/2))
    else: # Standard logic
        atoms_per_slice = int(round(num_atoms/args.thickness))

    # Check atom ordering
    if stability_data[atoms_per_slice//2,2] - zmin > 1:
        iup = 1
        while stability_data[iup,2] - zmin < 1 and iup < stability_data.shape[0]-1:
            iup += 1
        stability_data = np.vstack([stability_data[iup:], stability_data[:iup]])

    # Read dislocation positions
    dxa_df = pd.read_csv(args.outDXA, sep='\s+', header=None, comment='#')
    left_dislocation_pos = dxa_df.iloc[0,:2].values
    right_dislocation_pos = dxa_df.iloc[1,:2].values

    # Get core atom groups
    if args.pbc == 'true':
        # PBC logic
        left_cores_perb = [stability_data[i*atoms_per_slice:(i+1)*atoms_per_slice][:,[0,1,3]] for i in range(2*args.thickness)]
        right_cores_perb = [stability_data[i*atoms_per_slice:(i+1)*atoms_per_slice][:,[0,1,3]] for i in range(2*args.thickness)]
        left_cores_perb_filtered = []
        for core in left_cores_perb:
            mask = [np.linalg.norm(core[i,:2] - left_dislocation_pos) < core_radius for i in range(core.shape[0])]
            left_cores_perb_filtered.append(core[mask])
        left_cores_perb = left_cores_perb_filtered
        right_cores_perb_filtered = []
        for core in right_cores_perb:
            mask = [np.linalg.norm(core[i,:2] - right_dislocation_pos) < core_radius for i in range(core.shape[0])]
            right_cores_perb_filtered.append(core[mask])
        right_cores_perb = right_cores_perb_filtered
        if args.thickness == 1:
            ordering = [0]
        else:
            ordering = [i for i in range(0,2*args.thickness,2)]
        left_cores = [np.vstack([left_cores_perb[i], left_cores_perb[(i+1)%(2*args.thickness)]]) for i in ordering]
        right_cores = [np.vstack([right_cores_perb[i], right_cores_perb[(i+1)%(2*args.thickness)]]) for i in ordering]
    else:
        # Standard logic
        left_cores = [stability_data[i*atoms_per_slice:(i+1)*atoms_per_slice][:,[0,1,3]] for i in range(args.thickness)]
        right_cores = [stability_data[i*atoms_per_slice:(i+1)*atoms_per_slice][:,[0,1,3]] for i in range(args.thickness)]
        left_cores_filtered = []
        for core in left_cores:
            mask = [np.linalg.norm(core[i,:2] - left_dislocation_pos) < core_radius for i in range(core.shape[0])]
            left_cores_filtered.append(core[mask])
        left_cores = left_cores_filtered
        right_cores_filtered = []
        for core in right_cores: 
            mask = [np.linalg.norm(core[i,:2] - right_dislocation_pos) < core_radius for i in range(core.shape[0])]
            right_cores_filtered.append(core[mask])
        right_cores = right_cores_filtered

    # Fit Gaussian to each core
    left_core_fits = [fit_gaussian(core, left_dislocation_pos, vmin, vmax) for core in left_cores]
    right_core_fits = [fit_gaussian(core, right_dislocation_pos, vmin, vmax) for core in right_cores]
    all_core_fits = np.vstack([left_core_fits, right_core_fits])
    np.savetxt(args.outFitting, all_core_fits, fmt='%.4f')

if __name__ == "__main__":
    main() 