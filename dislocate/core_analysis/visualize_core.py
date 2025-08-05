import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.patches import FancyArrowPatch

class Core:
    def __init__(self, file):
        self.file = file
        self.left_dislocation_position_dxa = None 
        self.right_dislocation_position_dxa = None  
        self.left_dislocation_position_babel = None 
        self.right_dislocation_position_babel = None 
        self.energy = None
        self.stress = None
        self.n_atoms = None
        self.thickness = None
        self.burgers_vector = None
        self.core_angle = None
        self.core_ecc = None
        self.core_area = None
        self.core_angle_uncertainty = None
        self.core_ecc_uncertainty = None
        self.core_area_uncertainty = None
        self.r_squared = None
        self.data_map = None
        self._parse_data()

    def _parse_data(self):
        with open(self.file, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if "Cell characteristics" in line:
                self.n_atoms = int(lines[i+1].split()[0])
                self.thickness = int(lines[i+2].split()[0])
                self.burgers_vector = float(lines[i+3].split()[0]) # Burgers vector in A
            elif "Cell energy and stress" in line:
                self.energy = float(lines[i+1].split()[0]) # Energy in eV
                self.stress = [float(x) for x in lines[i+2].split()[:6]] # Stress in voigt notation (xx, yy, zz, yz, xz, xy) in MPa
            elif "Dislocation position dxa" in line:
                self.left_dislocation_position_dxa = np.array([float(x) for x in lines[i+1].split()[:2]])
                self.right_dislocation_position_dxa = np.array([float(x) for x in lines[i+2].split()[:2]])
            elif "Fitting data" in line:
                self.core_angle = [float(lines[i+2].split()[0]), float(lines[i+4].split()[0])]
                self.core_ecc = [float(lines[i+2].split()[1]), float(lines[i+4].split()[1])]
                self.core_area = [float(lines[i+2].split()[2]), float(lines[i+4].split()[2])]
                self.core_angle_uncertainty = [float(lines[i+2].split()[-3]), float(lines[i+4].split()[-3])]
                self.core_ecc_uncertainty = [float(lines[i+2].split()[-2]), float(lines[i+4].split()[-2])]
                self.core_area_uncertainty = [float(lines[i+2].split()[-1]), float(lines[i+4].split()[-1])]
                self.r_squared = [float(lines[i+2].split()[-4]), float(lines[i+4].split()[-4])]
                self.left_dislocation_position_babel = np.array([float(lines[i+2].split()[3]), float(lines[i+2].split()[4])])
                self.right_dislocation_position_babel = np.array([float(lines[i+4].split()[3]), float(lines[i+4].split()[4])])
            elif "Data from ovito and Babel" in line:
                self.data_map = pd.DataFrame([line.split() for line in lines[i+2:]], columns=lines[i+1].split()[1:]).astype(float)

    def find_closest_points(self, index):
        """
        Find the 6 closest points to a given point in the xy plane.
        
        Args:
            index: index of the point to find the closest points to
            
        Returns:
            indices: array of indices for the 6 closest points
        """
        # Calculate distances from point to all ref_points using only x,y coordinates
        distances = np.sqrt(np.sum((self.data_map[['ref_x', 'ref_y']].values - self.data_map[['ref_x', 'ref_y']].values[index])**2, axis=1))
        
        # Get indices of 6 closest points
        closest_indices = np.argsort(distances)[1:7]
        for i in closest_indices:
            if distances[i] > 1.5 * distances[closest_indices[0]]:
                closest_indices = closest_indices[:i]
                break
        
        return closest_indices

    def differential_displacement(self, i, j):
        dd = (self.data_map['ref_z'][i] - self.data_map['ref_z'][j]) - (self.data_map['dis_z'][i] - self.data_map['dis_z'][j])
        while abs(dd + self.burgers_vector) < abs(dd):
            dd += self.burgers_vector
        while abs(dd - self.burgers_vector) < abs(dd):
            dd -= self.burgers_vector
        return dd

    def arrow_coordinates(self, i, j, dd_ij):
        mid_point = (self.data_map[['ref_x', 'ref_y']].values[i] + self.data_map[['ref_x', 'ref_y']].values[j])/2
        direction = (self.data_map[['ref_x', 'ref_y']].values[i] - self.data_map[['ref_x', 'ref_y']].values[j])/np.linalg.norm(self.data_map[['ref_x', 'ref_y']].values[i] - self.data_map[['ref_x', 'ref_y']].values[j])
        end_point = mid_point + dd_ij * 2 /2 * direction
        start_point = mid_point - dd_ij * 2 /2 * direction
        return [start_point, end_point]

    def vitek_map(self):
        """Translate of VitekMap function"""
        
        neighbors = [self.find_closest_points(i) for i in range(len(self.data_map))]
        dd = [[self.differential_displacement(i, j) for j in neighbors[i]] for i in range(len(self.data_map))]
        
        arrows = []
        # Plot arrows
        for i in range(len(self.data_map)):
            for k, j in enumerate(neighbors[i]):
                arrow_size = abs(dd[i][k]) if abs(dd[i][k]) > self.burgers_vector / 10 else 0
                     
                if arrow_size > 0:
                    start_end = self.arrow_coordinates(i, j, dd[i][k])
                    arrows.append(FancyArrowPatch(
                        tuple(start_end[0]), tuple(start_end[1]),
                        arrowstyle='->', 
                        mutation_scale=arrow_size*15,
                        linewidth=1.5,
                        color='k'))
        return arrows

    def _density_plot(self, value_name, title=None, cmap='plasma', dislocation='left', x_range=10, y_range=5):
        # Add dislocation argument to specify which core to plot around
        if dislocation not in ['left', 'right']:
            raise ValueError("dislocation must be 'left' or 'right'")
            
        # Get reference position based on specified dislocation
        if dislocation == 'left':
            if self.left_dislocation_position_babel is not None:
                dis_pos = self.left_dislocation_position_babel
            elif self.left_dislocation_position_dxa is not None:
                dis_pos = self.left_dislocation_position_dxa
            else:
                raise ValueError("Left dislocation position not found")
        elif dislocation == 'right':
            if self.right_dislocation_position_babel is not None:
                dis_pos = self.right_dislocation_position_babel
            elif self.right_dislocation_position_dxa is not None:
                dis_pos = self.right_dislocation_position_dxa
            else:
                raise ValueError("Right dislocation position not found")
            
        # Filter atoms within +-10 in x and +-5 in y of specified dislocation
        mask = [(abs(self.data_map['ref_x'].values[i] - dis_pos[0]) <= x_range) & (abs(self.data_map['ref_y'].values[i] - dis_pos[1]) <= y_range) for i in range(len(self.data_map))]
        
        # Filter positions and values
        x = self.data_map[mask]['ref_x']
        y = self.data_map[mask]['ref_y']
        z = self.data_map[mask]['ref_z']
        v = self.data_map[mask][value_name]
        
        # Split atoms based on z coordinate
        z_mean = z.mean()
        below_mask = z < z_mean
        above_mask = z >= z_mean
        
        # Create grid
        xi = np.linspace(x.min(), x.max(), 200)
        yi = np.linspace(y.min(), y.max(), 200)
        xi, yi = np.meshgrid(xi, yi)
        
        # Interpolate
        vi = griddata((x, y), v, (xi, yi), method='cubic')
        plt.figure(figsize=(7, 6))
            
        plt.pcolormesh(xi, yi, vi, shading='auto', cmap=cmap)
        im = plt.gca().get_children()[0]
        plt.colorbar(im, label=title, shrink=0.5)  # Make colorbar smaller

        if value_name == "nye_zz":
            arrows = self.vitek_map()
            for arrow in arrows:
                plt.gca().add_patch(arrow)
                
        # Add contour lines if plotting elastic stability parameter
        plt.clabel(plt.contour(xi, yi, vi, colors='k', alpha=0.5, linewidths=1.0), inline=True, fmt='%1.2f')
            
        # Maintain aspect ratio and increase dot size
        plt.gca().set_aspect('equal')
        plt.scatter(x[below_mask], y[below_mask], c='r', s=40, alpha=0.8)
        plt.scatter(x[above_mask], y[above_mask], c='none', edgecolors='r', s=40, alpha=0.8)  
                
        plt.xlim(xi.min(), xi.max())
        plt.ylim(yi.min(), yi.max())
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(title)
        plt.gcf().set_size_inches(10, 8)  # Make plot larger
        plt.tight_layout()
        plt.show()

    def plot_dislocation_density(self, title='Dislocation density', dislocation='left', x_range=10, y_range=5):
        self._density_plot("nye_zz", title=title, dislocation=dislocation, x_range=x_range, y_range=y_range)

    def plot_edge_density_xz(self, title='Edge density xz component', dislocation='left', x_range=10, y_range=5):
        self._density_plot("nye_xz", title=title, dislocation=dislocation, x_range=x_range, y_range=y_range)

    def plot_edge_density_yz(self, title='Edge density yz component', dislocation='left', x_range=10, y_range=5):
        self._density_plot("nye_yz", title=title, dislocation=dislocation, x_range=x_range, y_range=y_range)

    def plot_elastic_stability(self, title='Elastic stability parameter', dislocation='left', x_range=10, y_range=5):
        self._density_plot("elastic_stability", title=title, dislocation=dislocation, x_range=x_range, y_range=y_range) 