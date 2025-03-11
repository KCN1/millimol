"""Creates a class for a molecule (atoms, coordinates, bonds) and draws the molecule using matplotlib.pyplot."""
from sys import argv
from math import dist
from typing import Tuple, List

from scipy.spatial import KDTree
import matplotlib.pyplot as plt


class Molecule:
    """
    Describes a molecule as a set of atoms and bonds (connectivity). Bonds are determined by interatomic distance.
    Virtual points in the middle of the bonds are added to color the bonds according to element colors.
    """
    _atom_numbers = (
        'Bq', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F',
        'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
        'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
        'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
        'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
        'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
        'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
        'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
        'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
        'Fm', 'Md', 'No', 'Lr'
    )
    _hyper_colors = {
            'C' : 'Cyan',
            'N' : 'Blue',
            'O' : 'Red',
            'F' : 'Yellow',
            'Na': 'Purple',
            'P' : 'Yellow',
            'S' : 'Yellow',
            'Cl': 'Yellow',
            'K' : 'Purple',
            'Fe': 'Red',
            'Co': 'Blue',
            'Cu': 'Green',
            'Br': 'Yellow',
            'I' : 'Red',
            'Au': 'Yellow',
            'Bq': 'Pink'
        }
    _radii = {
            'Bq': 0.00, 'H' : 0.35, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 
            'B' : 0.84, 'C' : 0.76, 'N' : 0.71, 'O' : 0.66, 'F' : 0.57, 
            'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P' : 1.07, 'S' : 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K' : 2.03, 
            'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V' : 1.53, 'Cr': 1.39, 
            'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 
            'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95, 'Y' : 1.90, 
            'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 
            'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 
            'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I' : 1.39, 'Xe': 1.40, 
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 
            'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 
            'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 
            'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W' : 1.62, 
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 
            'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 
            'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 
            'Th': 2.06, 'Pa': 2.00, 'U' : 1.96, 'Np': 1.90, 'Pu': 1.87, 
            'Am': 1.80, 'Cm': 1.69, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 
            'Fm': 1.67, 'Md': 1.73, 'No': 1.76, 'Lr': 1.61
        }

    def __init__(self, filename):
        self.element_colors = self._hyper_colors
        self.filename = filename
        self.connectivity = []
        self.virtual_connectivity = []
        self.virtual_points = []
        self.n, self.description, self.elements, self.points = 0, '', [], []
        self.unconnected = set()
        self.atom_colors = []
        self.read_xyz()

    def from_orca_gauss(self) -> Tuple[int, str, List[str], List[Tuple[float]]]:
        """Parses Orca, Gaussian and Priroda output formats."""
        n, description, p, el = 0, '', [], []
        flag = 'red'  # start/stop flag for reading the data
        log_format = ''
        geometries = 0  # to navigate through all geometries (for next versions)
        n0, p0, el0 = 0, [], []  # current number of atoms, coordinates, elements

        with open(self.filename, 'r') as file:

            for line in file:   # guess the log file format
                if 'Gaussian' in line:
                    log_format = 'gaussian'
                    break
                elif 'O   R   C   A' in line:
                    log_format = 'orca'
                    break
                elif 'Priroda' in line:
                    log_format = 'priroda'
                    break
                # TODO: add Gamess/Firefly format
                # elif 'GAMESS' in line:
                #     log_format = 'gamess'

            line = file.readline()
            if log_format == 'orca':    # check for a description
                while line and 'END OF INPUT' not in line:
                    if '|  1> #' in line:
                        description = line.lstrip('|  1> #').strip()
                    line = file.readline()
            elif log_format == 'gaussian':
                curr_line, prev_line, target_line = '', '', ''
                while line and 'Charge' not in line and 'orientation' not in line:
                    target_line = prev_line
                    prev_line = curr_line
                    curr_line = line
                    line = file.readline()
                if '-' in prev_line:
                    description = target_line.strip()
            elif log_format == 'priroda':
                while line and 'atoms' not in line:
                    if 'molecule input:' in line:
                        description = line.split("'")[1]
                    line = file.readline()
            else:
                description = 'Unknown format'

            for line in file:
                if flag == 'red':
                    if log_format == 'orca' and 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                        flag = 'yellow'
                    elif log_format == 'gaussian' and 'Coordinates (Angstroms)' in line:
                        flag = 'yellow'
                    elif log_format == 'priroda' and 'Atomic Coordinates:' in line:
                        flag = 'green'
                elif flag == 'yellow' and '-------' in line:
                    flag = 'green'
                elif flag == 'green':
                    if line.strip() and '-------' not in line and '#' not in line:
                        s = line.split()
                        p0.append((float(s[-3]), float(s[-2]), float(s[-1])))
                        if log_format == 'gaussian':
                            el0.append(self._atom_numbers[int(s[1])])
                        elif log_format == 'priroda':
                            el0.append(self._atom_numbers[int(s[0])])
                        elif log_format == 'orca':
                            el0.append(s[0])
                        n0 += 1
                    else:
                        flag = 'red'
                        n = n0  # last number of atoms
                        p = p0.copy()  # last coordinates
                        el = el0.copy()  # last list of els
                        n0, p0, el0 = 0, [], []
                        geometries += 1
        return n, description, el, p

    def from_trj_xyz(self):
        """Parses .xyz atom coordinate format."""
        n, description, p, el = 0, '', [], []
        geometries = 0
        n0, p0, el0 = 0, [], []  # current number of atoms, their coordinates, elements

        with open(self.filename, 'r') as file:

            line = file.readline()
            while line:
                n0 = int(line.strip())
                description = file.readline().strip()

                # create an element array and a coordinate array
                for i in range(n0):
                    s = file.readline().split()
                    p_i = (float(s[-3]), float(s[-2]), float(s[-1]))
                    p0.append(p_i)
                    el0.append(s[0])
                n = n0
                p = p0.copy()
                el = el0.copy()
                n0, p0, el0 = 0, [], []
                geometries += 1
                line = file.readline()
        return n, description, el, p

    def read_xyz(self):
        """Creates a connectivity list from a list of coordinates"""

        # check atomic coordinates file format
        if self.filename[-4:] == '.xyz':
            self.n, self.description, self.elements, self.points = self.from_trj_xyz()
        elif self.filename[-4:] in ('.log', '.out'):
            self.n, self.description, self.elements, self.points = self.from_orca_gauss()
        else:
            self.n, self.description, self.elements, self.points = 0, '', [], []

        # define element colors
        self.atom_colors = ['LightGray'] * self.n
        for i in range(self.n):
            if self.elements[i] in self.element_colors:
                self.atom_colors[i] = self.element_colors[self.elements[i]]

        # build a KD-Tree to search for pairs of atoms within (2 * max(radii) + 10%)
        if self.n:
            atom_tree = KDTree(self.points, leafsize=100)
            prelim_pairs = atom_tree.query_pairs(r = 2.2 * max((self._radii[i] for i in self.elements)))
        else:
            prelim_pairs = []

        # create connectivity array
        self.connectivity = []
        self.virtual_connectivity = []
        self.virtual_points = []
        self.unconnected = set(range(self.n))

        for (i, j) in prelim_pairs:
            if dist(self.points[i], self.points[j]) < 1.1 * (self._radii[self.elements[i]] + self._radii[self.elements[j]]):
                
                self.connectivity.append((i, j))

                # add virtual points in the middle of each bond
                pv = ((self.points[i][0] + self.points[j][0]) / 2, (self.points[i][1] + self.points[j][1]) / 2, (self.points[i][2] + self.points[j][2]) / 2)
                self.virtual_points.append(pv)
                
                # connect last pair of bonded atoms to the last virtual point
                self.virtual_connectivity.append((i, len(self.virtual_points) - 1))
                self.virtual_connectivity.append((j, len(self.virtual_points) - 1))

                # remove both from the list of unbound atoms
                self.unconnected.discard(i)
                self.unconnected.discard(j)


def main(molecule: Molecule):
    """Draws a molecule using pyplot"""

    # pyplot window settings
    plt.style.use('dark_background')
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    plt.subplots_adjust(0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
    fig.set_size_inches(6.4, 6.4)

    # plot bonds between atoms and middle points
    for (pi, pvi) in molecule.virtual_connectivity:
        x1, y1, z1 = molecule.points[pi]
        x2, y2, z2 = molecule.virtual_points[pvi]
        ax.plot((x1, x2), (y1, y2), (z1, z2), color=molecule.atom_colors[pi])

    # plot unbound atoms as points
    for pi in molecule.unconnected:
        x0, y0, z0 = molecule.points[pi]
        ax.scatter(x0, y0, z0, color=molecule.atom_colors[pi])

    # subplot settings
    ax.set_aspect('equal', 'box')
    ax._axis3don = False

    plt.show()


if __name__ == '__main__':
    # TODO: add async timer to load the file periodically
    if len(argv) >= 2:
        filename = argv[-1]
    else:
        filename = input('Enter filename\n')
    molecule = Molecule(filename)
    main(molecule)
