from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray
from scipy.spatial import KDTree
from periodic_table import PeriodicTable, Element


@dataclass
class Molecule:
    """
    Container class for atomic coordinates and bonds (pairs of atom indices).
    
    Accept element list as numpy array (N) or a list of integers / Element enums.
    Accept coordinates as numpy array (N, 3) or a list of tuples / lists.
    Accept bonds as numpy array of indices (atom numbering starts from 0).
    """
    elements: NDArray = np.zeros(0, dtype=np.int64)
    coords: NDArray = np.zeros((0, 3), dtype=np.float64)
    bonds: NDArray = np.zeros((0, 2), dtype=np.int64)
    
    def __post_init__(self):
        assert len(self.elements) == len(self.coords)
        if not isinstance(self.elements, np.ndarray):
            self.elements = np.fromiter(
                map(Element._map_.__getitem__, self.elements),
                count=len(self.elements),
                dtype=np.int64
            )
        self.coords = np.asarray(self.coords, dtype=np.float64)
    
    def __len__(self):
        return len(self.elements)
    
    @property
    def atoms(self) -> NDArray:
        """Return atoms as numpy records array (element number, x, y, z)"""
        return np.rec.fromarrays(
            (self.elements, *self.coords.T),
            names=('elements', 'x', 'y', 'z')
        )
    
    def calc_bonds(self) -> None:
        """Estimate bonds from covalent radii. Rewrites existing bonds!"""
        radii = PeriodicTable['radius'][self.elements]
        atom_tree = KDTree(self.coords)
        prelim_pairs = atom_tree.query_pairs(
            r = 2.2 * radii.max(),
            output_type='ndarray'
        )
        distances = np.linalg.norm(
            self.coords[prelim_pairs[:, 0]]
            - self.coords[prelim_pairs[:, 1]], 
            axis=1
        )
        self.bonds = prelim_pairs[
            distances < 1.1 * (
                radii[prelim_pairs[:, 0]]
                + radii[prelim_pairs[:, 1]]
            )
        ]
    
    def set_bonds(self, connectivity: dict) -> None:
        """
        Get bonds from connectivities (adjacency list).
        Atom numbering starts from 0!
        Rewrites existing bonds!
        """
        self.bonds = np.array([
            (min(u, v), max(u, v)) for u, v_list
            in connectivity.items() for v in v_list
        ], dtype=np.int64)

