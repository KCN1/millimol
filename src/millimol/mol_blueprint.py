import numpy as np
from numpy.typing import NDArray
from molecule import Molecule


class MoleculeBlueprint:
    """
    Prepare the molecule for drawing:
    Bisect bonds to draw each half in each atom's color.
    Make an array of unbound atoms to draw them separately.
    
    coords: coordinates of "real" points (atoms) and "virtual" points in the middle of each bond
    elements: array of element numbers (for "real" atoms only)
    n_atoms: number of atoms (excluding "virtual" middle points)
    halfbonds: pairs of indices connecting "real" and "virtual" middle points
    unbound: array of atoms without bonds
    """
    def __init__(self, molecule: Molecule):
        # make a copy of numpy attributes
        self.coords: NDArray = molecule.coords.copy()
        self.elements: NDArray = molecule.elements.copy()
        self.n_atoms: int = self.elements.shape[0]
        assert self.n_atoms == self.coords.shape[0]
        
        # add virtual points in the middle of each bond
        midpoints = self.coords[molecule.bonds].mean(axis=1)
        self.coords = np.vstack((self.coords, midpoints))
        midpoint_indices = np.arange(self.n_atoms, self.n_atoms + molecule.bonds.shape[0])
        
        # connect each pair of bound atoms (u, v) to corresponding midpoint
        self.halfbonds: NDArray = np.vstack((
            np.column_stack((molecule.bonds[:, 0], midpoint_indices)),
            np.column_stack((molecule.bonds[:, 1], midpoint_indices))
        ))
        
        # get indices of unbound atoms
        atom_indices = np.arange(self.n_atoms)
        self.unbound: NDArray = atom_indices[
            np.isin(atom_indices, molecule.bonds.ravel(), invert=True)
        ]
    
    def linear_as_unbound(self):
        """Find atoms with linear geometry and mark them as unbound"""
        pass

