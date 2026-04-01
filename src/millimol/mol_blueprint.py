import numpy as np
from molecule import Molecule


class MoleculeBlueprint:
    """
    Prepare molecule for drawing.
    Bisect bonds to draw each half in element color.
    Make list of unbound atoms to draw them separately.
    """
    def __init__(self, molecule: Molecule):
        # make a copy of numpy attributes
        self.coords = molecule.coords.copy()
        self.elements = molecule.elements.copy()
        self.n_atoms = self.elements.shape[0]
        assert self.n_atoms == self.coords.shape[0]
        
        # add virtual points in the middle of each bond
        midpoints = self.coords[molecule.bonds].mean(axis=1)
        self.coords = np.vstack((self.coords, midpoints))
        midpoint_indices = np.arange(self.n_atoms, self.n_atoms + molecule.bonds.shape[0])
        
        # connect each pair of bound atoms (u, v) to corresponding midpoint
        self.halfbonds = np.vstack((
            np.column_stack((molecule.bonds[:, 0], midpoint_indices)),
            np.column_stack((molecule.bonds[:, 1], midpoint_indices))
        ))
        
        # get indices of unbound atoms
        atom_indices = np.arange(self.n_atoms)
        self.unbound = atom_indices[
            np.isin(atom_indices, molecule.bonds.ravel(), invert=True)
        ]
    
    def linear_as_unbound(self):
        """Find atoms with linear geometry and mark them as unbound"""
        pass

