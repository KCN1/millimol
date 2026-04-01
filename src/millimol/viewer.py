#!/usr/bin/env python3

from pathlib import Path
from blueprint import MoleculeBlueprint
from plot_config import COLORMAP
from molecule import Molecule
from parser import MoleculeParser
from plot_pyvista import draw_molecule


def view_molecule(filename: Path) -> None:
    """
    Driver function of molecule viewer.
    
    1. Parse file to get last geometry.
    2. Prepare molecule for drawing.
    3. Draw molecule with Pyvista.
    """
    mol_parser = MoleculeParser(filename)
    molecule = mol_parser.parse(all_geoms=False)[-1]
    bp = MoleculeBlueprint(molecule)
    draw_molecule(
        Path(filename).name,
        bp.coords,
        bp.halfbonds,
        bp.unbound,
        bp.elements,
        COLORMAP,
        'Black'
    )


if __name__ == '__main__':
    from argparse import ArgumentParser
    
    arg_parser = ArgumentParser()
    arg_parser.add_argument("file_path", type=Path)
    filename = arg_parser.parse_args().file_path
    
    view_molecule(filename)

