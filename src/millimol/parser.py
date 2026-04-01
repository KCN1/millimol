import re
from typing import TextIO, List
from pathlib import Path
from parser_config import FileFormat
from molecule import Molecule
from periodic_table import Element


class MoleculeParser:
    """
    Parse files with Cartesian coordinates of atoms.
    
    Accept Gaussian / ORCA input and log files, xyz files.
    Return list of Molecule objects with geometries found.
    """
    def __init__(self, filename: Path):
        N_LINES = 3
        self.filename = filename
        self.file_format = FileFormat.XYZ
        with open(filename, 'r') as fp:
            for _, line in zip(range(N_LINES), fp):
                for ff in FileFormat:
                    if ff.format_trigger and (ff.format_trigger in line):
                        self.file_format = ff

    def _seek_xyz(self, fp: TextIO) -> None:
        if self.file_format.read_trigger:
            for line in fp:
                if self.file_format.read_trigger in line:
                    break

    def _parse_xyz(self, fp: TextIO) -> Molecule:
        seek_flag = True
        element_list = []
        xyz_list = []
        for line in fp:
            if match := re.match(self.file_format.regex_, line):
                element_list.append(
                    Element._map_[match.group(self.file_format.group_[0])]
                )
                xyz_list.append(
                    tuple(map(match.group, self.file_format.group_[1:]))
                )
                seek_flag = False
            elif seek_flag:
                continue
            else:
                break
        return Molecule(element_list, xyz_list)

    def parse(self, calc_bonds=True, all_geoms=True) -> List[Molecule]:
        """
        Main parsing method.
        
        calc_bonds: whether to calculate bonds for each geometry,
        all_geoms: whether to return all geometries or the last one.
        """
        geometries = []
        with open(self.filename, 'r') as fp:
            self._seek_xyz(fp)
            while molecule := self._parse_xyz(fp):
                if all_geoms:
                    if calc_bonds:
                        molecule.calc_bonds()
                    geometries.append(molecule)
                else:
                    geometries = [molecule]
                self._seek_xyz(fp)
        if calc_bonds and not all_geoms:
            geometries[-1].calc_bonds()
        return geometries

