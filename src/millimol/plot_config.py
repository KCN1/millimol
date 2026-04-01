from periodic_table import PeriodicTable, Element


HYPERCHEM_COLORS = {
    'Bq': 'Pink',
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
    'Au': 'Yellow'
}


COLORS = HYPERCHEM_COLORS
DEFAULT_COLOR = 'White'


COLORMAP = [DEFAULT_COLOR] * len(PeriodicTable)
for name, color in COLORS.items():
    COLORMAP[Element._map_[name]] = color

