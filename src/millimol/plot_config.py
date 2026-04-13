from periodic_table import PeriodicTable, Element


HYPERCHEM_COLORS = {
    'Bq': [255, 192, 203], # Pink
    'C' : [0, 255, 255],   # Cyan
    'N' : [0, 0, 255],     # Blue
    'O' : [255, 0, 0],     # Red
    'F' : [255, 255, 0],   # Yellow
    'Na': [128, 0, 128],   # Purple
    'P' : [255, 255, 0],   # Yellow
    'S' : [255, 255, 0],   # Yellow
    'Cl': [255, 255, 0],   # Yellow
    'K' : [128, 0, 128],   # Purple
    'Fe': [255, 0, 0],     # Red
    'Co': [0, 0, 255],     # Blue
    'Cu': [0, 128, 0],     # Green
    'Br': [255, 255, 0],   # Yellow
    'I' : [255, 0, 0],     # Red
    'Au': [255, 255, 0]    # Yellow
}


COLORS = HYPERCHEM_COLORS
DEFAULT_COLOR = [255, 255, 255] # White


COLORMAP = [DEFAULT_COLOR] * len(PeriodicTable)
for name, color in COLORS.items():
    COLORMAP[Element._map_[name]] = color

