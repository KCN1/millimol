import numpy as np
from numpy.typing import NDArray
from pyvista import Plotter, PolyData, ColorLike


PLOTTER_PARAMS = dict(
    rgb=True,
    render_points_as_spheres=True,
    render_lines_as_tubes=True,
)


def draw_molecule(
    mol_name: str,
    coords: NDArray[np.floating],
    pair_indices: NDArray[np.integer],
    point_indices: NDArray[np.integer],
    color_nums: NDArray[np.integer],
    color_map: list[ColorLike],
    background_color: ColorLike, # e.g. "red", [255, 0, 0], "#FF0000"
    *,
    point_size=6,
    line_width=3,
    lighting=True,
) -> None:
    """
    Main function to draw molecule (MoleculeBlueprint) with Pyvista.
    See also: mol_blueprint.MoleculeBlueprint.__doc__
    
    coords: coordinates of points ("real" atoms and "virtual" middle points),
    pair_indices: indices of bound pairs of points to draw them as lines,
    point_indices: indices of unbound points to draw them as spheres*,
    color_nums: color indices in color_map**,
    color_map: list of colors
    
    * To draw all atoms as spheres, include all atoms in point_indices.
    ** The line color is determined by its first point (a "real" atom).
    """
    mesh_params = dict(
        point_size=point_size,
        line_width=line_width,
        lighting=lighting,
    )
    
    cmap = np.asarray(color_map, dtype=np.uint8)
    # array of lines [2, u, v, ...]
    padding = np.full((pair_indices.shape[0], 1), 2)
    edges = np.hstack((padding, pair_indices)).ravel()
    # create mesh (PolyData)
    mesh = PolyData(coords, lines=edges)
    # create window
    plotter = Plotter()
    plotter.background_color = background_color
    # draw lines colored by 1st point
    plotter.add_mesh(
        mesh,
        scalars=cmap[color_nums[pair_indices[:, 0]]],
        **PLOTTER_PARAMS,
        **mesh_params
    )
    # draw points
    if len(point_indices):
        plotter.add_points(
            mesh.points[point_indices],
            scalars=cmap[color_nums[point_indices]],
            **PLOTTER_PARAMS,
            **mesh_params
        )
    plotter.show(title=mol_name)

