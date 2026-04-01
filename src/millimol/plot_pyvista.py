import numpy as np
from numpy.typing import NDArray
from pyvista import Plotter, PolyData, ColorLike


PLOTTER_PARAMS = dict(
    categories=True,
    interpolate_before_map=False,
    show_scalar_bar=False,
    render_points_as_spheres=True,
    render_lines_as_tubes=True,
)


def draw_molecule(
    mol_name: str,
    coords: NDArray[np.floating],
    pair_indices: NDArray[np.integer],
    point_indices: NDArray[np.integer],
    color_nums: NDArray[np.integer],
    color_map: list,
    background_color: ColorLike, # e.g. "red", [255, 0, 0], "#FF0000"
    *,
    point_size=6,
    line_width=3,
    lighting=True,
) -> None:
    """
    Main function to draw molecule (MoleculeBlueprint) with Pyvista.
    
    coords: coordinates of points (atoms and virtual middle points),
    pair_indices: indices of bound pairs of points to draw lines,
    point_indices: indices of unbound points to draw points as spheres,
    color_nums: color indices in color_map
    color_map: colormap to use, e.g. {1: "red", 2: [255, 0, 0], 3: "#FF0000"}
    To draw all atoms as spheres, point_indices should contain array of atoms.
    """
    MESH_PARAMS = dict(
        cmap=color_map,
        clim=(-0.5, len(color_map) - 0.5),
        n_colors=len(color_map),
        point_size=point_size,
        line_width=line_width,
        lighting=lighting,
    )
    # array of lines [2, u, v, ...]
    # 2 is the number of line ends
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
        scalars=color_nums[pair_indices[:, 0]],
        **PLOTTER_PARAMS,
        **MESH_PARAMS
    )
    # draw points
    if len(point_indices):
        plotter.add_points(
            mesh.points[point_indices],
            scalars=color_nums[point_indices],
            **PLOTTER_PARAMS,
            **MESH_PARAMS
        )
    plotter.show(title=mol_name)

