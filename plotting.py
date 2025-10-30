import pyvista as pv
from stl import mesh
import numpy as np

def plot_mesh_with_track(stl_mesh: mesh.Mesh, track_points: np.ndarray):
    """

    :param stl_mesh:
    :param track_points:
    :return:
    """
    pv_mesh = stlmesh_to_pvmesh(stl_mesh)

    # create plotter
    plotter = pv.Plotter()
    plotter.add_mesh(pv_mesh, color="lightgray", opacity=0.5)

    # sample track data: assume numpy arrays x, y, z
    plotter.add_lines(track_points, color="red", width=3)  # or add_points

    plotter.show()  # opens interactive window (rotate/zoom/pan)


def stlmesh_to_pvmesh(stl_mesh):
    """
    Convert an `stl.mesh.Mesh` object (numpy-stl) to a PyVista PolyData mesh.

    Parameters:
        stl_mesh: instance of stl.mesh.Mesh

    Returns:
        pv.PolyData
    """
    verts = stl_mesh.vectors.reshape(-1, 3)

    faces = np.hstack(
        np.column_stack([
            np.full(len(stl_mesh.vectors), 3),
            np.arange(len(verts)).reshape(-1, 3)
        ])
    )

    return pv.PolyData(verts, faces)

def plot_mesh_with_track_w_sliders(stl_mesh, track_points: np.ndarray):
    pv_mesh = stlmesh_to_pvmesh(stl_mesh)

    plotter = pv.Plotter()
    plotter.add_mesh(pv_mesh, color="lightgray", opacity=0.5)

    # store original track so we can reapply shifts
    base = pv.PolyData(track_points)
    actor = plotter.add_mesh(base, color="red", line_width=3, name="track_actor")

    shift = {"x": 0.0, "y": 0.0, "z": 0.0}

    def update():
        tx, ty, tz = shift["x"], shift["y"], shift["z"]
        new_ = base.copy()
        new_.points = new_.points + np.array([tx, ty, tz])
        plotter.remove_actor("track_actor")
        plotter.add_mesh(new_, color="red", line_width=3, name="track_actor")
        plotter.render()

    plotter.add_slider_widget(lambda v: _set_and_update(shift, "x", v, update),
                              rng=[-5.0, 5.0], title="Shift X",
                              pointa=(0.02, 0.90),
                              pointb=(0.32, 0.90))
    plotter.add_slider_widget(lambda v: _set_and_update(shift, "y", v, update),
                              rng=[-5.0, 5.0], title="Shift Y",
                              pointa=(0.02, 0.83),
                              pointb=(0.32, 0.83))
    plotter.add_slider_widget(lambda v: _set_and_update(shift, "z", v, update),
                              rng=[-5.0, 5.0], title="Shift Z",
                              pointa=(0.02, 0.76),
                              pointb=(0.32, 0.76))

    plotter.show()

    final_shift = np.array([shift["x"], shift["y"], shift["z"]])
    track_points += final_shift


def _set_and_update(shift, axis, value, update_fn):
    shift[axis] = value
    update_fn()



