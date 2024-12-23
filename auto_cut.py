import numpy as np
from stl import mesh
from scipy.spatial import cKDTree
import pyvista as pv

def point_line_segment_distance_3d(p, a, b):
    ap = p - a
    ab = b - a
    ab_norm_sq = np.dot(ab, ab)
    if ab_norm_sq == 0.0:
        return np.linalg.norm(ap)
    t = np.dot(ap, ab) / ab_norm_sq
    if t < 0.0:
        closest = a
    elif t > 1.0:
        closest = b
    else:
        closest = a + t * ab
    return np.linalg.norm(p - closest)

def point_polyline_distance_3d(p, polyline):
    min_dist = float('inf')
    for i in range(len(polyline) - 1):
        dist_seg = point_line_segment_distance_3d(
            p,
            np.array(polyline[i]),
            np.array(polyline[i+1])
        )
        if dist_seg < min_dist:
            min_dist = dist_seg
    return min_dist

def smootherstep(x):
    """ Smoother step in [0,1]. """
    return 6*x**5 - 15*x**4 + 10*x**3

def smooth_falloff(distance, max_dist, max_drop):
    """
    Returns the z-drop for a given distance, using smoothstep.

    distance : float
        Distance from the track.
    max_dist : float
        Distance beyond which there is no drop.
    max_drop : float
        Max drop at distance=0.
    """
    if distance >= max_dist:
        return 0.0
    # Scale distance into [0, 1]
    x = distance / max_dist
    # We can pick either smoothstep or smootherstep
    # Here, let's use standard smoothstep:
    return max_drop * (1.0 - smootherstep(x))

def half_circle_drop(distance, radius):
    """
    distance : float
        Distance from the track.
    radius : float
        Radius of the half-circle cross section.

    Returns
    -------
    float
        The z-drop. Max drop is 'radius' when distance=0.
    """
    if distance >= radius:
        return 0.0
    return np.sqrt(radius**2 - distance**2)

def locally_refine(mesh, track_points_kdtree):
    # Create an array of distances for all mesh vertices
    distances = np.zeros(mesh.n_points, dtype=np.float64)
    for i in range(mesh.n_points):
        p = mesh.points[i]  # (x, y, z)
        distances[i], _ = track_points_kdtree.query(p)

    # Attach distances to the mesh as a "point_data" array
    mesh.point_data["distances"] = distances

    # 3. Extract (mask) the region near the track
    #    We'll threshold to keep points where distance <= max_dist
    submesh_near = mesh.threshold(
        value=max_dist,
        scalars="distances",
        invert=False  # False => keep those <= value
    )
    # `submesh_near` now contains only the cells whose *all* points
    # or partial points are within max_dist, depending on thresholding logic.
    # By default, thresholding in PyVista tries to keep entire cells
    # for which the scalar is within the threshold. You can tweak
    # how partial cells are handled with `contouring`, etc.

    # 4. Subdivide only that sub-mesh
    #    (Pick your subdivision method; e.g. 'loop')
    #refined_submesh = submesh_near.subdivide_tetra()
    # Increase `nsub` to subdivide more aggressively.

    # 5. Extract the "far" region (i.e., outside max_dist)
    #    so we can merge it back.  We'll do an inverted threshold
    submesh_far = mesh.threshold(
        value=max_dist,
        scalars="distances",
        invert=True  # True => keep those > value
    )

    # 6. Combine the refined submesh with the unrefined remainder
    combined = submesh_far + submesh_near

    # 7. Save the merged result
    combined.save("locally_refined.stl")

if __name__ == "__main__":
    track_csv_path = r"D:\Projects\Python\GpxTracks3d\data\jess_boulder_5_peaks\jess_boulder_5_peaks_AUTOCUT.csv"
    stl_path = r'D:\Projects\Python\GpxTracks3d\data\loyalsock\small_version\10m_-76.55_41.48_tile_1_1.STL'
    out_path = r'D:\Projects\Python\GpxTracks3d\data\jess_boulder_5_peaks\jess_boulder_5_peaks_autocut.stl'

    # Your track in 3D, for example:
    track_points = np.loadtxt(track_csv_path, delimiter=",")
    kdtree = cKDTree(track_points)

    max_dist = 0.6   # Distance threshold
    z_drop   = 1    # Amount to drop the Z-coordinate

    # 0. up-sample the STL:
    print(f"Refining STL...")
    pv_mesh = pv.read(stl_path)
    below, above = pv_mesh.clip(return_clipped=True, normal='z', origin=(0, 0, 1.0))
    refined = above.subdivide(nsub=2, subfilter='loop')  # subdivide each triangle
    # refined.save("refined.stl")
    #locally_refine(pv_mesh, kdtree)
    remade = below + refined
    remade.save("test.stl")
    print(f"Saved refined STL")

    # # 1. Read the STL file
    # my_mesh = mesh.Mesh.from_file("refined.stl")
    #
    # # 2. Modify Z of vertices that are within `max_dist` of the track
    # n_vectors = len(my_mesh.vectors)
    # i_to_print = n_vectors // 100
    # for i in range(len(my_mesh.vectors)):
    #     percent_done = i / n_vectors * 100
    #     if i % i_to_print == 0:
    #         print(f"processing vertices... {round(percent_done, 0)}% complete")
    #     # Each `vectors[i]` is a triangle with 3 vertices
    #     for j in range(3):
    #         vertex = my_mesh.vectors[i][j]  # [x, y, z]
    #         dist_to_track, idx = kdtree.query(vertex)
    #         # dist_to_track = point_polyline_distance_3d(vertex, track_points)
    #         vertex[2] -= half_circle_drop(dist_to_track, max_dist)  # reduce Z
    #         # if dist_to_track <= max_dist:
    #         #     vertex[2] -= smooth_falloff(dist_to_track, max_dist, z_drop)  # reduce Z
    #
    # # 3. Optional: Update normals if you want them recalculated
    # my_mesh.update_normals()
    #
    # # 4. Save the modified mesh back to file
    # my_mesh.save(out_path)
