import numpy as np
from stl import mesh
from scipy.spatial import cKDTree
import pyvista as pv


def densify_track_linear(track_points, step=0.01):
    """
    track_points: (N,3) array of existing points
    step: desired spacing for new points (float, mm)

    Returns: (M,3) array with new points inserted
    """
    dense_pts = []
    for i in range(len(track_points)-1):
        start = track_points[i]
        end   = track_points[i+1]
        seg_vec = end - start
        seg_len = np.linalg.norm(seg_vec)
        #print(f"Segment length: {seg_len} mm")
        if seg_len == 0:
            # Degenerate segment, skip or just add 'start'
            dense_pts.append(start)
            continue

        # Number of subdivisions on this segment
        n_sub = int(seg_len // step)
        # Parametric increments
        for s in range(n_sub):
            t = s * step / seg_len
            point = start + t * seg_vec
            dense_pts.append(point)

    return np.array(dense_pts)


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


def polynomial_drop(distance, radius):
    """
    distance : float
        Distance from the track.
    radius : float
        Radius of the half-circle-ish cross section with flat derivatives
        at the ends

    Returns
    -------
    float
        The z-drop. Max drop is 'radius' when distance=0. Min drop is 0 when
        distance > radius. d(z-drop)/d(distance) == 0 when distance == radius.
    """
    if distance >= radius:
        return 0.0
    n = 5
    alpha = 3
    cut_amt = radius * pow(1 - pow(distance / radius, n), alpha)
    return cut_amt


def is_near_track(point, track_kdtree, threshold):
    dist, _ = track_kdtree.query(point)
    return dist < threshold


def get_distances_from_track(np_mesh, track_kdtree):
    distances = np.zeros((len(np_mesh.vectors), 3), dtype=bool)

    n_vectors = len(np_mesh.vectors)
    i_to_print = n_vectors // 100
    for i in range(n_vectors):
        percent_done = i / n_vectors * 100
        if i % i_to_print == 0:
            print(f"  finding distances of vertices to track ... {round(percent_done, 0)}% complete")

        tri_verts = np_mesh.vectors[i]  # shape (3,3)
        for j in range(3):
            distances[i][j] = track_kdtree.query(tri_verts[j])

    return distances

def create_subdivide_mask(np_mesh, track_kdtree, dist_to_track):
    """
    Return a mask that is true for triangles that have any vertex within
    dist_to_track of the track, else false.
    :param np_mesh:
    :param track_kdtree:
    :param dist_to_track:
    :return:
    """
    # Suppose we have a function that sets up `subdivide_mask`
    # with True/False per triangle:
    subdivide_mask = np.zeros(len(np_mesh.vectors), dtype=bool)
    #distances_to_track = get_distances_from_track(np_mesh, track_kdtree)

    n_vectors = len(np_mesh.vectors)
    i_to_print = n_vectors // 100
    for i in range(n_vectors):
        percent_done = i / n_vectors * 100
        if i % i_to_print == 0:
            print(f"  finding points near track ... {round(percent_done, 0)}% complete")

        tri_verts = np_mesh.vectors[i]  # shape (3,3)

        # If ANY vertex is within threshold => subdivide
        if any(is_near_track(v, track_kdtree, dist_to_track) for v in tri_verts):
            subdivide_mask[i] = True

    return subdivide_mask


def construct_new_triangles(np_mesh, subdivide_mask):
    new_triangles = []  # will be a list of shape (n_new, 3, 3)

    n_vectors = len(np_mesh.vectors)
    i_to_print = n_vectors // 100
    for i in range(n_vectors):
        percent_done = i / n_vectors * 100
        if i % i_to_print == 0:
            print(f"  creating new triangles ... {round(percent_done, 0)}% complete")

        A, B, C = np_mesh.vectors[i]  # shape (3,)

        if not subdivide_mask[i]:
            # Keep original triangle
            new_triangles.append([A, B, C])
        else:
            # Subdivide
            AB_mid = 0.5*(A + B)
            BC_mid = 0.5*(B + C)
            CA_mid = 0.5*(C + A)

            # Add the 4 new triangles
            new_triangles.append([A, AB_mid, CA_mid])
            new_triangles.append([B, BC_mid, AB_mid])
            new_triangles.append([C, CA_mid, BC_mid])
            new_triangles.append([AB_mid, BC_mid, CA_mid])

    return np.array(new_triangles)  # shape (N_new, 3, 3)


def refine_along_track(np_mesh, track_kdtree, dist_to_track):
    print(f"Refining mesh near the GPX track...")
    subdivide_mask = create_subdivide_mask(np_mesh, track_kdtree, dist_to_track)
    new_triangles = construct_new_triangles(np_mesh, subdivide_mask)
    refined_mesh = mesh.Mesh(np.zeros(len(new_triangles), mesh.Mesh.dtype))
    refined_mesh.vectors = new_triangles
    return refined_mesh


def cut_along_track(np_mesh, track_kdtree, cut_radius_mm):
    """
    Perform a cut along the track: reduce the z coordinate of all STL
    vertices within cut_radius_mm of the track so that they form a circular
    cut cross-section
    """
    n_vectors = len(np_mesh.vectors)
    i_to_print = n_vectors // 100
    for i in range(len(np_mesh.vectors)):
        percent_done = i / n_vectors * 100
        if i % i_to_print == 0:
            print(f"Cutting out the track ... {round(percent_done, 0)}% complete")

        # Each `vectors[i]` is a triangle with 3 vertices
        for j in range(3):
            vertex = np_mesh.vectors[i][j]  # [x, y, z]
            dist_to_track, idx = track_kdtree.query(vertex)
            vertex[2] -= polynomial_drop(dist_to_track, cut_radius_mm)

    # 3. Optional: Update normals if you want them recalculated
    np_mesh.update_normals()
    return np_mesh


if __name__ == "__main__":
    track_csv_path = r"D:\Projects\Python\GpxTracks3d\data\jess_boulder_5_peaks\jess_boulder_5_peaks_AUTOCUT.csv"
    stl_path = r'D:\Projects\Python\GpxTracks3d\data\jess_boulder_5_peaks\10m_-105.30_39.99_tile_1_1.STL'
    out_path = r'D:\Projects\Python\GpxTracks3d\data\jess_boulder_5_peaks\jess_boulder_5_peaks_autocut.stl'

    # Load the track into a KDTree for easy nearest-neighbor lookups:
    track_points = np.loadtxt(track_csv_path, delimiter=",")
    track_points = densify_track_linear(track_points, step=0.01)
    track_kdtree = cKDTree(track_points)

    cut_radius_mm = 0.6
    dist_to_refine = 1.5 * cut_radius_mm # refine triangles within this distance

    # Load the STL, subdivide the triangles that are near the track, then
    # perform the cut
    np_mesh = mesh.Mesh.from_file(stl_path)
    refined_mesh = refine_along_track(np_mesh, track_kdtree, dist_to_refine)
    refined_mesh = refine_along_track(refined_mesh, track_kdtree, dist_to_refine)
    cut_mesh = cut_along_track(refined_mesh, track_kdtree, cut_radius_mm)

    # 4. Save the modified mesh back to file
    cut_mesh.save(out_path)
