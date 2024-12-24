import numpy as np
from stl import mesh
from scipy.spatial import cKDTree


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
    """
    For each triangle in the numpy stl mesh, store the distance from each point in the
    triangle to the track. Returns Nx3, where N is the number of triangles, and each
    of the three elements are the distances from each triangle vertex to the track.
    """
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


def construct_new_triangles(np_mesh, subdivide_mask, n_sub=1):
    new_triangles = []  # will be a list of shape (n_new, 3, 3)
    is_subdivided_mask = [] # false for triangles that are the same as original, true for subdivided ones

    n_vectors = len(np_mesh.vectors)
    i_to_print = n_vectors // 100
    for i in range(n_vectors):
        percent_done = i / n_vectors * 100
        if i % i_to_print == 0:
            print(f"  creating new triangles ... {round(percent_done, 0)}% complete")

        triangle = np_mesh.vectors[i]  # shape (3,)

        if not subdivide_mask[i]:
            # Keep original triangle
            new_triangles.append(triangle)
            is_subdivided_mask.append(False)
        else:
            new_sub_tris = subdivide_triangle(triangle, n_sub)
            for tri in new_sub_tris:
                new_triangles.append(tri)
                is_subdivided_mask.append(True)

    return np.array(new_triangles), np.array(is_subdivided_mask)  # shape (N_new, 3, 3)


def subdivide_triangle(triangle, n_sub=1):
    """
    Recursively subdivide a triangle n_sub times. Each subdivision splits one triangle into 4.

    Parameters
    ----------
    triangle : array-like of shape (3, 3)
        The coordinates of the triangle's vertices.
    n_sub : int
        Number of subdivisions. For example:
          - n_sub=0 => return the original triangle (no subdivision).
          - n_sub=1 => return 4 sub-triangles.
          - n_sub=2 => each of those 4 is again subdivided => 16, etc.

    Returns
    -------
    list of np.ndarray
        A list of sub-triangles. Each sub-triangle is a (3, 3) NumPy array of vertices.
    """
    # If no more subdivisions left, just return the triangle itself
    if n_sub <= 0:
        return [triangle]

    A, B, C = triangle

    # Compute midpoints
    AB_mid = 0.5 * (A + B)
    BC_mid = 0.5 * (B + C)
    CA_mid = 0.5 * (C + A)

    # Perform one level of subdivision => 4 new triangles
    new_triangles = [
        np.array([A, AB_mid, CA_mid]),
        np.array([B, BC_mid, AB_mid]),
        np.array([C, CA_mid, BC_mid]),
        np.array([AB_mid, BC_mid, CA_mid])
    ]

    # If n_sub=1, we're done. Otherwise, subdivide each of these further.
    if n_sub == 1:
        return new_triangles
    else:
        final = []
        for tri in new_triangles:
            final.extend(subdivide_triangle(tri, n_sub - 1))
        return final


def upsample_along_track(np_mesh, track_kdtree, dist_to_track, n_sub=1):
    print(f"Refining mesh near the GPX track...")
    subdivide_mask = create_subdivide_mask(np_mesh, track_kdtree, dist_to_track)
    new_triangles, is_subdivided_mask = construct_new_triangles(np_mesh, subdivide_mask, n_sub)
    refined_mesh = mesh.Mesh(np.zeros(len(new_triangles), mesh.Mesh.dtype))
    refined_mesh.vectors = new_triangles
    return refined_mesh, is_subdivided_mask


def cut_along_track(np_mesh, track_kdtree, is_upsampled_mask, cut_radius_mm):
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

        if not is_upsampled_mask[i]:
            # non-upsampled triangles won't participate in the cut; they've been determined
            # earlier to be too far from the track
            continue

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
    refined_mesh, is_subdivided_mask = upsample_along_track(np_mesh, track_kdtree, dist_to_refine, n_sub=2)
    cut_mesh = cut_along_track(refined_mesh, track_kdtree, is_subdivided_mask, cut_radius_mm)

    # 4. Save the modified mesh back to file
    cut_mesh.save(out_path)
