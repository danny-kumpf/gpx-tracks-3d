import json
import os
import numpy as np
from stl import mesh
from scipy.spatial import cKDTree

from gpx_coordinate_transforms import convert_gpx_track_to_stl_coordinates
import auto_cut
from plotting import plot_mesh_with_track

# In fusion 360, we (use to) have to:
# - load the STL (import mesh)
# - load the snapped track that this script creates (utilities, add-ins, scripts and add-ins, ImportSplineCSV)
#   (this script aligns the coordinate systems; no manual movement in f360 needed)
# - create a plane perpendicular to the track (solid, construct, plane along path)
# - in that plane, create a 1.2mm circle cross-section centered on the path
# - create a swept SURFACE (cant do solid bc path intersects itself)
# - create end caps: surface, patch, click the end caps of the new surface pipe
# - fill the surface: solid, create, boundary fill
# - convert the surface to mesh: mesh, tesselate (can do low quality)
# - merge the two meshes as a cut: mesh, modify, combine (has to cook for a few min)
# (Now, autocut.py does the job that fusion 360 used to do)

def load_from_json(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data


def clean_up_track(enu_points, hike_type):
    """
    Currently performs the following clean-up steps:
    - For hike_type OUT_AND_BACK, removes points after the highest elevation
    -
    """

    if hike_type == "OUT_AND_BACK":
        print("OUT_AND_BACK, removing points after highest elevation")
        enu_points = remove_after_highest_elev(enu_points)

    # Dont really need to downsample with new method, but could
    #enu_points = downsample_to(enu_points, 100)
    #plot_3d(enu_points)

    print("Removing identical points...")
    enu_points = filter_identical(enu_points)

    if hike_type == "LOOP":
        print("LOOP, duplicating first point at end")
        # Add a duplicate of the first point at the end, so that we get a closed
        # loop. Saves some steps in fusion360.
        enu_points = np.vstack((enu_points, enu_points[0]))
    #plot_3d(enu_points)

    return enu_points


def remove_after_highest_elev(points):
    index_max = max(range(len(points)), key=lambda idx: points[idx][2])
    print(f"Index of maximum elevation: {index_max} (elevation = {points[index_max][2]} m")
    return points[:index_max]


def filter_identical(enu_points):
    """
    Return only the unique points
    :param enu_points (np.ndarray):
    :param min_dist_m (float):
    :return:
    """
    _, idx = np.unique(enu_points, axis=0, return_index=True)
    return enu_points[np.sort(idx)]


if __name__ == "__main__":
    ## Change the json path to change inputs
    json_inputs = load_from_json(os.path.join("data", "mr_dlp_ironman", "config.json"))

    ## -- Dont need to change below here
    in_stl_mesh = mesh.Mesh.from_file(json_inputs["in_stl_path"])

    # allow the user to specify one or multiple paths
    if not isinstance(json_inputs["gpx_path"], list):
        json_inputs["gpx_path"] = [json_inputs["gpx_path"]]

    track_kdtrees = []
    for gpx_path in json_inputs["gpx_path"]:
        track_points = convert_gpx_track_to_stl_coordinates(
            gpx_path,
            in_stl_mesh,
            json_inputs["box_upper_right"],
            json_inputs["box_lower_left"]
        )
        track_points = clean_up_track(track_points, json_inputs["hike_type"])
        plot_mesh_with_track(in_stl_mesh, track_points)

        track_points = auto_cut.densify_track_linear(track_points, step=0.01)
        track_kdtrees.append(cKDTree(track_points))

    cut_radius_mm = 0.9
    dist_to_refine = 1.2 * cut_radius_mm  # upsample triangles within this distance from the track

    # subdivide the triangles that are near the track, then perform the cut
    mesh_upsampled_near_track, is_upsampled_mask = auto_cut.upsample_along_track(
        in_stl_mesh, track_kdtrees, dist_to_refine, n_sub=2)
    cut_mesh = auto_cut.cut_along_track(mesh_upsampled_near_track, track_kdtrees, is_upsampled_mask, cut_radius_mm)

    # 4. Save the modified mesh back to file
    out_path = json_inputs["out_stl_path"]
    cut_mesh.save(out_path)
    print(f"Saved new STL file: {out_path}")
