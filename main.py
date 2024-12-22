import gpxpy
import csv
import math
import json
import os
import pyproj
from pyproj import Transformer
import numpy as np
from stl import mesh
from scipy.spatial import KDTree
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import requests

def read_stl_vertices(stl_file_path):
    """
    Reads an STL file and extracts unique vertices.

    Parameters:
        stl_file_path (str): Path to the STL file.

    Returns:
        vertices (numpy.ndarray): Array of unique vertices (Nx3).
    """
    # Load the STL file
    stl_mesh = mesh.Mesh.from_file(stl_file_path)
    # Concatenate all vertices from the triangles
    all_vertices = np.concatenate((stl_mesh.v0, stl_mesh.v1, stl_mesh.v2))
    # Get unique vertices
    vertices = np.unique(all_vertices, axis=0)
    return vertices

def find_nearest_vertices(points, vertices):
    """
    Finds the closest mesh vertex to each point.

    Parameters:
        points (numpy.ndarray): Array of track points (Mx3).
        vertices (numpy.ndarray): Array of mesh vertices (Nx3).

    Returns:
        nearest_vertices (numpy.ndarray): Array of closest vertices (Mx3).
    """
    # Build KD-tree from the mesh vertices
    tree = KDTree(vertices)
    # Query the tree for nearest neighbors
    distances, indices = tree.query(points)
    nearest_vertices = vertices[indices]

    # If the original track (before the snap-to-vertices) was some fixed offset
    # above/below the STL surface, then we'd like to first raise/lower the whole
    # track until it's closer to the surface, and then re-try the snap. We don't
    # want to snap with the track too far away from the surface, because we
    # could snap to the wrong points
    avg_z_diff_mmm = np.mean(points[:,2] - nearest_vertices[:,2])
    print(f"Average height of track above STL: {avg_z_diff_mmm} mmm")
    if abs(avg_z_diff_mmm) > 0.01:
        points[:, 2] -= avg_z_diff_mmm
        return find_nearest_vertices(points, vertices)

    return nearest_vertices

def read_gpx_file(gpx_file_path):
    """
    Reads a GPX file and returns a GPX object.

    Parameters:
        gpx_file_path (str): Path to the GPX file.

    Returns:
        gpx (gpxpy.gpx.GPX): Parsed GPX object.
    """
    with open(gpx_file_path, 'r') as gpx_file:
        gpx = gpxpy.parse(gpx_file)
    return gpx

def get_first_track_point(gpx):
    """
    Extracts the first track point to use as the origin.

    Parameters:
        gpx (gpxpy.gpx.GPX): Parsed GPX object.

    Returns:
        origin (tuple): Tuple containing (lon_deg, lat_deg, altitude_m).
    """
    for track in gpx.tracks:
        for segment in track.segments:
            if segment.points:
                first_point = segment.points[0]
                lon_deg = first_point.longitude  # Degrees
                lat_deg = first_point.latitude    # Degrees
                altitude_m = first_point.elevation if first_point.elevation is not None else 0.0  # Meters
                return (lon_deg, lat_deg, altitude_m)
    # If no tracks with points are found, raise an error
    raise ValueError("No track points found in the GPX file.")

def create_transformer_lla_to_ecef():
    """
    Creates and returns the LLA (Latitude, Longitude, Altitude) to ECEF (Earth-Centered, Earth-Fixed) transformer.

    Returns:
        transformer (pyproj.Transformer): Transformer object for LLA to ECEF conversion.
    """
    # Initialize the transformer once for efficiency
    transformer = Transformer.from_crs("epsg:4979", "epsg:4978", always_xy=True)
    return transformer

def create_transformer_ecef_to_lla():
    """
    Creates and returns the ECEF (Earth-Centered, Earth-Fixed) to LLA (Latitude, Longitude, Altitude) transformer.

    Returns:
        transformer (pyproj.Transformer): Transformer object for ECEF to LLA conversion.
    """
    # Initialize the transformer once for efficiency
    transformer = Transformer.from_crs("epsg:4978", "epsg:4979", always_xy=True)
    return transformer

def precompute_origin_parameters(lat_deg, lon_deg):
    """
    Precomputes sine and cosine of the origin's latitude and longitude in radians.

    Parameters:
        lat_deg (float): Latitude of the origin in degrees.
        lon_deg (float): Longitude of the origin in degrees.

    Returns:
        sin_lat (float): Sine of the origin's latitude in radians.
        cos_lat (float): Cosine of the origin's latitude in radians.
        sin_lon (float): Sine of the origin's longitude in radians.
        cos_lon (float): Cosine of the origin's longitude in radians.
    """
    # Convert latitude and longitude to radians
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    # Precompute sines and cosines
    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)
    return sin_lat, cos_lat, sin_lon, cos_lon

def extract_points(gpx):
    """
    Extracts all points from the GPX tracks.

    Parameters:
        gpx (gpxpy.gpx.GPX): Parsed GPX object.

    Returns:
        points (np.ndarray): Nx3 lon_deg, lat_deg, altitude_m
    """
    points = []
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                lon_deg = point.longitude  # Degrees
                lat_deg = point.latitude    # Degrees
                altitude_m = point.elevation if point.elevation is not None else 0.0  # Meters
                points.append((lon_deg, lat_deg, altitude_m))

    return np.array(points)


def remove_after_highest_elev(points):
    index_max = max(range(len(points)), key=lambda idx: points[idx][2])
    print(f"Index of maximum elevation: {index_max} (elevation = {points[index_max][2]} m")
    return points[:index_max]


def convert_points_to_enu(points, origin, transformer):
    """
    Converts a list of geodetic coordinates to ENU coordinates relative to the origin.

    Parameters:
        points (np.ndarray): Nx3 lon_deg, lat_deg, altitude_m.
        origin (tuple): Tuple containing (lon_deg, lat_deg, altitude_m).
        transformer (pyproj.Transformer): Transformer for LLA to ECEF conversion.

    Returns:
        enu_points (np.ndarraay): Nx3 east_m, north_m, up_m
    """
    # Unpack origin coordinates
    lon_origin_deg, lat_origin_deg, altitude_origin_m = origin

    # Precompute sine and cosine values for the origin's latitude and longitude
    sin_lat0, cos_lat0, sin_lon0, cos_lon0 = precompute_origin_parameters(
        lat_origin_deg, lon_origin_deg
    )

    # Convert origin to ECEF coordinates once
    x0_m, y0_m, z0_m = transformer.transform(
        lon_origin_deg, lat_origin_deg, altitude_origin_m
    )

    enu_points = np.ndarray(points.shape)
    for idx in range(points.shape[0]):
        lon_deg, lat_deg, altitude_m = points[idx]

        # Convert current point to ECEF coordinates
        x_m, y_m, z_m = transformer.transform(
            lon_deg, lat_deg, altitude_m
        )
        # Compute differences in ECEF coordinates (meters)
        dx_m = x_m - x0_m
        dy_m = y_m - y0_m
        dz_m = z_m - z0_m

        # Compute ENU coordinates (meters)
        east_m = -sin_lon0 * dx_m + cos_lon0 * dy_m
        north_m = (
                -sin_lat0 * cos_lon0 * dx_m
                - sin_lat0 * sin_lon0 * dy_m
                + cos_lat0 * dz_m
        )
        up_m = (
                cos_lat0 * cos_lon0 * dx_m
                + cos_lat0 * sin_lon0 * dy_m
                + sin_lat0 * dz_m
        )

        enu_points[idx] = np.array([east_m, north_m, up_m])

    return enu_points


def filter_min_distance(enu_points, min_dist_m):
    """
    Remove a point if it is not at least min_dist_m away from the previous point
    :param enu_points (np.ndarray):
    :param min_dist_m (float):
    :return:
    """
    if len(enu_points) == 0:
        return enu_points

    # Initialize the list with the first point
    filtered_points = [enu_points[0]]

    for point in enu_points[1:]:
        last_point = filtered_points[-1]
        if distance(point, last_point) >= min_dist_m:
            filtered_points.append(point)

    return np.array(filtered_points)


def filter_identical(enu_points):
    """
    Remove a point if it is not at least min_dist_m away from the previous point
    :param enu_points (np.ndarray):
    :param min_dist_m (float):
    :return:
    """
    _, idx = np.unique(enu_points, axis=0, return_index=True)
    return enu_points[np.sort(idx)]
    # if len(enu_points) == 0:
    #     return enu_points
    #
    # # Initialize the list with the first point
    # filtered_points = [enu_points[0]]
    #
    # for point in enu_points[1:]:
    #     last_point = filtered_points[-1]
    #     if not np.array_equal(point, last_point):
    #         filtered_points.append(point)
    #
    # return np.array(filtered_points)


def write_csv(csv_file_path, enu_points):
    """
    Writes the list of ENU coordinates to a CSV file.

    Parameters:
        csv_file_path (str): Path to the output CSV file.
        enu_points (list): List of tuples containing (east_m, north_m, up_m).
    """
    with open(csv_file_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header with units specified
        csv_writer.writerow(['east_m', 'north_m', 'up_m'])
        # Write data rows
        for east_m, north_m, up_m in enu_points:
            csv_writer.writerow([east_m, north_m, up_m])


def distance(point_1, point_2):
    return np.linalg.norm(point_1 - point_2)


FUSION_360_MODEL_MM_PER_CSV = 10


def snap_gpx_to_stl(gpx_file_path, stl_path, box_upper_right, box_lower_left, hike_type):
    """
    Main function to convert a GPX file to a CSV file with ENU coordinates relative to the first track point.

    Parameters:
        gpx_file_path (str): Path to the input GPX file.
        csv_file_path (str): Path to the output CSV file.
    """
    # Read the GPX file
    gpx = read_gpx_file(gpx_file_path)

    # Read the STL file
    vertices = read_stl_vertices(stl_path)
    x_range_mmm = (vertices[:, 0].min(), vertices[:, 0].max())
    y_range_mmm = (vertices[:, 1].min(), vertices[:, 1].max())
    z_range_mmm = (vertices[:, 2].min(), vertices[:, 2].max())
    z_max_mmm = z_range_mmm[1]

    print(f"STL Number of vertices: {len(vertices)}")
    print(f"STL X range: {x_range_mmm[0]} to {x_range_mmm[1]}")
    print(f"STL Y range: {y_range_mmm[0]} to {y_range_mmm[1]}")
    print(f"STL Z range: {z_range_mmm[0]} to {z_range_mmm[1]}")

    # Create the LLA to ECEF transformer once
    transformer = create_transformer_lla_to_ecef()

    box_enu = convert_points_to_enu(np.array([box_upper_right, box_lower_left]),
                                    origin=get_first_track_point(gpx),
                                    transformer=transformer)
    box_diag_m = distance(box_enu[0], box_enu[1])
    print(f"Box diagonal meters = {box_diag_m}")

    model_width_mmm = x_range_mmm[1] - x_range_mmm[0]
    model_height_mmm = y_range_mmm[1] - y_range_mmm[0]
    model_diag_mm = math.sqrt(model_width_mmm ** 2 + model_height_mmm ** 2)
    mmm_per_m = model_diag_mm / box_diag_m  # model mm per meter
    print(f"Meters-to-model-mm scale factor: {mmm_per_m}")

    # Get the origin point (lon, lat, alt of lower left corner of model)

    # this elevation isn't super accurate, but we correct for it later by
    # iteratively snapping to STL vertices
    origin_elev = get_elevation_open_elevation(box_lower_left[1], box_lower_left[0])

    print(f"Elevation of bottom-left corner ground surface: ~{origin_elev}")
    _, origin_top_pt = highest_z_near_origin(vertices)
    origin_elev -= origin_top_pt[2] / mmm_per_m
    origin = (box_lower_left[0],
              box_lower_left[1],
              origin_elev)
    print(f"Elevation of bottom-left model corner (m): ~{origin[2]}")

    # Extract all points from the GPX file
    points = extract_points(gpx)
    if hike_type == "OUT_AND_BACK":
        points = remove_after_highest_elev(points)

    # Convert points to ENU coordinates relative to the origin
    enu_points = convert_points_to_enu(points, origin, transformer)
    #plot_3d(enu_points)
    enu_points = filter_min_distance(enu_points, 30)
    enu_points *= mmm_per_m  # convert to "model millimeters"
    #plot_3d(enu_points)

    # TODO: downsample alg that keeps detail where needed adaptively
    enu_points = downsample_to(enu_points, 100)
    #plot_3d(enu_points)

    snapped_to_vertices = find_nearest_vertices(enu_points, vertices)
    snapped_to_vertices = filter_identical(snapped_to_vertices)
    #snapped_to_vertices = enu_points

    if hike_type == "LOOP":
        print("LOOP, putting first point at beginning")
        # Add a duplicate of the first point at the end, so that we get a closed
        # loop. Saves some steps in fusion360.
        print(f"First point: {snapped_to_vertices[0]}")
        snapped_to_vertices = np.vstack((snapped_to_vertices, snapped_to_vertices[0]))
        print(f"Last point: {snapped_to_vertices[-1]}")
    #plot_3d(snapped_to_vertices)


    return snapped_to_vertices


def downsample_to(points, n_desired):
    n = points.shape[0]
    take_every_nth = n // n_desired
    return points[::take_every_nth]


# In fusion 360, we have to:
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

def load_from_json(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data


def plot_3d(np_array):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot for individual points
    ax.scatter(np_array[:,0], np_array[:,1], np_array[:,2], c='r', marker='o', label='Points')

    # Line plot connecting the points
    ax.plot(np_array[:,0], np_array[:,1], np_array[:,2], color='blue', label='Line')

    # Add some labels for clarity
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.legend()
    plt.tight_layout()
    plt.show()

def highest_z_near_origin(points, radius=0.00001):
    """
    Find the index and coordinates of the point with the highest z-value
    among points whose (x,y) are within `radius` of (0,0).

    Parameters
    ----------
    points : np.ndarray
        Nx3 array of coordinates, where each row is [x, y, z].
    radius : float, optional
        Maximum distance from (0,0) in xy-plane to consider a point "close".
        Default is 0.1.

    Returns
    -------
    max_idx : int or None
        The index of the point with the highest z among those close to (0,0).
        Returns None if no points meet the criteria.
    max_coord : np.ndarray or None
        The coordinates [x, y, z] of the selected point.
        Returns None if no points meet the criteria.
    """
    # Compute distances from origin in xy-plane
    distances = np.sqrt(points[:,0]**2 + points[:,1]**2)

    # Filter points that are within the given radius
    close_mask = distances <= radius
    close_points = points[close_mask]

    if close_points.size == 0:
        # No points are within the specified radius
        return None, None

    # Find the index of the max z among the filtered points
    z_values = close_points[:, 2]
    local_max_idx = np.argmax(z_values)

    # Get the global index relative to the original array
    global_indices = np.where(close_mask)[0]
    max_idx = global_indices[local_max_idx]

    max_coord = points[max_idx]

    return max_idx, max_coord

# def enu_to_lla(point, origin_lon_lat_alt):
#     # Setup pyproj transformations
#     wgs84_geod = pyproj.Geod(ellps='WGS84')
#     lla = pyproj.Proj(proj='latlong', datum='WGS84')
#     ecef = pyproj.Proj(proj='geocent', datum='WGS84', ellps='WGS84')
#
#     origin_lon_deg = origin_lon_lat_alt[0]
#     origin_lat_deg = origin_lon_lat_alt[1]
#     origin_alt_m = origin_lon_lat_alt[2]
#
#     # Convert LLA origin to ECEF
#     x_origin_ecef, y_origin_ecef, z_origin_ecef = pyproj.transform(
#         lla, ecef, origin_lon_deg, origin_lat_deg, origin_alt_m, radians=False)
#
#     # Compute rotation matrix from ENU to ECEF
#     lat_rad = np.radians(origin_lat_deg)
#     lon_rad = np.radians(origin_lon_deg)
#
#     R = np.array([
#         [-np.sin(lon_rad),                  np.cos(lon_rad),                 0],
#         [-np.sin(lat_rad)*np.cos(lon_rad), -np.sin(lat_rad)*np.sin(lon_rad), np.cos(lat_rad)],
#         [ np.cos(lat_rad)*np.cos(lon_rad),  np.cos(lat_rad)*np.sin(lon_rad), np.sin(lat_rad)]
#     ])
#
#     # ENU to ECEF increment
#     dx, dy, dz = R.T.dot(np.array([point[0], point[1], point[2]]))
#
#     # Target ECEF coordinates
#     x_ecef = x_origin_ecef + dx
#     y_ecef = y_origin_ecef + dy
#     z_ecef = z_origin_ecef + dz
#
#     # Convert ECEF back to LLA
#     lon_target, lat_target, alt_target = pyproj.transform(ecef, lla, x_ecef, y_ecef, z_ecef, radians=False)
#     return np.array([lon_target, lat_target, alt_target])

def get_elevation_open_elevation(lat, lon):
    """
    Query the Open-Elevation API for the elevation of a given lat/lon.
    Returns elevation in meters, or None if the query fails.
    """
    url = "https://api.open-elevation.com/api/v1/lookup"
    params = {"locations": f"{lat}, {lon}"}
    try:
        response = requests.get(url, params=params)
        data = response.json()
        elevation = data["results"][0]["elevation"]  # in meters
        return elevation
    except Exception as e:
        print(f"Error retrieving elevation: {e}")
        return None

if __name__ == "__main__":

    ###### Inputs

    # ## South Arapaho
    # stl_path = r'D:\Projects\Python\GpxTracks3d\data\south_arapaho_2\10m_-105.64_40.02_tile_1_1.STL'
    # gpx_path = r'D:\Projects\Python\GpxTracks3d\data\south_arapaho_2\South_Arapaho.gpx'
    # csv_out_path = r'south_arapaho2.csv'
    # box_upper_right = (-105.60270900306992, 40.04201533473237, 0)
    # box_lower_left = (-105.6771528202086, 39.99179468640927, 0)
    # model_w_h_mm = (70, 61.4)
    # max_elev_m = 4117.2
    # out_n_back = True
    # ##

    # ## Loyalsock
    # stl_path = r'D:\Projects\Python\GpxTracks3d\data\loyalsock\10m_-76.55_41.48_tile_1_1.STL'
    # gpx_path = r'D:\Projects\Python\GpxTracks3d\data\loyalsock\Loyalsock_Link_Trail_Overnight.gpx'
    # csv_out_path = r'loyalsock.csv'
    # box_upper_right = (-76.50209028784886, 41.517332857084725, 0)
    # box_lower_left = (-76.60119612280978, 41.441583047037696, 0)
    # model_w_h_mm = (50, 50.8)
    # max_elev_m = 584.9
    # out_n_back = False
    # ##

    ## Boulder Three Peaks
    # stl_path = r'D:\Projects\Python\GpxTracks3d\data\boulder_three_peaks\10m_-105.29_39.97_tile_1_1.STL'
    # gpx_path = r'D:\Projects\Python\GpxTracks3d\data\boulder_three_peaks\boulder_three_peaks.gpx'
    # csv_out_path = r'boulder_three_peaks.csv'
    # box_upper_right = (-105.25157829589844, 40.01993951046981, 0)
    # box_lower_left = (-105.33517738647461, 39.92443275264491, 0)
    # model_w_h_mm = (100, 148.4)
    # max_elev_m = 2603
    ##

    ####### End Inputs
    json_inputs = load_from_json(os.path.join("data", "jess_boulder_5_peaks", "config.json"))
    snapped_points = snap_gpx_to_stl(json_inputs["gpx_path"],
                                     json_inputs["stl_path"],
                                     json_inputs["box_upper_right"],
                                     json_inputs["box_lower_left"],
                                     json_inputs["hike_type"])

    # Fusion 360 interprets 0.1 in the csv to mean 1mm. So we need to account
    # for this in the csv that we write for fusion360.
    scaled_snapped = snapped_points * 1.0/FUSION_360_MODEL_MM_PER_CSV
    write_csv(json_inputs["csv_out_path"], scaled_snapped)