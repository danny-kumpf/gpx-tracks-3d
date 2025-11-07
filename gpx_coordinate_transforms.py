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
from plotting import plot_mesh_with_track, plot_mesh_with_track_w_sliders


def convert_gpx_track_to_stl_coordinates(gpx_file_path, stl_mesh, box_upper_right, box_lower_left):
    """
    Given a GPX file path and STL file path, this function returns a 3d
    track (Nx3 ndarray) created by converting the GPX LLA coordinates
    into the STL's coordinate frame.

    Parameters:
        gpx_file_path (str): Path to the input GPX file.
        stl_mesh (numpy-stl mesh.Mesh): stl mesh
        box_upper_right (): lon, lat, alt of upper-right model corner
        box_lower_left (): lon, lat, alt of lower left model corner
    """

    # Read the STL file
    vertices = get_stl_vertices(stl_mesh)

    # Create the LLA to ECEF transformer once
    transformer = create_transformer_lla_to_ecef()

    # Get the scale factor between real meters and "model mmm"
    mmm_per_m = get_model_scale(box_upper_right, box_lower_left, transformer, vertices)

    # Get the origin point (lon, lat, alt of lower left corner of model)
    origin = get_model_origin_lla(box_lower_left, mmm_per_m, vertices)

    # Extract all points from the GPX file
    gpx = read_gpx_file(gpx_file_path)
    lla_track_points = extract_points(gpx)

    summit_lat, summit_lon, summit_alt = highest_point_lat_lon(lla_track_points)
    print(f"Summit lat: {summit_lat}, Summit lon: {summit_lon}, Summit alt: {summit_alt}")

    # Convert points to ENU coordinates relative to the origin
    enu_points = convert_points_to_enu(lla_track_points, origin, transformer)
    enu_points = filter_min_distance(enu_points, 10)
    enu_points *= mmm_per_m  # convert to "model millimeters"
    #plot_3d(enu_points)
    plot_mesh_with_track_w_sliders(stl_mesh, enu_points)

    # move the points vertically until they sit as close to the surface as possible
    return move_closer_to_surface(enu_points, vertices)

def highest_point_lat_lon(lla_track_points: np.ndarray):
    """
    lla_track_points: Nx3 array [lat, lon, alt]
    Returns (lat, lon) at max altitude.
    """
    # alt is column 2
    idx = np.argmax(lla_track_points[:, 2])
    lon, lat, alt = lla_track_points[idx]
    return lat, lon, alt

def get_model_scale(box_upper_right, box_lower_left, transformer, stl_vertices):
    """
    Return model mm per actual real-life meters
    """
    # Only need these in ENU to measure distance; origin doesnt matter
    box_enu = convert_points_to_enu(np.array([box_upper_right, box_lower_left]),
                                    origin=box_lower_left,
                                    transformer=transformer)
    box_diag_m = distance(box_enu[0], box_enu[1])
    print(f"Box diagonal meters = {box_diag_m}")

    x_range_mmm = (stl_vertices[:, 0].min(), stl_vertices[:, 0].max())
    y_range_mmm = (stl_vertices[:, 1].min(), stl_vertices[:, 1].max())
    z_range_mmm = (stl_vertices[:, 2].min(), stl_vertices[:, 2].max())

    print(f"STL Number of vertices: {len(stl_vertices)}")
    print(f"STL X range: {x_range_mmm[0]} to {x_range_mmm[1]}")
    print(f"STL Y range: {y_range_mmm[0]} to {y_range_mmm[1]}")
    print(f"STL Z range: {z_range_mmm[0]} to {z_range_mmm[1]}")

    model_width_mmm = x_range_mmm[1] - x_range_mmm[0]
    model_height_mmm = y_range_mmm[1] - y_range_mmm[0]
    model_diag_mm = math.sqrt(model_width_mmm ** 2 + model_height_mmm ** 2)
    mmm_per_m = model_diag_mm / box_diag_m  # model mm per meter
    print(f"Model-mm-to-Meters scale factor: {1.0 / mmm_per_m}")
    mmm_per_m = 1.0 / 220.45111
    print(f"HARDCODING Model-mm-to-Meters scale factor: {1.0 / mmm_per_m}")
    # TODO: suspicious that this doesn't agree with the logfile. Not a big
    #       enough diff to cause the shifts I'm seeing, apparently.
    return mmm_per_m


def get_model_origin_lla(box_lower_left, mmm_per_m, stl_vertices):
    """
    The [0,0,0] point on the model is the lower-left corner. Find the
    lat-long-alt coordinates of that point in real-space.
    """
    # this elevation isn't super accurate, but we correct for it later by
    # iteratively snapping to STL vertices
    #origin_elev = get_elevation_open_elevation(box_lower_left[1], box_lower_left[0])
    origin_elev = get_elevation_opentopodata(box_lower_left[1], box_lower_left[0])

    print(f"Elevation of bottom-left corner ground surface: ~{origin_elev}")
    _, origin_top_pt = highest_z_near_origin(stl_vertices)
    origin_elev -= origin_top_pt[2] / mmm_per_m
    origin = (box_lower_left[0],
              box_lower_left[1],
              origin_elev)
    print(f"Elevation of bottom-left model corner (m): ~{origin[2]}")
    return origin


def get_stl_vertices(stl_mesh):
    """
    Reads an STL file and extracts unique vertices.

    Parameters:
        stl_file_path (str): Path to the STL file.

    Returns:
        vertices (numpy.ndarray): Array of unique vertices (Nx3).
    """
    # Concatenate all vertices from the triangles
    all_vertices = np.concatenate((stl_mesh.v0, stl_mesh.v1, stl_mesh.v2))
    # Get unique vertices
    vertices = np.unique(all_vertices, axis=0)
    return vertices


def move_closer_to_surface(points: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """
    "Snap" each track point to the surface in the Z direction.

    For each point (x,y,z), find the mesh vertex with nearest (X,Y)
    and set z := Z of that vertex.

    Parameters:
        points   : (M,3) float array
        vertices : (N,3) float array

    Returns:
        adjusted_points : (M,3) float array
    """
    tree_xy = KDTree(vertices[:, :2])               # build on XY only
    _, idx = tree_xy.query(points[:, :2])           # nearest in XY
    adjusted = points.copy()
    adjusted[:, 2] = vertices[idx, 2]               # z to surface z
    return adjusted

# def shift_track_to_surface(track_xyz: np.ndarray, stl_vertices: np.ndarray) -> np.ndarray:
#     """
#     track_xyz: (N,3) array of track points in same coordinate space as STL
#     stl_vertices: (M,3) array of vertex positions from the STL mesh
#
#     Returns shifted track (N,3) with one global offset applied.
#     """
#     # kd-tree on STL verts
#     tree = KDTree(stl_vertices)
#
#     # nearest vertex index for each track point
#     _, idx = tree.query(track_xyz)
#
#     # nearest vertex positions
#     nearest = stl_vertices[idx]
#
#     # vectors from vertex -> track point
#     vectors = track_xyz - nearest
#
#     # average vector
#     avg_vec = vectors.mean(axis=0)
#     print(f"Average offset track to surface: {avg_vec}")
#
#     # shift entire track
#     return track_xyz + avg_vec



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
                lon_deg = point.longitude # + 0.0015  # Degrees
                lat_deg = point.latitude   # Degrees
                altitude_m = point.elevation if point.elevation is not None else 0.0  # Meters
                points.append((lon_deg, lat_deg, altitude_m))# - 17.0))

    return np.array(points)


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


def distance(point_1, point_2):
    return np.linalg.norm(point_1 - point_2)


def downsample_to(points, n_desired):
    n = points.shape[0]
    take_every_nth = n // n_desired
    return points[::take_every_nth]


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
    print(f"Number of close points to origin: {len(close_points)}")
    print(f"They are: ")
    for point in close_points:
        print(f"  {point}")

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
    print(f"params: {params}")
    try:
        response = requests.get(url, params=params)
        print(f"response: {response}")
        data = response.json()
        print(f"data: {data}")
        elevation = data["results"][0]["elevation"]  # in meters
        return elevation
    except Exception as e:
        print(f"Error retrieving open_elevation elevation: {e}")
        raise e

def get_elevation_opentopodata(lat, lon, dataset="srtm90m"):
    """
    Query the OpenTopoData API for elevation in meters.
    Returns elevation or None on failure.
    """
    url = f"https://api.opentopodata.org/v1/{dataset}"
    params = {"locations": f"{lat},{lon}"}

    try:
        r = requests.get(url, params=params, timeout=10)
        data = r.json()
        return data["results"][0]["elevation"]
    except Exception as e:
        print(f"Error retrieving opentopodata elevation: {e}")
        raise e



