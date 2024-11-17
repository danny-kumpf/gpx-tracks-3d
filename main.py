import gpxpy
import csv
import math
from pyproj import Transformer
import numpy as np
from stl import mesh
from scipy.spatial import KDTree

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

def get_origin(gpx):
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

def create_transformer():
    """
    Creates and returns the LLA (Latitude, Longitude, Altitude) to ECEF (Earth-Centered, Earth-Fixed) transformer.

    Returns:
        transformer (pyproj.Transformer): Transformer object for LLA to ECEF conversion.
    """
    # Initialize the transformer once for efficiency
    transformer = Transformer.from_crs("epsg:4979", "epsg:4978", always_xy=True)
    return transformer

def precompute_origin_parameters(lat_deg, lon_deg):
    """
    Precomputes sine and cosine of the origin's latitude and longitude in radians.

    Parameters:
        lat_deg (float): Latitude of the origin in degrees.
        lon_deg (float): Longitude of the origin in degrees.

    Returns:
        sin_latitude (float): Sine of the origin's latitude in radians.
        cos_latitude (float): Cosine of the origin's latitude in radians.
        sin_longitude (float): Sine of the origin's longitude in radians.
        cos_longitude (float): Cosine of the origin's longitude in radians.
    """
    # Convert latitude and longitude to radians
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    # Precompute sines and cosines
    sin_latitude = math.sin(lat_rad)
    cos_latitude = math.cos(lat_rad)
    sin_longitude = math.sin(lon_rad)
    cos_longitude = math.cos(lon_rad)
    return sin_latitude, cos_latitude, sin_longitude, cos_longitude

def extract_points(gpx):
    """
    Extracts all points from the GPX tracks.

    Parameters:
        gpx (gpxpy.gpx.GPX): Parsed GPX object.

    Returns:
        points (list): List of tuples containing (lon_deg, lat_deg, altitude_m).
    """
    points = []
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                lon_deg = point.longitude  # Degrees
                lat_deg = point.latitude    # Degrees
                altitude_m = point.elevation if point.elevation is not None else 0.0  # Meters
                points.append((lon_deg, lat_deg, altitude_m))
    return points

def convert_points_to_enu(points, origin, transformer):
    """
    Converts a list of geodetic coordinates to ENU coordinates relative to the origin.

    Parameters:
        points (list): List of tuples containing (lon_deg, lat_deg, altitude_m).
        origin (tuple): Tuple containing (lon_deg, lat_deg, altitude_m).
        transformer (pyproj.Transformer): Transformer for LLA to ECEF conversion.

    Returns:
        enu_points (list): List of tuples containing (east_m, north_m, up_m).
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
    enu_points = []

    for lon_deg, lat_deg, altitude_m in points:
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

        enu_points.append((east_m, north_m, up_m))
    return enu_points

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
    dx = point_1[0] - point_2[0]
    dy = point_1[1] - point_2[1]
    dz = point_1[2] - point_2[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


FUSION_360_MODEL_MM_PER_CSV = 10


def snap_gpx_to_stl(gpx_file_path, stl_path, box_upper_right, box_lower_left, model_w_h_mm):
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
    z_max = vertices[:, 2].max()
    print(f"STL Number of vertices: {len(vertices)}")
    print(f"STL X range: {vertices[:, 0].min()} to {vertices[:, 0].max()}")
    print(f"STL Y range: {vertices[:, 1].min()} to {vertices[:, 1].max()}")
    print(f"STL Z range: {vertices[:, 2].min()} to {vertices[:, 2].max()}")

    # Get the origin point (lon, lat, alt of lower left corner of model)
    origin = (box_lower_left[0], box_lower_left[1], 1479)

    # Create the LLA to ECEF transformer once
    transformer = create_transformer()

    # Extract all points from the GPX file
    points = extract_points(gpx)

    # Convert points to ENU coordinates relative to the origin
    enu_points = convert_points_to_enu(points, origin, transformer)

    # Fusion 360 interprets csv valuesu as "0.1 in csv = 1mm in fusion360"
    # 8000m real -> 100 mm on model (approx.)
    # 100 mm in fusion360 = 10.0 in the csv
    # so 8000.0 in the csv (real) -> 10.0 in the csv
    # scale factor is divide by 800
    box_enu = convert_points_to_enu([box_upper_right, box_lower_left], origin, transformer)
    box_diag_m = distance(box_enu[0], box_enu[1])
    print(f"Box diagonal meters = {box_diag_m}")
    model_diag_mm = math.sqrt(model_w_h_mm[0] ** 2 + model_w_h_mm[1] ** 2)
    # print(f"Model diagonal mm = {model_diag_mm}")
    # csv_diag = model_diag_mm / FUSION_360_MODEL_MM_PER_CSV
    # print(f"CSV diagonal needs to be: {csv_diag}")
    # scale_factor = csv_diag / box_diag_m
    scale_factor = model_diag_mm / box_diag_m
    print(f"Meters-to-model-mm scale factor: {scale_factor}")
    enu_points = scale(enu_points, scale_factor)

    # TODO: downsample alg that keeps detail where needed adaptively
    enu_points = downsample_to(enu_points, 100)

    enu_points_np = np.ndarray((len(enu_points), 3))
    for idx in range(len(enu_points)):
        enu_points_np[idx][0] = enu_points[idx][0]
        enu_points_np[idx][1] = enu_points[idx][1]
        enu_points_np[idx][2] = enu_points[idx][2]

    closest_vertices = find_nearest_vertices(enu_points_np, vertices)
    closest_vertices = scale(closest_vertices, 1.0/FUSION_360_MODEL_MM_PER_CSV)

    return closest_vertices


def downsample_to(points, n_desired):
    n = len(points)
    take_every_nth = n // n_desired
    new_points = []
    for i in range(n):
        if i % take_every_nth == 0:
            new_points.append(points[i])
    return new_points


def scale(enu_points, scale_factor):
    """
    0.1 in the csv = 1mm in Fusion360
    :return:
    """
    # 8000m real -> 100 mm on model
    # 100 mm in fusion360 = 10.0 in the csv
    # so 8000.0 in the csv (real) -> 10.0 in the csv
    # scale factor is divide by 800
    enu_scaled = []
    for east_m, north_m, up_m in enu_points:
        enu_scaled.append((east_m*scale_factor,
                           north_m*scale_factor,
                           up_m*scale_factor))
    return enu_scaled


# Example usage:
# Replace 'input.gpx' with your GPX file path
# Replace 'output.csv' with your desired CSV output path
# gpx_to_enu_csv('input.gpx', 'output.csv')


# In fusion 360, we have to:
# - load the STL (import mesh)
# - load the snapped track that this script creates (utilities, add-ins, scripts and add-ins, ImportSplineCSV)
#   (this script aligns the coordinate systems; no manual movement in f360 needed)
# - create a plane perpendicular to the track (solid, construct, plane along path)
# - in that plane, create the cut cross-section centered on the path
# - create a swept SURFACE (cant do solid bc path intersects itself)
# - create end caps: surface, patch, click the end caps of the new surface pipe
# - fill the surface: solid, create, boundary fill
# - convert the surface to mesh: mesh, tesselate (can do low quality)
# - merge the two meshes as a cut: mesh, modify, combine (has to cook for a few min)


if __name__ == "__main__":

    ## Inputs

    stl_path = r'D:\Projects\Python\GpxTracks3d\data\south_arapaho\10m_-105.65_40.01_tile_1_1.STL'
    gpx_path = r'D:\Projects\Python\GpxTracks3d\data\south_arapaho\South_Arapaho.gpx'
    csv_out_path = r'south_arapaho.csv'

    # Get these values from expanding the "Area Selection Box" field of touchterrain
    # data/boulder_three_peaks
    # box_upper_right = (-105.25157829589844, 40.01993951046981, 0)
    # box_lower_left = (-105.33517738647461, 39.92443275264491, 0)
    # model_w_h_mm = (100, 148.4)

    # data/south_arapaho
    box_upper_right = (-105.6259717590332, 40.03538452772285, 0)
    box_lower_left = (-105.67412277526856, 39.99094728785393, 0)
    model_w_h_mm = (30, 36)

    ## End Inputs

    snapped_points = snap_gpx_to_stl(gpx_path,
                                     stl_path,
                                     box_upper_right,
                                     box_lower_left,
                                     model_w_h_mm)

    snapped_points_csv = []
    for row in snapped_points:
        snapped_points_csv.append((row[0], row[1], row[2]))
    write_csv(csv_out_path, snapped_points_csv)