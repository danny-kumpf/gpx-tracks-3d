import pymeshfix

in_file = r'D:\Projects\Python\GpxTracks3d\refined.stl'
out_file = r'D:\Projects\Python\GpxTracks3d\refined_fixed.stl'

# Read mesh from infile and output cleaned mesh to outfile
pymeshfix.clean_from_file(in_file, out_file)
