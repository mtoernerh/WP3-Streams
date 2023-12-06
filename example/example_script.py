# main_script.py
import geopandas as gpd
import os
from code.WP3_stream_merge import fix_stream_network

def main():
    # Set the working directory
    working_directory = os.environ.get("WORKING_DIRECTORY", r"M:\VP3_Fyn\Renew\Python")
    os.chdir(working_directory)

    # Read in stream shapefile
    streams = gpd.read_file('data/VP3_Fyn_streams.shp')
    streams = streams.rename(columns={"Shape_Leng": "Length"})

    # Read in terminus shapefile (Coastline)
    terminus = gpd.read_file('data/VP3_Fyn_coastline.shp')
    terminus = terminus.rename(columns={"Shape_Leng": "Length"})

    # Add buffer to the terminus geometry
    terminus['buffer'] = terminus['geometry'].buffer(100)
    terminus['geometry'] = terminus['buffer']
    terminus = terminus.drop(columns=['buffer'])

    # Set output folder path
    outpath = 'output/'
    is_exist = os.path.exists(outpath)

    # Create output folder if it does not already exist
    if not is_exist:
        os.makedirs(outpath)
        print("The folder '" + str(outpath) + "' is created!")

    # Run upstream merge function
    fix_stream_network(streams, terminus, outpath, topo_id='VP3')

if __name__ == "__main__":
    main()
