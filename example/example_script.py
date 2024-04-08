# main_script.py
import geopandas as gpd
import sys, os
sys.path.append(os.getcwd())
from WP3_stream_merge import fix_stream_network, split_lines_at_junctions


def main():
    # Set the working directory
    wd = os.path.dirname(os.getcwd())

    # Read in stream shapefile
    streams = gpd.read_file(fr'{wd}\Data\VP3_Fyn_streams.shp')
    streams_pp = split_lines_at_junctions(streams)
    
    # Read in terminus shapefile (Coastline)
    terminus = gpd.read_file(fr'{wd}\Data\coastline_buffer.shp')
    terminus = terminus.rename(columns={"Shape_Leng": "Length"})

    # Add buffer to the terminus geometry
    terminus['buffer'] = terminus['geometry'].buffer(50)
    terminus['geometry'] = terminus['buffer']
    terminus = terminus.drop(columns=['buffer'])

    # Set output folder path
    outpath = fr'{wd}\Output/'
    is_exist = os.path.exists(outpath)

    # Create output folder if it does not already exist
    if not is_exist:
        os.makedirs(outpath)
        print("The folder '" + str(outpath) + "' is created!")

    # Run upstream merge function
    fix_stream_network(stream_in = streams_pp, 
                       terminus  = terminus, 
                       outpath   = outpath, 
                       topo_id   = 'VP3')

if __name__ == "__main__":
    main()






