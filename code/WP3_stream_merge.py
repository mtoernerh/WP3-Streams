# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 13:05:08 2022

@author: mfth
"""
# =============================================================================#
# Upstream_merge function #
# =============================================================================#
import geopandas
import pandas as pd
from geopandas.tools import overlay
from shapely.ops import linemerge
import shapely.ops
from shapely.wkt import loads
from shapely.geometry import LineString, Point, MultiLineString
from shapely.ops import split
import warnings
import numpy as np
from shapely.ops import polygonize
import shapely.geometry as geom
import os
import sys
import glob
from tqdm import tqdm

warnings.filterwarnings("ignore")


def reverse_geom(geom):
    def _reverse(x, y, z=None):
        if z:
            return x[::-1], y[::-1], z[::-1]
        return x[::-1], y[::-1]

    return shapely.ops.transform(_reverse, geom)


def merge_streams_upstream(stream_input, name_string, terminus):
    """
    Parameters
    ----------
    stream_input : Polyline shapefile (.shp)
        A polyline feature class (.shp) that must have a Name/ID (string) column, a Length
        (float/int) column and lines should be topologically connected and diveded in different
        features at each intersection.
    name_string : String (Str)
        A string that can be anything, although it is appended to fields in the Name/ID columm.
        An example can be MS 'Mainstem', SS1 'substem 1' or SS2 'substem 2' which represents
        different stream orders.
    terminus : shapefile (.shp)
        At least one shapefile (.shp) that intersects with furthest downstream segment of Stream_input
        must be defined. That can be already defined streams that Stream_input feeds into, a coastline
        or lakes. This feature is used to identify the furtherst downstream part of Stream_input.
        It is important that the stream_input is topologically touching the terminus feature otherwise
        the function will not work properly. The shapefile can be points, lines or polygons or a
        combination although polygons have been the most succusful.


    Returns : GeoDataFrame (gdf)
        A GeoDataFrame that contains processed Stream_input features.
        This function identifies the fursterst downstream segment of a stream input
        that intersects with specified terminus feature(s) and merges with the next upstream
        segment, if there is one, based on total upstream length. If this function
        is run in a loop, where the symmetric difference between the output, which are
        single merged lines intersecting with e.g. already established streams, and the input:
        any sub-tributaries that are feeding into the output streams can be isolated and used as input
        for the next iteration, where the output from last iteration is then used as terminus.
        This way every sub-tributary can be defined and exported in its own shapefile.
        See the following PowerPoint for a visual example of what the function does.
        \\netapp2p\Dkmodel-hydro\projekter\Vandlb_VPIII\python\Merge_upstream\Merge_upstream_function.pptx
    ----------

    """
    stream_input_crs = (
        stream_input.crs
    )  # Extracts CRS/SRS string to respecify later for output
    stream_input_srs = stream_input_crs.srs
    Stream_out = pd.DataFrame(columns=stream_input.columns).astype(
        stream_input.dtypes
    )  # Output GeoDataFrame is defined as an empty copy of the input while retaining the same data format
    Temp_1 = pd.DataFrame(columns=stream_input.columns).astype(
        stream_input.dtypes
    )  # Temporary GeoDataFrame for merging stream segments
    Temp_2 = pd.DataFrame(columns=stream_input.columns).astype(
        stream_input.dtypes
    )  # Temporary GeoDataFrame for merging stream segments

    stream_input = stream_input.reset_index(
        drop=True
    )  # Index order is reset for stream input, only matters if function is used in a loop that might alter index
    object_column_name = (stream_input.select_dtypes(include=[object])).columns[
        0
    ]  # Isolating the first input columnn that contain string data

    counter = 0  # Iterable counter used for various indexing
    # =============================================================================#
    # Iterating over each row (line feature) in Stream input  #
    # =============================================================================#
    for index, row in tqdm(
        stream_input.iterrows(),
        total=len(stream_input),
        desc="Merging " + str(name_string[1 : len(name_string)]),
        colour="#8dd3c7",
    ):
        Stream_table = stream_input.copy(
            deep=True
        )  # Make a copy of the input data that resets after each iteration
        current_stream = stream_input.iloc[
            [counter]
        ]  # Select the current stream feature as an isolated GeoDataFrame
        Stream_i_index = (
            counter  # Index of the current stream feature before the counter increases
        )
        Stream_table.loc[
            index, "geometry"
        ] = None  # Set the geometry of the current stream feature to None in the copy of the input data
        # =============================================================================#
        # Intersection between Stream feature i and termini features #
        # =============================================================================#
        terminus_intersect = overlay(
            current_stream, terminus, how="intersection", keep_geom_type=False)  # Intersection between the current stream feature and terminus features
        counter += 1  # Increment the counter
        touch_check = 1  # Set a flag indicating that the current stream feature touches one of the terminus features
        # =============================================================================#
        # Merging between stream i and upstream segments #
        # =============================================================================#
        if (
            len(terminus_intersect) > 0
        ):  # Loop starts if stream i intersects a termini feature
            while (
                touch_check > 0
            ):  # While loop is active as long as stream i touches a new upstream feature
                upstream_features = Stream_table.loc[
                    Stream_table.touches(current_stream.loc[index, "geometry"])
                ]  # Returns upstream features touching stream i as a GeoDataFrame
                upstream_features_temp = upstream_features.copy(
                    deep=True
                )  # Creating a copy of touching upstream features, which is used a temporary geoproccessing step
                touch_feature = 0
                for (
                    index2,
                    row,
                ) in (
                    upstream_features_temp.iterrows()
                ):  # Iterating over number of touching upstream segments
                    upstream_terminus_intersect = overlay(
                        upstream_features_temp.iloc[[touch_feature]],
                        terminus,
                        how="intersection",
                        keep_geom_type=False,
                    )  # Checking if upstream touching segment is also touching a termini feature
                    touch_feature = touch_feature + 1
                    if (
                        len(upstream_terminus_intersect) > 0
                    ):  # If upstream touching segment is also touching terminus feature then...
                        upstream_features = upstream_features.drop(
                            labels=index2, axis=0
                        )  # Remove that feature from GeoDataFrame containing touching upstream streams
                        Stream_table.loc[
                            index2, "geometry"
                        ] = None  # Removing the given feature from temporary stream input, so it cannot be encountered again during this iteration

                upstream_indexes = list(upstream_features.index.values)
                upstream_geometries = upstream_features[
                    "geometry"
                ]  # Returns GeoDataFrame only containing geometries for upstream features
                upstream_lengths = upstream_features[
                    "Length"
                ]  # Creating a list of lengths of upstream touching streams, only used to check if there are any upstream segments
                # =============================================================================#
                # Two different scenarios, either there are upstream features or not #
                # =============================================================================#
                if (len(upstream_lengths) == 0):  # If Length is 0, then Stream i does not touch another upstream segment
                    stream_input.loc[index, object_column_name] = (stream_input.loc[index, object_column_name] + name_string)  # Adding pre-defined name string to ID name for current_stream
                    Stream_out = pd.concat([Stream_out, stream_input.iloc[[index]]])
                    touch_check = len(
                        Stream_table.loc[
                            Stream_table.touches(current_stream.loc[index, "geometry"])
                        ]
                    )  # Touch check to see if there are any touching upstream segments (only NEW touching), if 0 the while loop ends here
                    Stream_table = stream_input.copy(deep=True)  # Maybe not neccessary
                # =============================================================================#
                # The total upstream lengths are compared between touching upstream segments #
                # =============================================================================#
                else:  # When there is at least 1 upstream touching feature, this loop starts
                    total_upstream_length = upstream_features["Length"].tolist()  # Returns a list of the length of each touching upstream segment
                    touch_feature = 0  # Counter to iterate over each tóuching upstream segment
                    for touch_feature in range(len(upstream_features)):
                        Stream_length_table = Stream_table.copy(deep=True)  # A temporary copy of already temp Stream_table used in this loop only
                        
                        upstream_indexes = list(upstream_features.index.values)  # List of index value for touching upstream segments
                        
                        Stream_length_table.loc[upstream_indexes, "geometry"] = None  # Touching upstream segments geometries are set to None so they cannot be encountered again
                       
                        #current_upstream_feature = upstream_features.iloc[[touch_feature]]  # Selecting ith upstream touching feature
                        #current_upstream_feature = (current_upstream_feature.unary_union)  # Might not be neccesary
                        current_upstream_feature = upstream_features.iloc[[touch_feature]]["geometry"].unary_union
                        while_check = 1  # Another abitrary while check, that can use any value above 0
                        
                        while (while_check > 0):  # A while loop is active as long as there is still another new upstream touching feature
                            Stream_upstream_touch = Stream_length_table.loc[Stream_length_table.touches(current_upstream_feature)]  # The next touching upstream feature is found
                            
                            #current_upstream_feature = shapely.ops.unary_union(Stream_upstream_touch["geometry"])  # Neccessary in case multiple upstream features, because these operation has to work on single features
                           
                            upstream_indexes = list(Stream_upstream_touch.index.values)  # Removing the already encountered features from future encounters
                            
                            Stream_length_table.loc[upstream_indexes, "geometry"] = None
                           
                            total_upstream_length[touch_feature] = (total_upstream_length[touch_feature]+ Stream_upstream_touch["Length"].sum())  # Adding the upstream length to the first encountered upstream feature
                            
                            geometry_upstream = Stream_upstream_touch["geometry"]
                            
                            upstream_geometries = pd.concat([upstream_geometries, geometry_upstream])
                            
                            current_upstream_feature = (geometry_upstream).unary_union
                            
                            if len(Stream_upstream_touch) == 0:
                                while_check = 0
                            else:
                                while_check = len(Stream_length_table.loc[Stream_length_table.touches(current_upstream_feature)])
                            
                        #touch_feature = touch_feature + 1

                    geometry_touch_list = upstream_geometries.tolist()
                    geometry_touch_list = shapely.ops.unary_union(geometry_touch_list)
                    upstream_geomtries = shapely.ops.polygonize_full(
                        geometry_touch_list
                    )
                    loop_poly = geom.MultiPolygon(upstream_geomtries[0])
                    loop_poly_buffer = loop_poly.buffer(0.1)

                    upstream_max_length_id = np.argmax(total_upstream_length)
                    Stream_upstream_merge = upstream_features.iloc[
                        [upstream_max_length_id]
                    ]
                    Stream_id_merge = (Stream_upstream_merge.index)[0]

                    geom_string = str(
                        Stream_upstream_merge.loc[Stream_id_merge, "geometry"]
                    )
                    geom_string_wkt = shapely.wkt.loads(geom_string)

                    within_loop = geom_string_wkt.within(loop_poly_buffer)

                    if within_loop == True:
                        upstream_max_length_id = np.argmin(
                            upstream_features["Length"]
                        )  # CHECK IF WORKS #
                        Stream_upstream_merge = upstream_features.iloc[
                            [upstream_max_length_id]
                        ]
                        Stream_id_merge = (Stream_upstream_merge.index)[0]

                    Temp_1 = pd.concat([Temp_1, current_stream.iloc[[0]]])
                    Temp_2 = pd.concat([Temp_2, stream_input.iloc[[Stream_id_merge]]])
                    Temp_1[
                        "dissolvefield"
                    ] = 1  # Create abitrary common field, neccessary for merging
                    Temp_2[
                        "dissolvefield"
                    ] = 1  # Create abitrary common field, neccessary for merging
                    Temp_join = geopandas.GeoDataFrame(
                        pd.concat([Temp_1, Temp_2])
                    )  # The two temporary GdFs are combined
                    Stream_merge = Temp_join.dissolve(
                        "dissolvefield"
                    )  # The two features are merged into one
                    Stream_merge = pd.DataFrame(
                        columns=stream_input.columns, data=Stream_merge
                    )  # Reordering of columns to match current_stream
                    Stream_merge.index = [
                        index
                    ]  # Redifining index to match current_stream index, in case it changes during merging
                    current_stream = Stream_merge.copy(
                        deep=True
                    )  # current_stream is replaced with new merged feature
                    Temp_1 = Temp_1[0:0]  # Reset to an empty gdf
                    Temp_2 = Temp_2[0:0]  # Reset to an empty gdf
                    # =============================================================================#
                    # After merging current_stream with an upstream segment, any encountered upstream
                    # segments are removed from Stream_table and the new merged feature is checked
                    # to see if new upstream segments are encountered. In which case, the while loop
                    # repeats with above steps except current_stream is replaced with the merged feature.
                    # =============================================================================#
                    Stream_touching_index = (
                        upstream_features.index.tolist()
                    )  # List of indicies for any touching upstream segment
                    Stream_table.loc[
                        Stream_touching_index, "geometry"
                    ] = None  # Set Geometry of touching upstream segments to None, so they are not encountered in next touch_check
                    stream_input.loc[Stream_id_merge, "geometry"] = None
                    upstream_lengths = [1]
                    touch_check = len(
                        Stream_table.loc[
                            Stream_table.touches(
                                current_stream.loc[Stream_i_index, "geometry"]
                            )
                        ]
                    )  # Checks if merged feature touches any new upstream segments, if not the while loop ends, otherwise above steps are repeated for the merged feature
            Stream_table = stream_input.copy(deep=True)
            if (
                len(terminus_intersect) > 0 and len(upstream_lengths) != 0
            ):  # When while loop ends, run this part if a feature was merged
                Stream_merge.loc[index, object_column_name] = (
                    Stream_merge.loc[index, object_column_name] + name_string
                )  # Add name string to name ID
                Stream_out = pd.concat([Stream_out, Stream_merge.iloc[[0]]])
                touch_check = len(
                    Stream_table.loc[
                        Stream_table.touches(
                            current_stream.loc[Stream_i_index, "geometry"]
                        )
                    ]
                )  # Add merged feature to Stream_out
                Stream_table = stream_input.copy(deep=True)
    Stream_out = geopandas.GeoDataFrame(
        Stream_out, crs=stream_input_srs, geometry="geometry"
    )
    Stream_out = Stream_out.reset_index(drop=True)
    return Stream_out


"""
 Example of using the 'Merge_upstream' function iterately, with the output of
 one iteration as input for the next iteration. The result are several .shp files
 that each contain merged stream features, and together they form the complete
 branching stream network of the initial stream input.
"""


def fix_stream_network(stream_in, terminus, outpath, topo_id):
    """


    Parameters
    ----------
    stream_in : GeoDataBase
        A polyline feature class (.shp) that must have a Name/ID (string) column, a Length
        (float/int) column and lines should be topologically connected and diveded in different
        features at each intersection. Errors arise if these requirements are not met.
    terminus : GeoDataBase
        A shapefile (.shp) that intersects with furthest downstream segment of stream_in
        must be defined. That can be already defined streams that stream_in feeds into, a coastline
        or lakes. These features are used to identify the furtherst downstream part of stream_in.
        It is important that the stream_in is topologically touching the terminus features, otherwise
        the function will not work properly. The shapefile can be points, lines or polygons or a
        combination although polygons have been the most succusful..
    outpath : String
        The name of output path e.g., .../output/, in which output files are saved.
    topo_id : String
        An ID string that identifies the dataset such as VP3 (Vandplan 3) that is added to the filename,
        in the TopoID column and to the feature names in the object ID column.

    Returns
    -------
    GeoDatabase/Shp files.

    """
    object_column_name = (stream_in.select_dtypes(include=[object])).columns[0]
    string_list = ["_MS", "_SS1", "_SS2", "_SS3", "_SS4", "_SS5", "_SS6", "SS7", "SS8"]
    merged_streams_list = []

    i = 0
    for n in string_list:
        if len(stream_in) == 0:
            break
        try:
            Stream_out = merge_streams_upstream(
                stream_input=stream_in, name_string=n, terminus=terminus
            )
            print("\n Postprocessing " + str(n[1 : len(n)]) + "!\n")
            # =============================================================================#
            # Creating a copy of Stream_out, which will be edited for exporting to .shp #
            # Stream_out should not be edited as it must retain the same column structure
            # for the next iteration.
            # =============================================================================#
            Stream_out_PP = Stream_out.copy(deep=True)
            Stream_out_PP = Stream_out_PP.reset_index(drop=True)
            Stream_out_PP["BR_StartCh"] = 0
            Stream_out_PP["BR_EndCh"] = Stream_out_PP.geometry.length
            Stream_out_PP["BR_EndCh"] = (Stream_out_PP["BR_EndCh"] + 0.5).astype(int)
            Stream_out_PP["BR_FlowDir"] = 0
            Stream_out_PP["BR_Type"] = 2
            Stream_out_PP["LeakCoef"] = None
            Stream_out_PP["TopoID"] = topo_id
            Stream_out_PP["StreamOrd"] = str(i + 1)

            # Stop if there are no features left to process
            if len(Stream_out_PP.geometry) == 0:
                break

            for index, row in tqdm(
                Stream_out_PP.iterrows(),
                total=len(Stream_out_PP),
                desc="Linemerging " + str(n[1 : len(n)]),
                colour="#feb24c",
            ):
                Stream_out_PP.loc[index, object_column_name] = (
                    Stream_out_PP.loc[index, object_column_name]
                    + "_"
                    + str(Stream_out_PP.loc[index, "BR_EndCh"])
                )
                geom_string = str(
                    Stream_out_PP.loc[index, "geometry"]
                )  # selects geometry string of ith feature
                geom_string_wkt = shapely.wkt.loads(geom_string)

                try:
                    Stream_out_PP.loc[index, "geometry"] = linemerge(geom_string_wkt)
                except ValueError:
                    continue

            for index, row in tqdm(
                Stream_out_PP.iterrows(),
                total=len(Stream_out_PP),
                desc="Direction reversal " + str(n[1 : len(n)]),
                colour="#feb24c",
            ):
                Stream_table_i = Stream_out_PP.copy(deep=True)
                current_stream = Stream_out_PP.iloc[[index]]
                Stream_table_i.loc[index, "geometry"] = None

                Start_vertex = Point((Stream_out_PP.loc[index, "geometry"]).coords[0])
                Start_vertex = geopandas.GeoDataFrame(
                    index=[0], crs="epsg:32632", geometry=[Start_vertex]
                )

                End_vertex = Point((Stream_out_PP.loc[index, "geometry"]).coords[-1])
                End_vertex = geopandas.GeoDataFrame(
                    index=[0], crs="epsg:32632", geometry=[End_vertex]
                )

                Terminus_isect = overlay(
                    Start_vertex, terminus, how="intersection", keep_geom_type=False
                )  # Intersection between start vertex and MHYDRO
                Stream_isect = overlay(
                    Start_vertex,
                    Stream_table_i,
                    how="intersection",
                    keep_geom_type=False,
                )  # Intersection between start vertex and line features in VP3 copy
                End_isect = overlay(
                    Start_vertex, End_vertex, how="intersection", keep_geom_type=False
                )  # Intersection between start vertex and end vertex
                Stream_table_i = Stream_out_PP.copy(deep=True)
                if (
                    len(Terminus_isect) + len(Stream_isect) + len(End_isect)
                ) > 0:  # If VP3 start vertex intersects any of above features, then reverse order of coordinates as below
                    geom_string = str(
                        current_stream.loc[index, "geometry"]
                    )  # selects geometry string of ith feature
                    geom_string_wkt = shapely.wkt.loads(
                        geom_string
                    )  # Converts to wkt Shapely file format
                    geom_string_reverse = str(
                        reverse_geom(geom_string_wkt)
                    )  # Applies geometry coordinate reversal function
                    Stream_out_PP.loc[index, "geometry"] = loads(
                        geom_string_reverse
                    )  # Replaces old geom with new reversed geom in a VP3 output table
            if n != string_list[0]:
                list_of_files = glob.glob(outpath + "*.shp")
                latest_file = max(list_of_files, key=os.path.getmtime)
                test_int_dir = geopandas.GeoDataFrame.from_file(latest_file)
                for index, row in Stream_out_PP.iterrows():
                    Stream_table_i = Stream_out_PP.copy(deep=True)
                    current_stream = Stream_out_PP.iloc[[index]]
                    Stream_table_i.loc[index, "geometry"] = None

                    Start_vertex = Point(
                        (Stream_out_PP.loc[index, "geometry"]).coords[0]
                    )
                    Start_vertex = geopandas.GeoDataFrame(
                        index=[0], crs="epsg:32632", geometry=[Start_vertex]
                    )
                    End_vertex = Point(
                        (Stream_out_PP.loc[index, "geometry"]).coords[-1]
                    )
                    End_vertex = geopandas.GeoDataFrame(
                        index=[0], crs="epsg:32632", geometry=[End_vertex]
                    )

                    Terminus_isect = overlay(
                        Start_vertex,
                        test_int_dir,
                        how="intersection",
                        keep_geom_type=True,
                    )  # Intersection between start vertex and MHYDRO
                    if (
                        len(Terminus_isect)
                    ) > 0:  # If VP3 start vertex intersects any of above features, then reverse order of coordinates as below
                        test = test_int_dir.loc[
                            test_int_dir[object_column_name]
                            == Terminus_isect.loc[0, object_column_name]
                        ]
                        test_id = (test.index)[0]
                        Terminus_end_vertex = Point(
                            (test.loc[test_id, "geometry"]).coords[-1]
                        )
                        Terminus_end_vertex = geopandas.GeoDataFrame(
                            index=[0], crs="epsg:32632", geometry=[Terminus_end_vertex]
                        )

                        End_end_distance = int(Terminus_end_vertex.distance(End_vertex))
                        End_start_distance = int(
                            Terminus_end_vertex.distance(Start_vertex)
                        )
                        if End_end_distance > End_start_distance:
                            geom_string = str(
                                current_stream.loc[index, "geometry"]
                            )  # selects geometry string of ith feature
                            geom_string_wkt = shapely.wkt.loads(
                                geom_string
                            )  # Converts to wkt Shapely file format
                            geom_string_reverse = str(
                                reverse_geom(geom_string_wkt)
                            )  # Applies geometry coordinate reversal function
                            Stream_out_PP.loc[index, "geometry"] = loads(
                                geom_string_reverse
                            )

            for index, row in Stream_out_PP.iterrows():
                Stream_out_PP.loc[index, object_column_name] = (
                    Stream_out_PP.loc[index, object_column_name] + "_" + str(index)
                )

            Stream_out_save = Stream_out_PP.copy(deep=True)
            Stream_out_save = Stream_out_save.drop(["Length"], axis=1)
            merged_streams_list.append(Stream_out_save)

            Stream_out_PP.drop(
                columns=Stream_out_PP.columns.difference(
                    ["Length", "geometry", object_column_name]
                ),
                inplace=True,
            )
            terminus = Stream_out_PP.copy(
                deep=True
            )  # The intersection feature is replaced with the output GeoDataFrame

            terminus = terminus.buffer(0.1)
            ii = 0
            Stream_in_copy = stream_in.copy(deep=True)
            for index, row in Stream_in_copy.iterrows():
                for ii in range(0, len(terminus)):
                    within_loop = terminus.loc[ii].covers(
                        Stream_in_copy.loc[index, "geometry"]
                    )
                    ii = ii + 1
                    if within_loop == True:
                        try:
                            stream_in = stream_in.drop(labels=index, axis=0)
                        except KeyError:
                            continue
            terminus = terminus.to_frame(name="geometry")
            i = i + 1
            print("\n")
        except ValueError:
            break
    print(
        "\n A total of "
        + str(i)
        + " iterations were merged and saved into a single shapefile!"
    )
    merged_streams_gdf = geopandas.GeoDataFrame(
        pd.concat(merged_streams_list, ignore_index=True)
    )
    merged_streams_gdf.to_file(
        f"{outpath}streams_{topo_id}.shp", driver="ESRI Shapefile"
    )


def explode_multilinestrings(geometry):
    if isinstance(geometry, LineString):
        return [geometry]
    elif isinstance(geometry, MultiLineString):
        return list(geometry)
    else:
        raise ValueError("Input geometry must be LineString or MultiLineString")

def split_lines_at_junctions(gdf, outpath = None):
    exploded_lines = geopandas.GeoDataFrame(columns=gdf.columns)
   # encountered_geometries = set()  # Track encountered geometries to remove duplicates
    encountered_geometries = []
    for idx, row in gdf.iterrows():
        exploded_geoms = explode_multilinestrings(row.geometry)
        for geom in exploded_geoms:
            # Check if geometry is already encountered
            exploded_lines = exploded_lines.append({'geometry': geom}, ignore_index=True)

    uniqline = geopandas.GeoDataFrame(columns=gdf.columns)
    
    for line in exploded_lines["geometry"]:
        if not any(p.equals(line) for p in uniqline["geometry"]):
            uniqline = uniqline.append({'geometry': line}, ignore_index=True)
    
    split_lines = geopandas.GeoDataFrame(columns=gdf.columns)
    
    for index, line in uniqline.iterrows():
        start_point = line.geometry.coords[0]
        end_point = line.geometry.coords[-1]
    
        # Check if start or end points of the line intersect with any other line
        start_intersect = uniqline[uniqline.geometry.touches(Point(start_point)) & (uniqline.index != index)]
        end_intersect = uniqline[uniqline.geometry.touches(Point(end_point)) & (uniqline.index != index)]
        
        if not start_intersect.empty or not end_intersect.empty:
            # Split the line at start or end vertices
            split_line = [line.geometry]
            for idx, intersect_line in start_intersect.iterrows():
                split_line.append(intersect_line.geometry)
            for idx, intersect_line in end_intersect.iterrows():
                split_line.append(intersect_line.geometry)
            
            # Create separate lines for each segment
            for segment in split_line:
                split_lines = split_lines.append({'geometry': segment}, ignore_index=True)
        else:
            split_lines = split_lines.append(line)
    
    uniqline = geopandas.GeoDataFrame(columns=gdf.columns)
    
    for line in split_lines["geometry"]:
        if not any(p.equals(line) for p in uniqline["geometry"]):
            uniqline = uniqline.append({'geometry': line}, ignore_index=True)
    
    uniqline['Length'] = uniqline['geometry'].length
    uniqline['ID'] = ""
    uniqline.crs = gdf.crs
    if outpath is not None:
        uniqline.to_file(outpath)
    return uniqline