"""
This script uses libraries in arcgis to create connections from
a series of points (buildings) to the closest street
"""

import os
import cea.globalvar
import cea.inputlocator
from cea.interfaces.arcgis.modules import arcpy
import cea.config
import os

__author__ = "Jimeno A. Fonseca"
__copyright__ = "Copyright 2017, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Jimeno A. Fonseca"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


def calc_connectivity_network(path_arcgis_db, path_streets_shp, path_connection_point_buildings_shp,
                              path_potential_network):
    """
    This script outputs a potential network connecting a series of building points to the closest street network
    the street network is assumed to be a good path to the district heating or cooling network

    :param path_arcgis_db: path to default ArcGIS database
    :param path_streets_shp: path to street shapefile
    :param path_connection_point_buildings_shp: path to substations in buildings (or close by)
    :param path_potential_network: output path shapefile
    :return:
    """
    ## first add distribution network to each building form the roads

    arcpy.env.overwriteOutput = True
    spatialReference = arcpy.Describe(path_connection_point_buildings_shp).spatialReference
    # create temporary layers
    memorybuildings = path_arcgis_db + "\\" + "points"
    merge = path_arcgis_db + "\\" + "merge"
    Newlines = path_arcgis_db + "\\" + "linesToerase"
    Finallines = path_arcgis_db + "\\" + "final_line"
    near_table = path_arcgis_db + "\\" + "near_table"
    Newpoints = path_arcgis_db + '\\' + "New_points"
    # write building nodes to points
    arcpy.CopyFeatures_management(path_connection_point_buildings_shp, memorybuildings)
    # arcpy.Near_analysis(memorybuildings, path_streets_shp, location=True, angle=True)
    rank_count = 1 # TODO: set to 3 to include 3 surrounding streets
    arcpy.GenerateNearTable_analysis(memorybuildings, path_streets_shp, near_table, location=True, closest='ALL',
                                     closest_count=rank_count)
    arcpy.JoinField_management(near_table, "IN_FID", memorybuildings, "OBJECTID", ["Name"])
    arcpy.MakeXYEventLayer_management(near_table, "NEAR_X", "NEAR_Y", "New_Points_Layer", spatialReference)
    arcpy.FeatureClassToFeatureClass_conversion("New_Points_Layer", path_arcgis_db, "New_points")
    # FIXME: add lines from substations to streets by NEAR_RANK
    for i in range(rank_count):
        memorybuildings_base = path_arcgis_db + "\\" + "points_base"
        arcpy.CopyFeatures_management(path_connection_point_buildings_shp, memorybuildings_base)
        rank = i+1 #fixme: i + 1
        lines_to_substations = path_arcgis_db + "\\" + "line_to_substations_%s" % rank
        new_points_rank = path_arcgis_db + "\\" + "new_points_rank_%s" % rank

        arcpy.MakeFeatureLayer_management(Newpoints, "POINTS_layer")
        arcpy.SelectLayerByAttribute_management("POINTS_layer", "NEW_SELECTION", '"NEAR_RANK"=%s' %rank)
        arcpy.CopyFeatures_management("POINTS_layer", new_points_rank)
        arcpy.Append_management(new_points_rank, memorybuildings_base, "No_Test")
        arcpy.MakeFeatureLayer_management(memorybuildings_base, "POINTS_layer")
        arcpy.env.workspace = path_arcgis_db
        arcpy.PointsToLine_management(memorybuildings_base, lines_to_substations, "Name", "#", "NO_CLOSE") # FIXME: the problem is that when i = 2, the lines_to_substations are still messed up
        arcpy.Merge_management([path_streets_shp, lines_to_substations], merge)

    arcpy.FeatureToLine_management(merge, path_potential_network)  # necessary to match vertices


def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(scenario=config.scenario)

    path_streets_shp = locator.get_street_network()  # shapefile with the stations
    path_connection_point_buildings_shp = locator.get_connection_point()  # substation, it can be the centroid of the building
    path_potential_network = locator.get_connectivity_potential()  # shapefile, location of output.
    path_default_arcgis_db = os.path.expanduser(os.path.join('~', 'Documents', 'ArcGIS', 'Default.gdb'))
    calc_connectivity_network(path_default_arcgis_db, path_streets_shp, path_connection_point_buildings_shp,
                              path_potential_network)


if __name__ == '__main__':
    main(cea.config.Configuration())
