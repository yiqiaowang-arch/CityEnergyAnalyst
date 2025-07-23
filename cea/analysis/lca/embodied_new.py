from __future__ import annotations

import os
from functools import reduce
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import geopandas as gpd

import cea.config
import cea.inputlocator
from cea.constants import (
    SERVICE_LIFE_OF_BUILDINGS,
    SERVICE_LIFE_OF_TECHNICAL_SYSTEMS,
    CONVERSION_AREA_TO_FLOOR_AREA_RATIO,
    EMISSIONS_EMBODIED_TECHNICAL_SYSTEMS,
)
from cea.utilities.standardize_coordinates import (
    get_lat_lon_projected_shapefile,
    get_projected_coordinate_system,
)


def lca_embodied(year_to_calculate: int, locator: cea.inputlocator.InputLocator) -> None:
    """_summary_

    :param year_to_calculate: The "current year" when this LCA is calculated. Buildings components
        older than its lifetime but is still used will have no contribution to the embodied emission. 

        For example, if current year is 2025, a building is built in 1900, and its original roof has a 
        lifetime of 60 years.
        Since 60 < 2025-1900 = 125, the embodied emission of the roof is not anymore accounted in this building.

    :type year_to_calculate: int
    :param locator: the locator object.
    :type locator: cea.inputlocator.InputLocator
    """
    # 1. load building data and prepare inputs
    zone_df = get_building_geometry(locator)
    architecture_df = pd.read_csv(locator.get_building_architecture())
    component_df = read_components(locator)
    material_df = pd.read_csv(os.path.join(os.getcwd(), "data", "materials.csv"))

    return

def get_building_geometry(locator: cea.inputlocator.InputLocator) -> pd.DataFrame:
    """This function reads and reprojects the building geometry.

    :param locator: the cea locator object.
    :type locator: cea.inputlocator.InputLocator
    :return: A dataframe containing the geometrical information including 
        footprint area and perimeter.
    :rtype: pd.DataFrame
    """
    zone_gdf = gpd.read_file(locator.get_zone_geometry())
    lat, lon = get_lat_lon_projected_shapefile(zone_gdf)
    zone_gdf_proj = zone_gdf.to_crs(get_projected_coordinate_system(lat, lon))
    zone_df = zone_gdf_proj.assign(
        footprint=lambda d: d.area,
        perimeter=lambda d: d.length,
    ).drop(columns=["geometry"])

    return pd.DataFrame(zone_df)

def read_components(locator: cea.inputlocator.InputLocator) -> pd.DataFrame:
    """reads the composition of components from the database: 
    roof, wall, window, base, floor, partition.
    They are read separately from different files and output to one single dataframe.

    :param locator: the cea locator object
    :type locator: cea.inputlocator.InputLocator
    :return: one dataframe containing all types of component construction. 
        For each construction, one or more rows of materials, 
        as well as their thickness and emission intensity in different phases 
        are stored as separate columns.
    :rtype: pd.DataFrame
    """
    roofs = pd.read_csv(locator.get_database_assemblies_envelope_roof())
    walls = pd.read_csv(locator.get_database_assemblies_envelope_wall())
    windows = pd.read_csv(locator.get_database_assemblies_envelope_window())
    floors = pd.read_csv(locator.get_database_assemblies_envelope_floor())
    # add all components together to a long dataframe
    


def main(config: "cea.config.Configuration") -> None:  # pragma: no cover
    assert os.path.exists(config.scenario), f"Scenario not found: {config.scenario}"
    locator = cea.inputlocator.InputLocator(scenario=config.scenario)

    print("Running embodied_modular with scenario =", config.scenario)
    print("Year to calculate =", config.emissions.year_to_calculate)

    lca_embodied(year_to_calculate=config.emissions.year_to_calculate, locator=locator)


if __name__ == "__main__":  # pragma: no cover
    main(cea.config.Configuration())