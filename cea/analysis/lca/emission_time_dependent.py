import os
import warnings

import pandas as pd
from pandas.errors import EmptyDataError, ParserError

from cea.analysis.lca.emission_timeline_backend import (
    EmissionTimelineBackend,
    feedstock_policies_from_config,
    resolve_emissions_end_year,
)
from cea.analysis.lca.emission_timeline import (
    NormalTimelineFrontend,
    load_building_properties_for_emissions,
)
from cea.analysis.lca.hourly_operational_emission import OperationalHourlyTimeline
from cea.config import Configuration
from cea.inputlocator import InputLocator

__author__ = "Yiqiao Wang, Zhongming Shi"
__copyright__ = "Copyright 2025, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Yiqiao Wang", "Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Reynold Mok"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


def _load_grid_emission_intensity_override(config: Configuration):
    """Load and validate an optional external CSV for grid carbon intensity.

    Returns a tuple (override, values) where:
    - override: bool indicating whether to override GRID intensity timeline
    - values: numpy array of length 8760 if override is True; otherwise None

    Rules:
    - If grid_carbon_intensity_dataset_csv is not provided, returns (False, None)
    - If provided, csv_carbon_intensity_column_name must also be provided
    - Accepts 8760 rows; if 8784 rows, drops Feb 29 (hours 1416..1439) and warns
    - Raises with clear messages on missing file/column, parse errors, or NaNs
    """
    emissions_cfg = getattr(config, 'emissions')
    intensity_csv_path = getattr(emissions_cfg, 'grid_carbon_intensity_dataset_csv', None)
    intensity_column_name = getattr(emissions_cfg, 'csv_carbon_intensity_column_name', None)

    if not intensity_csv_path:
        return False, None
    if not intensity_column_name:
        raise ValueError(
            "If grid_carbon_intensity_dataset_csv is provided, csv_carbon_intensity_column_name must also be provided."
        )

    try:
        series: pd.Series = pd.read_csv(
            intensity_csv_path,
            usecols=[intensity_column_name],
            dtype={intensity_column_name: "float64"},
        )[intensity_column_name]
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"Could not find the provided CSV file '{intensity_csv_path}'."
        ) from e
    except PermissionError as e:
        raise PermissionError(
            f"Permission denied reading '{intensity_csv_path}'. It may be open in another program."
        ) from e
    except (EmptyDataError, ParserError, OSError, ValueError) as e:
        raise ValueError(
            f"Could not parse '{intensity_csv_path}' with column '{intensity_column_name}': {e}"
        ) from e

    n = len(series)
    if n == 8760:
        pass
    elif n == 8784:
        # Drop Feb 29 (hours 1416..1439) for a non-leap-year 8760-hour year
        series = series.drop(index=range(1416, 1440))
        warnings.warn(
            "Emission intensity CSV has 8784 rows; dropped Feb 29 to produce 8760 rows.",
            RuntimeWarning,
        )
    else:
        raise ValueError(
            f"Emission intensity dataset CSV file must have 8760 or 8784 rows, but has {n} rows."
        )

    if series.isna().any():
        na_count = int(series.isna().sum())
        raise ValueError(
            f"Emission intensity contains {na_count} NaN values; please clean or impute the data."
        )

    return True, series.to_numpy(dtype=float)


def operational_hourly(config: Configuration) -> None:
    locator = InputLocator(config.scenario)
    emissions_cfg = config.emissions
    buildings = emissions_cfg.buildings

    # Validate zone geometry EARLY before expensive calculations
    print("Validating zone geometry...")
    import geopandas as gpd
    from cea.utilities.standardize_coordinates import get_lat_lon_projected_shapefile
    try:
        zone_gdf = gpd.read_file(locator.get_zone_geometry())
        # This will raise ValueError with detailed message if geometries are invalid
        get_lat_lon_projected_shapefile(zone_gdf)
        print("Zone geometry validation passed.")
    except ValueError as e:
        print(f"ERROR: {e}")
        raise

    # Check PV requirements BEFORE processing any buildings
    consider_pv = emissions_cfg.include_pv
    pv_codes: list[str] = [] # prevent unbound variable error later in apply_pv_offsetting
    if consider_pv:
        pv_codes = emissions_cfg.pv_codes
        first_building = buildings[0]

        # Check which panels are missing
        missing_panels = []
        for pv_code in (pv_codes if pv_codes else []):
            pv_path = locator.PV_results(first_building, pv_code)
            if not os.path.exists(pv_path):
                missing_panels.append(pv_code)

        if missing_panels:
            missing_list = ', '.join(missing_panels)
            error_msg = (
                f"PV electricity results missing for panel type(s): {missing_list}. "
                f"Please run the 'photovoltaic (PV) panels' script first to generate PV potential results for these panel types."
            )
            print(f"ERROR: {error_msg}")
            raise FileNotFoundError(error_msg)

    building_properties = load_building_properties_for_emissions(locator, buildings)
    results: list[tuple[str, pd.DataFrame]] = []
    # Load optional GRID carbon intensity override once for all buildings
    override_grid_emission, grid_emission_final_g = _load_grid_emission_intensity_override(config)
    grid_emission_final = grid_emission_final_g / 1000.0 if grid_emission_final_g is not None else None  # convert g to kg
    for building in buildings:
        bpr = building_properties[building]
        hourly_timeline = OperationalHourlyTimeline(locator, bpr)

        if override_grid_emission and grid_emission_final is not None:
            hourly_timeline.emission_intensity_timeline["GRID"] = grid_emission_final

        hourly_timeline.calculate_operational_emission()

        if consider_pv:
            hourly_timeline.apply_pv_offsetting(pv_codes)

        hourly_timeline.save_results()
        print(
            f"Hourly operational emissions for {building} calculated and saved in: {locator.get_lca_operational_hourly_building(building)}."
        )
        results.append((building, hourly_timeline.operational_emission_timeline_extended))

    # df_by_building = to_ton(sum_by_building(results))
    backend = EmissionTimelineBackend()
    df_by_building = backend.aggregate_by_building(results)
    # df_by_hour = to_ton(sum_by_index([df for _, df in results]))
    df_by_hour = backend.aggregate_by_index([df for _, df in results])
    df_by_building.to_csv(locator.get_total_yearly_operational_building(), float_format='%.2f')
    df_by_hour.to_csv(locator.get_total_yearly_operational_hour(), index=False, float_format='%.2f')
    print(
        f"District-level operational emissions saved in: {locator.get_lca_emissions_results_folder()}"
    )


def total_yearly(config: Configuration) -> None:
    locator = InputLocator(scenario=config.scenario)
    emissions_cfg = config.emissions
    buildings = emissions_cfg.buildings
    end_year = resolve_emissions_end_year(emissions_cfg)

    # Validate zone geometry EARLY before expensive calculations
    print("Validating zone geometry...")
    import geopandas as gpd
    from cea.utilities.standardize_coordinates import get_lat_lon_projected_shapefile
    try:
        zone_gdf = gpd.read_file(locator.get_zone_geometry())
        # This will raise ValueError with detailed message if geometries are invalid
        get_lat_lon_projected_shapefile(zone_gdf)
        print("Zone geometry validation passed.")
    except ValueError as e:
        print(f"ERROR: {e}")
        raise

    # Check PV requirements BEFORE processing any buildings
    consider_pv: bool = getattr(emissions_cfg, "include_pv", False)
    pv_codes: list[str] = []

    if consider_pv:
        # Get PV codes from configuration (user must specify which PV types to include)
        pv_codes = getattr(emissions_cfg, "pv_codes", [])

        if not pv_codes:
            print("  Warning: include-pv is True but no PV codes specified in pv-codes parameter.")
            print("           No PV embodied emissions will be included.")
            consider_pv = False
        else:
            # Validate that results exist for all configured PV codes
            missing_panels = []
            for pv_code in pv_codes:
                pv_total_path = locator.PV_total_buildings(pv_code)
                if not os.path.exists(pv_total_path):
                    missing_panels.append(pv_code)

            if missing_panels:
                missing_list = ', '.join(missing_panels)
                error_msg = (
                    f"PV electricity results missing for panel type(s): {missing_list}. "
                    f"Please run the 'photovoltaic (PV) panels' script first to generate PV potential results for these panel types."
                )
                print(f"ERROR: {error_msg}")
                raise FileNotFoundError(error_msg)

            print(f"  Including PV life cycle emissions for panel types: {', '.join(pv_codes)}")

    building_properties = load_building_properties_for_emissions(locator, buildings)
    backend = EmissionTimelineBackend(
        feedstock_policies=feedstock_policies_from_config(emissions_cfg),
    )
    results: list[tuple[str, pd.DataFrame]] = []
    for building in buildings:
        frontend = NormalTimelineFrontend(
            building_properties=building_properties,
            building_name=building,
            locator=locator,
            end_year=end_year,
            pv_codes=pv_codes if consider_pv else [],
        )
        timeline = backend.build(frontend)
        timeline_to_save = timeline.copy()
        timeline_to_save["name"] = building
        timeline_to_save.to_csv(locator.get_lca_timeline_building(building), float_format='%.2f')
        print(
            f"Emission timeline for {building} calculated and saved in: {locator.get_lca_timeline_building(building)}."
        )
        results.append((building, timeline))
    #
    # df_by_building = to_ton(sum_by_building(results))
    df_by_building = backend.aggregate_by_building(results)
    # df_by_year = to_ton(sum_by_index([df for _, df in results]))
    df_by_year = backend.aggregate_by_index([df for _, df in results])
    df_by_building.to_csv(locator.get_total_emissions_building_year_end(year_end=end_year), float_format='%.2f')
    df_by_year.to_csv(locator.get_total_emissions_timeline_year_end(year_end=end_year), index=False, float_format='%.2f')
    print(
        f"District-level total emissions saved in: {locator.get_lca_timeline_folder()}"
    )


def to_ton(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a dataframe in kgCO2e to tonCO2e by dividing all values by 1000, and also rename the columns by changing 'kgCO2e' to 'tonCO2e'.

    :param df: A dataframe with values in kgCO2e.
    :type df: pd.DataFrame
    :return: A dataframe with values in tonCO2e.
    :rtype: pd.DataFrame
    """
    df_ton = df / 1000.0
    df_ton.columns = df_ton.columns.str.replace("kgCO2e", "tonCO2e")
    return df_ton


def main(config: Configuration) -> None:
    operational_hourly(config)
    total_yearly(config)


if __name__ == "__main__":
    main(Configuration())
