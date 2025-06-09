"""embodied_modular.py
 A refactored, modular version of the CEA embodied-energy / grey-emissions workflow.

 The public entry-point is :pyfunc:`lca_embodied`.
 """

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

__all__ = [
    "lca_embodied",
]

# -----------------------------------------------------------------------------
# 1. I/O HELPERS
# -----------------------------------------------------------------------------

def load_building_inputs(locator: cea.inputlocator.InputLocator) -> Tuple[pd.DataFrame, gpd.GeoDataFrame]:
    """Load *architecture* and *zone geometry* tables."""
    architecture = pd.read_csv(locator.get_building_architecture())
    zone = gpd.read_file(locator.get_zone_geometry())
    return architecture, zone


def project_and_prepare_zone(zone: gpd.GeoDataFrame) -> pd.DataFrame:
    """Project *zone* GeoDataFrame to an equal-area CRS and add *footprint* & *perimeter* columns."""
    lat, lon = get_lat_lon_projected_shapefile(zone)
    projected = zone.to_crs(get_projected_coordinate_system(lat, lon))
    result = projected.assign(
        footprint=lambda d: d.area,
        perimeter=lambda d: d.length,
    ).drop(columns=["geometry"])
    return result


def load_surface_databases(locator: "cea.inputlocator.InputLocator") -> Dict[str, pd.DataFrame]:
    """Read the four envelope assembly catalogues into a dictionary."""
    return {
        "win": pd.read_csv(locator.get_database_assemblies_envelope_window()),
        "roof": pd.read_csv(locator.get_database_assemblies_envelope_roof()),
        "wall": pd.read_csv(locator.get_database_assemblies_envelope_wall()),
        "floor": pd.read_csv(locator.get_database_assemblies_envelope_floor()),
    }


# -----------------------------------------------------------------------------
# 2. TRANSFORMATIONS
# -----------------------------------------------------------------------------

_ENVELOPE_MAP = {
    # key     arch-col     ghg_col                bio_ghg_col                  life_col
    "win":  ("type_win",  "GHG_win_kgCO2m2",     "GHG_biogenic_win_kgCO2m2",  "Service_Life_win"),
    "roof": ("type_roof", "GHG_roof_kgCO2m2",    "GHG_biogenic_roof_kgCO2m2", "Service_Life_roof"),
    "wall": ("type_wall", "GHG_wall_kgCO2m2",    "GHG_biogenic_wall_kgCO2m2", "Service_Life_wall"),
    "floor":("type_floor","GHG_floor_kgCO2m2",   "GHG_biogenic_floor_kgCO2m2","Service_Life_floor"),
    # Special cases that reuse *wall* / *floor* catalogues but need renaming
    "base": ("type_base", "GHG_floor_kgCO2m2",   "GHG_biogenic_floor_kgCO2m2","Service_Life_floor"),
    "part": ("type_part", "GHG_wall_kgCO2m2",    "GHG_biogenic_wall_kgCO2m2", "Service_Life_wall"),
}


def build_surface_properties(
    architecture: pd.DataFrame, db: Dict[str, pd.DataFrame]
) -> pd.DataFrame:
    """Join *architecture* with the envelope catalogues in a data-driven way."""

    frames = []
    for key, (arch_col, ghg, bio, life) in _ENVELOPE_MAP.items():
        catalogue_key = (
            "wall"  if key == "part"
            else "floor" if key == "base"
            else key
        )
        merged = architecture.merge(
            db[catalogue_key],            # e.g. window catalogue
            left_on=arch_col,             # e.g. 'type_win'
            right_on="code"               # catalogue code
        )
        rename_dict = {
            ghg:  f"GHG_{key}_kgCO2m2",
            bio:  f"GHG_biogenic_{key}_kgCO2m2",
            life: f"Service_Life_{key}"
        }
        # leave only building name and lca intensity of "key" component in df
        df = (
            merged.rename(columns=rename_dict)
            [["name", f"GHG_{key}_kgCO2m2", f"GHG_biogenic_{key}_kgCO2m2", f"Service_Life_{key}"]]
        )
        frames.append(df)

    surface_props = reduce(lambda l, r: l.merge(r, on="name"), frames)
    return surface_props


# -----------------------------------------------------------------------------
# 3. GEOMETRY CALCULATIONS
# -----------------------------------------------------------------------------

def add_building_geometry(df: pd.DataFrame) -> pd.DataFrame:
    """Vectorised calculation of façade, window, and floor areas."""

    # Average WWR across four orientations
    mean_wwr = df[["wwr_south", "wwr_north", "wwr_west", "wwr_east"]].mean(axis=1)

    df = df.assign(
        windows_ag=lambda d: mean_wwr * d.perimeter * d.height_ag, # total window area
        area_walls_ext_ag=lambda d: d.perimeter * d.height_ag - d.windows_ag,
    )

    # Adjust for void decks
    empty_ratio: pd.Series = 1 - (
        (df.void_deck * (df.height_ag / df.floors_ag))
        / (df.area_walls_ext_ag + df.windows_ag)
    )
    df["windows_ag"] *= empty_ratio
    df["area_walls_ext_ag"] *= empty_ratio

    # Remaining geometry
    df = df.assign(
        area_walls_ext_bg=lambda d: d.perimeter * d.height_bg,
        floor_area_ag=lambda d: d.footprint * d.floors_ag,
        floor_area_bg=lambda d: d.footprint * d.floors_bg,
    )
    df["GFA_m2"] = df.floor_area_ag + df.floor_area_bg
    return df


# -----------------------------------------------------------------------------
# 4. LCA CORE (ported, lightly cleaned)
# -----------------------------------------------------------------------------

def calculate_contributions(df: pd.DataFrame, year: int) -> pd.DataFrame:  # noqa: C901 – long but self-contained
    """Compute yearly embodied GHG emissions & uptake for every building row."""

    # Years since construction & existence mask (<= 60y → 1, else 0)
    df = df.copy()
    df["delta_year"] = year - df["year"]
    df["confirm"] = (df["delta_year"] <= SERVICE_LIFE_OF_BUILDINGS).astype(int)

    # Helper to vectorise component formulae
    def _component(col_ghg: str, area: pd.Series, life: pd.Series):
        yearly: pd.Series = (
            (df[col_ghg] * area * np.ceil(SERVICE_LIFE_OF_BUILDINGS / life)) # ceil() in order not to underestimate impact
            / 1000.0  # kg → ton
            * df["confirm"] # if building is older than service life, the emission will be fully discounted.
        )
        return yearly

    # Areas
    wall_area_total = df["area_walls_ext_ag"] + df["area_walls_ext_bg"]
    floor_area_total = df["floor_area_ag"] + df["floor_area_bg"]

    # Embodied & biogenic carbon uptake per component
    df["emb_win"] = _component("GHG_win_kgCO2m2", df["windows_ag"], df["Service_Life_win"])
    df["upt_win"] = _component("GHG_biogenic_win_kgCO2m2", df["windows_ag"], df["Service_Life_win"])

    df["emb_wall"] = _component("GHG_wall_kgCO2m2", wall_area_total, df["Service_Life_wall"])
    df["upt_wall"] = _component("GHG_biogenic_wall_kgCO2m2", wall_area_total, df["Service_Life_wall"])

    df["emb_part"] = _component(
        "GHG_part_kgCO2m2",
        floor_area_total * CONVERSION_AREA_TO_FLOOR_AREA_RATIO,
        df["Service_Life_part"],
    )
    df["upt_part"] = _component(
        "GHG_biogenic_part_kgCO2m2",
        floor_area_total * CONVERSION_AREA_TO_FLOOR_AREA_RATIO,
        df["Service_Life_part"],
    )

    df["emb_floor"] = _component("GHG_floor_kgCO2m2", df["floor_area_ag"], df["Service_Life_floor"])
    df["upt_floor"] = _component("GHG_biogenic_floor_kgCO2m2", df["floor_area_ag"], df["Service_Life_floor"])

    df["emb_base"] = _component("GHG_base_kgCO2m2", df["floor_area_bg"], df["Service_Life_base"])
    df["upt_base"] = _component("GHG_biogenic_base_kgCO2m2", df["floor_area_bg"], df["Service_Life_base"])

    df["emb_roof"] = _component("GHG_roof_kgCO2m2", df["footprint"], df["Service_Life_roof"])
    df["upt_roof"] = _component("GHG_biogenic_roof_kgCO2m2", df["footprint"], df["Service_Life_roof"])

    # Technical systems – constant factors
    df["emb_system"] = (
        floor_area_total
        * EMISSIONS_EMBODIED_TECHNICAL_SYSTEMS
        * np.ceil(SERVICE_LIFE_OF_BUILDINGS / SERVICE_LIFE_OF_TECHNICAL_SYSTEMS)
        / 1000.0
        * df["confirm"]
    )

    # Aggregation – divide by 60-year amortisation
    emb_cols = [c for c in df.columns if c.startswith("emb_")]
    upt_cols = [c for c in df.columns if c.startswith("upt_")]

    df["GHG_sys_embodied_tonCO2yr"] = df[emb_cols].sum(axis=1) / SERVICE_LIFE_OF_BUILDINGS
    df["GHG_sys_uptake_tonCO2yr"] = df[upt_cols].sum(axis=1) / SERVICE_LIFE_OF_BUILDINGS

    df["GHG_sys_embodied_kgCO2m2yr"] = (
        df["GHG_sys_embodied_tonCO2yr"] * 1000.0 / df["GFA_m2"]
    )
    df["GHG_sys_uptake_kgCO2m2yr"] = (
        df["GHG_sys_uptake_tonCO2yr"] * 1000.0 / df["GFA_m2"]
    )

    keep = [
        "name",
        "GFA_m2",
        "GHG_sys_embodied_tonCO2yr",
        "GHG_sys_uptake_tonCO2yr",
        "GHG_sys_embodied_kgCO2m2yr",
        "GHG_sys_uptake_kgCO2m2yr",
    ] + emb_cols + upt_cols + ["emb_system"]
    return df[keep]


# -----------------------------------------------------------------------------
# 5. PUBLIC ORCHESTRATOR
# -----------------------------------------------------------------------------

def lca_embodied(year_to_calculate: int, locator: "cea.inputlocator.InputLocator") -> None:
    """Run the full embodied-GHG pipeline and save results to *locator.get_lca_embodied()*."""

    # 1. Load & prepare inputs
    architecture_df, zone_gdf = load_building_inputs(locator)
    zone_df = project_and_prepare_zone(zone_gdf)
    surface_db = load_surface_databases(locator)

    # 2. Join tables
    surface_props = build_surface_properties(architecture_df, surface_db)
    df = (
        zone_df
        .merge(surface_props, on="name")
        .merge(architecture_df, on="name")
        .pipe(add_building_geometry)
    )

    # 3. LCA core
    results = calculate_contributions(df, year_to_calculate)

    # 4. Write output
    results.to_csv(
        locator.get_lca_embodied(),
        index=False,
        float_format="%.2f",
        na_rep="nan",
    )

    print("Embodied-GHG calculation completed →", locator.get_lca_embodied())


# -----------------------------------------------------------------------------
# 6. CLI ENTRY-POINT (kept for compatibility)
# -----------------------------------------------------------------------------

def main(config: "cea.config.Configuration") -> None:  # pragma: no cover
    assert os.path.exists(config.scenario), f"Scenario not found: {config.scenario}"
    locator = cea.inputlocator.InputLocator(scenario=config.scenario)

    print("Running embodied_modular with scenario =", config.scenario)
    print("Year to calculate =", config.emissions.year_to_calculate)

    lca_embodied(year_to_calculate=config.emissions.year_to_calculate, locator=locator)


if __name__ == "__main__":  # pragma: no cover
    main(cea.config.Configuration())
