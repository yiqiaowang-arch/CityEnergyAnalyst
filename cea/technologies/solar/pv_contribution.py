import os
import warnings

import numpy as np
import pandas as pd

from cea.config import Configuration
from cea.constants import HOURS_IN_YEAR
from cea.inputlocator import InputLocator

_GRID_COMPONENTS: list[list[str]] = [
    [  # components in E_sys
        "GRID_a",  # Appliances
        "GRID_l",  # Lighting
        "GRID_v",  # Ventilation
        "GRID_ve",  # Elevators
        "GRID_data",  # Data centers
        "GRID_pro",  # Industrial processes
        "GRID_aux",  # Auxiliary
    ],
    [  # components in Qcs_sys
        "GRID_cdata",  # Data center cooling
        "GRID_cre",  # Refrigeration cooling
        "GRID_cs",  # Space cooling
    ],
    ["GRID_hs"],  # Heating (Qhs_sys)
    ["GRID_ww"],  # Domestic hot water (Qww_sys)
]
_GRID_COMPONENTS_FLAT: list[str] = [
    comp for sublist in _GRID_COMPONENTS for comp in sublist
]


def log_pv_contribution(
    locator: InputLocator,
    building_name: str,
    pv_codes: list[str],
    allowed_components: list[str] | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Calculate PV emission offsets for GRID electricity components and grid export.

    PV is allocated to building GRID components in the order specified by allowed_components.
    Emission offsets are calculated per component and for grid export.

    :param pv_codes: List of PV panel codes to process
    :type pv_codes: list[str]
    :param offset_building: If True, PV offsets building GRID demand before being exported.
        If False, all PV is treated as export to grid (no building offset).
    :type offset_building: bool
    :param allowed_components: Ordered list of GRID components that PV can offset.
        Valid values: ['GRID_a', 'GRID_l', 'GRID_v', 'GRID_ve', 'GRID_data', 'GRID_pro',
                        'GRID_aux', 'GRID_cdata', 'GRID_cre', 'GRID_hs', 'GRID_cs', 'GRID_ww']
        The ORDER determines allocation priority. PV offsets the first component fully before
        moving to the second, and so on.
        Only used when offset_building=True. If None or empty when offset_building=True,
        all GRID components will be used in default order.
    :type allowed_components: list[str] | None

    Notes
    -----
    - All emission offsets are negative by convention (emissions avoided)
    - Components are only offset if they exist in the demand file and have non-zero demand
    - Per-component offsets are stored in PV_{code}_offset_{component}_kgCO2e columns
    - Grid export offset is stored in PV_{code}_offset_grid_export_kgCO2e column
    - Total offset is sum of all per-component and export offsets

    Examples
    --------
    >>> # Offset lighting first, then appliances, export surplus
    >>> hourly_timeline.log_pv_contribution(
    ...     locator,
    ...     building_name='Building_1',
    ...     pv_codes=['PV1'],
    ...     offset_building=True,
    ...     allowed_components=['GRID_l', 'GRID_a'],
    ... )

    >>> # Export all PV to grid without offsetting building
    >>> hourly_timeline.log_pv_contribution(
    ...     locator,
    ...     building_name='Building_1',
    ...     pv_codes=['PV1'],
    ...     offset_building=False,
    ... )

    >>> # Offset building with all components, don't credit export (conservative)
    >>> hourly_timeline.log_pv_contribution(
    ...     locator,
    ...     building_name='Building_1',
    ...     pv_codes=['PV1'],
    ...     offset_building=True,
    ...     allowed_components=None,
    ... )
    """
    if len(pv_codes) == 0:
        raise ValueError("pv_codes must contain at least one PV panel code.")

    # Validate and prepare allowed components (preserves order)
    # Note: CEA's MultiChoiceParameter returns ALL choices when parameter is empty
    # So allowed_components will be a full list when user leaves pv-offset-allowance empty
    if allowed_components is None:
        ordered_components = list(_GRID_COMPONENTS_FLAT)
    else:
        # Keep only valid components, preserving user-specified order
        ordered_components = [
            c for c in allowed_components if c in _GRID_COMPONENTS_FLAT
        ]

        # Warn about invalid components
        invalid = set(allowed_components) - set(_GRID_COMPONENTS_FLAT)
        if invalid:
            warnings.warn(
                f"Invalid GRID components in pv-offset-allowance will be ignored: {invalid}. "
                f"Valid components are: {_GRID_COMPONENTS_FLAT}",
                RuntimeWarning,
            )

    demand_timeseries = load_demand_timeseries(locator, building_name)
    contribution_by_pv: dict[str, pd.DataFrame] = {}
    for pv_code in pv_codes:
        pv_hourly = load_pv_hourly(locator, building_name, pv_code)
        pv_offset = calculate_pv_offset(
            demand_timeseries, pv_hourly, ordered_components
        )
        contribution_by_pv[pv_code] = pv_offset
    return contribution_by_pv


def calculate_pv_offset(
    demand_timeseries: pd.DataFrame,
    pv_hourly: pd.DataFrame,
    ordered_components: list[str],
) -> pd.DataFrame:
    # check that pv_hourly has correct length as demand_timeseries
    if len(pv_hourly) != len(demand_timeseries):
        raise ValueError(
            f"PV hourly data length {len(pv_hourly)} does not match demand timeseries length {len(demand_timeseries)}."
        )
    if set(ordered_components) - set(_GRID_COMPONENTS_FLAT):
        raise ValueError(
            f"ordered_components contains invalid components: {set(ordered_components) - set(_GRID_COMPONENTS_FLAT)}. "
            f"Valid components are: {_GRID_COMPONENTS_FLAT}"
        )

    pv_offset = pd.DataFrame(index=demand_timeseries.index)

    pv_gen = pv_hourly["E_PV_gen_kWh"].to_numpy(dtype=float)
    pv_remaining = pv_gen.copy()

    for comp in ordered_components:
        comp_col_demand = f"{comp}_kWh"
        comp_col_pv_offset = f"E_PV_offset_{comp}_kWh"

        if comp_col_demand not in demand_timeseries.columns:
            continue

        comp_demand = demand_timeseries[comp_col_demand].to_numpy(dtype=float)
        pv_to_comp = np.minimum(pv_remaining, comp_demand)
        pv_offset.loc[:, comp_col_pv_offset] = pv_to_comp
        pv_remaining -= pv_to_comp

    pv_offset.loc[:, "E_PV_export_to_grid_kWh"] = pv_remaining
    pv_offset.loc[:, "E_PV_offset_total_kWh"] = pv_gen - pv_remaining
    return pv_offset


def load_pv_hourly(
    locator: InputLocator, building_name: str, pv_code: str
) -> pd.DataFrame:
    """Load hourly PV generation data for a given building and PV panel code, index named by 'hour'."""
    pv_results_path = locator.PV_results(building_name, pv_code)
    if not os.path.exists(pv_results_path):
        raise FileNotFoundError(
            f"PV results file not found for {pv_code} at {pv_results_path}."
        )
    df = pd.read_csv(pv_results_path, index_col=None)
    df.index.set_names("hour", inplace=True)
    if len(df) != HOURS_IN_YEAR:
        raise ValueError(
            f"PV results for {pv_code} should have {HOURS_IN_YEAR} rows, but got {len(df)} rows."
        )
    if "E_PV_gen_kWh" not in df.columns:
        raise ValueError(
            f"Expected column 'E_PV_gen_kWh' not found in PV results for {pv_code}."
        )
    return df


def load_demand_timeseries(locator: InputLocator, building_name: str) -> pd.DataFrame:
    """Load hourly demand timeseries for a given building, index named by 'hour'."""
    demand_results_path = locator.get_demand_results_file(building_name)
    if not os.path.exists(demand_results_path):
        raise FileNotFoundError(
            f"Demand results file not found at {demand_results_path}."
        )
    df = pd.read_csv(demand_results_path, index_col=None)
    df.index.set_names("hour", inplace=True)
    if len(df) != HOURS_IN_YEAR:
        raise ValueError(
            f"Demand results should have {HOURS_IN_YEAR} rows, but got {len(df)} rows."
        )
    return df


def main(config: Configuration) -> None:
    locator = InputLocator(config.scenario)
    buildings = config.PV_offset.buildings
    pv_codes = config.PV_offset.pv_codes
    allowed_components = config.PV_offset.pv_offset_allowance

    for building_name in buildings:
        contribution = log_pv_contribution(
            locator, building_name, pv_codes, allowed_components
        )
        
