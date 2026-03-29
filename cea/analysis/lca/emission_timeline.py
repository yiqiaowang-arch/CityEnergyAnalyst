from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

import pandas as pd

from cea.analysis.lca.emission_timeline_backend import (
    COMPONENT_TO_SRC_COMPONENT,
    TECHNICAL_SYSTEM_COMPONENTS,
    LifecycleChange,
    OperationalSegment,
    ResolvedComponent,
    TimelineFrontendBase,
    _coerce_float,
    envelope_intensities_per_m2,
)
from cea.analysis.lca.hourly_operational_emission import OperationalHourlyTimeline
from cea.constants import (
    CONVERSION_AREA_TO_FLOOR_AREA_RATIO,
    EMISSIONS_EMBODIED_TECHNICAL_SYSTEMS,
    SERVICE_LIFE_OF_TECHNICAL_SYSTEMS,
)
from cea.datamanagement.database.envelope_lookup import EnvelopeLookup
from cea.demand.building_properties import BuildingProperties
from cea.utilities import epwreader

__author__ = "Yiqiao Wang, Zhongming Shi"
__copyright__ = "Copyright 2025, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Yiqiao Wang", "Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Reynold Mok"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


if TYPE_CHECKING:
    from cea.inputlocator import InputLocator


_COMPONENT_KEYS_BY_SRC_COMPONENT: dict[str, tuple[str, ...]] = {
    src_component: tuple(
        component
        for component, mapped_src_component in COMPONENT_TO_SRC_COMPONENT.items()
        if mapped_src_component == src_component
    )
    for src_component in dict.fromkeys(COMPONENT_TO_SRC_COMPONENT.values())
}


def get_component_quantities(
    building_properties: BuildingProperties,
    building_name: str,
) -> dict[str, float]:
    name = str(building_name)
    rc_model_props = building_properties.rc_model[name]
    envelope_props = building_properties.envelope[name]
    geometry_props = building_properties.geometry[name]

    surface_area: dict[str, float] = {}
    surface_area["Awall_ag"] = float(envelope_props["Awall_ag"])
    surface_area["Awall_bg"] = float(geometry_props["perimeter"]) * float(geometry_props["height_bg"])
    surface_area["Awall_part"] = float(rc_model_props["GFA_m2"]) * float(CONVERSION_AREA_TO_FLOOR_AREA_RATIO)
    surface_area["Awin_ag"] = float(envelope_props.get("Awin_ag", 0.0))
    surface_area["Aroof"] = float(envelope_props["Aroof"])
    surface_area["Aupperside"] = float(envelope_props.get("Aupperside", 0.0))
    surface_area["Aunderside"] = float(envelope_props.get("Aunderside", 0.0))

    if float(geometry_props["floors_bg"]) == 0 and float(geometry_props.get("void_deck", 0)) > 0:
        area_base = 0.0
    else:
        area_base = float(rc_model_props["footprint"])

    surface_area["Abase"] = float(area_base)
    surface_area["Afloor"] = max(
        0.0,
        float(rc_model_props["GFA_m2"]) - float(surface_area["Aunderside"]) - float(surface_area["Abase"]),
    )
    surface_area["Atechnical_systems"] = float(rc_model_props["GFA_m2"])
    return surface_area


def load_building_properties_for_emissions(
    locator: "InputLocator",
    buildings: list[str],
) -> BuildingProperties:
    weather_path = locator.get_weather_file()
    weather_data = epwreader.epw_reader(weather_path)[
        ["year", "drybulb_C", "wetbulb_C", "relhum_percent", "windspd_ms", "skytemp_C"]
    ]
    return BuildingProperties(locator, weather_data, buildings)


class NormalTimelineFrontend(TimelineFrontendBase):
    """Building-resolved frontend for the shared emission timeline backend.

    Normal keeps the building as the source of truth. It resolves the current
    building's envelope codes, areas, service lives, PV setup, and one
    carried-forward operational profile, then hands those as canonical changes
    and segments to the backend.
    """

    def __init__(
        self,
        *,
        building_properties: BuildingProperties,
        building_name: str,
        locator: "InputLocator",
        end_year: int,
        pv_codes: Sequence[str] | None = None,
    ) -> None:
        self.name = str(building_name)
        self.locator = locator
        self.envelope_lookup = EnvelopeLookup.from_locator(locator)
        self.geometry = building_properties.geometry[self.name]
        self.typology = building_properties.typology[self.name]
        self.envelope = building_properties.envelope[self.name]
        self.surface_area = get_component_quantities(building_properties, self.name)
        self._construction_year = int(self.typology["year"])
        self._end_year = int(end_year)
        if self._construction_year >= self._end_year:
            raise ValueError("The starting year must be less than the ending year.")
        self._window_code = str(self.envelope["type_win"])
        self._tech_emission_per_quarter = (
            float(EMISSIONS_EMBODIED_TECHNICAL_SYSTEMS) / len(TECHNICAL_SYSTEM_COMPONENTS)
        )
        self._pv_entries = self._load_pv_entries(list(pv_codes or []))

    def building_name(self) -> str:
        return self.name

    def start_year(self) -> int:
        return self._construction_year

    def end_year(self) -> int:
        return self._end_year

    def initial_change(self) -> LifecycleChange | None:
        components: list[ResolvedComponent] = []
        for src_component in ("wall", "base", "part", "roof", "floor"):
            components.extend(self._resolved_components_for_src_component(src_component))
        components.extend(self._resolved_window_components())
        components.extend(self._resolved_technical_components(TECHNICAL_SYSTEM_COMPONENTS))
        components.extend(self._resolved_pv_components())
        return LifecycleChange(
            year=self._construction_year,
            mode="add",
            components=tuple(components),
            note="Constructed",
        )

    def initial_due_years(self) -> dict[str, int]:
        due_years: dict[str, int] = {}
        for src_component in ("wall", "base", "part", "roof", "floor"):
            total_area = sum(
                float(self.surface_area.get(f"A{component}", 0.0))
                for component in _COMPONENT_KEYS_BY_SRC_COMPONENT[src_component]
            )
            if total_area <= 0.0:
                continue
            due_years[src_component] = self._construction_year + self._service_life_for_src_component(src_component)

        if float(self.surface_area.get("Awin_ag", 0.0)) > 0.0:
            due_years["win"] = self._construction_year + self._service_life_for_src_component("win")

        if float(self.surface_area.get("Atechnical_systems", 0.0)) > 0.0:
            for component in TECHNICAL_SYSTEM_COMPONENTS:
                due_years[component] = self._construction_year + int(SERVICE_LIFE_OF_TECHNICAL_SYSTEMS)

        for pv_key, entry in self._pv_entries.items():
            if float(entry["area_m2"]) > 0.0:
                due_years[pv_key] = self._construction_year + int(entry["lifetime"])

        return due_years

    def authored_changes_by_year(self) -> dict[int, tuple[LifecycleChange, ...]]:
        return {}

    def replacement_change(self, key: str, year: int) -> LifecycleChange | None:
        if key in TECHNICAL_SYSTEM_COMPONENTS:
            return LifecycleChange(
                year=year,
                mode="production_only",
                components=tuple(self._resolved_technical_components([key])),
                note="Service life reached: technical_systems",
            )

        if key.startswith("pv:"):
            entry = self._pv_entries.get(key)
            if entry is None or float(entry["area_m2"]) <= 0.0:
                return None
            pv_code = str(entry["pv_code"])
            return LifecycleChange(
                year=year,
                mode="add",
                components=(
                    ResolvedComponent(
                        component=str(entry["component"]),
                        area_m2=float(entry["area_m2"]),
                        production_per_area=float(entry["production_per_area"]),
                    ),
                ),
                note=f"Service life reached: PV_{pv_code} ({pv_code})",
            )

        if key == "win":
            components = self._resolved_window_components()
            if not components:
                return None
            return LifecycleChange(
                year=year,
                mode="replace",
                components=tuple(components),
                note=f"Service life reached: win_ag ({self._window_code})",
            )

        if key not in _COMPONENT_KEYS_BY_SRC_COMPONENT:
            return None

        components = self._resolved_components_for_src_component(key)
        if not components:
            return None
        code = str(self.envelope[f"type_{key}"])
        note = " | ".join(
            f"Service life reached: {component.component} ({code})"
            for component in components
        )
        return LifecycleChange(
            year=year,
            mode="replace",
            components=tuple(components),
            note=note,
        )

    def next_due_year(self, key: str, year: int) -> int | None:
        if key in TECHNICAL_SYSTEM_COMPONENTS:
            return int(year) + int(SERVICE_LIFE_OF_TECHNICAL_SYSTEMS)
        if key.startswith("pv:"):
            entry = self._pv_entries.get(key)
            if entry is None:
                return None
            return int(year) + int(entry["lifetime"])
        if key in _COMPONENT_KEYS_BY_SRC_COMPONENT or key == "win":
            return int(year) + self._service_life_for_src_component("win" if key == "win" else key)
        return None

    def operational_segments(self) -> tuple[OperationalSegment, ...]:
        operational = OperationalHourlyTimeline.from_result(self.locator, self.name)
        operational_timeseries = operational.operational_emission_timeline.copy()
        if "date" in operational_timeseries.columns:
            operational_timeseries = operational_timeseries.drop(columns=["date"])
        yearly_sum = operational_timeseries.sum(axis=0)
        column_values = {
            str(column): float(value)
            for column, value in yearly_sum.items()
            if isinstance(column, str) and column.endswith("_kgCO2e")
        }
        if not column_values:
            return ()
        return (
            OperationalSegment(
                start_year=self._construction_year,
                end_year=self._end_year,
                column_values=column_values,
            ),
        )

    def demolition_change(self) -> LifecycleChange | None:
        return None

    def _service_life_for_src_component(self, src_component: str) -> int:
        if src_component == "technical_systems":
            return int(SERVICE_LIFE_OF_TECHNICAL_SYSTEMS)
        code = self._window_code if src_component == "win" else str(self.envelope[f"type_{src_component}"])
        lifetime_any = self.envelope_lookup.get_item_value(code=code, field="Service_Life")
        if lifetime_any is None:
            raise ValueError(f"Envelope database returned None for Service_Life for item {code}.")
        lifetime = int(lifetime_any)
        if lifetime <= 0:
            raise ValueError(f"Lifetime must be positive for item {code}.")
        return lifetime

    def _resolved_components_for_src_component(self, src_component: str) -> list[ResolvedComponent]:
        if src_component == "technical_systems":
            return self._resolved_technical_components(TECHNICAL_SYSTEM_COMPONENTS)

        code = self._window_code if src_component == "win" else str(self.envelope[f"type_{src_component}"])
        production, demolition, biogenic = envelope_intensities_per_m2(
            self.envelope_lookup,
            code=code,
        )
        return [
            ResolvedComponent(
                component=component,
                area_m2=float(self.surface_area.get(f"A{component}", 0.0)),
                production_per_area=production,
                demolition_per_area=demolition,
                biogenic_per_area=biogenic,
            )
            for component in _COMPONENT_KEYS_BY_SRC_COMPONENT[src_component]
            if float(self.surface_area.get(f"A{component}", 0.0)) > 0.0
        ]

    def _resolved_window_components(self) -> list[ResolvedComponent]:
        return self._resolved_components_for_src_component("win")

    def _resolved_technical_components(
        self,
        components: Sequence[str],
    ) -> list[ResolvedComponent]:
        area = float(self.surface_area.get("Atechnical_systems", 0.0))
        if area <= 0.0:
            return []
        return [
            ResolvedComponent(
                component=str(component),
                area_m2=area,
                production_per_area=self._tech_emission_per_quarter,
            )
            for component in components
        ]

    def _resolved_pv_components(self) -> list[ResolvedComponent]:
        return [
            ResolvedComponent(
                component=str(entry["component"]),
                area_m2=float(entry["area_m2"]),
                production_per_area=float(entry["production_per_area"]),
            )
            for entry in self._pv_entries.values()
            if float(entry["area_m2"]) > 0.0
        ]

    def _load_pv_entries(self, pv_codes: list[str]) -> dict[str, dict[str, float | int | str]]:
        if not pv_codes:
            return {}

        pv_db = pd.read_csv(
            self.locator.get_db4_components_conversion_conversion_technology_csv("PHOTOVOLTAIC_PANELS"),
            index_col="code",
        )
        embodied_column = (
            "module_embodied_kgCO2m2"
            if "module_embodied_kgCO2m2" in pv_db.columns
            else "module_embodied_kgco2m2"
        )
        entries: dict[str, dict[str, float | int | str]] = {}
        for pv_code in pv_codes:
            if pv_code not in pv_db.index:
                raise ValueError(f"PV type {pv_code} not found in the PV database.")

            district_pv_area = pd.read_csv(self.locator.PV_total_buildings(pv_code), index_col="name")
            try:
                pv_area = cast(float, district_pv_area.at[self.name, "area_PV_m2"])
            except KeyError:
                pv_area = 0.0

            lifetime = cast(int, pv_db.loc[pv_code, "LT_yr"])
            if int(lifetime) <= 0:
                raise ValueError(f"PV lifetime must be positive for panel type {pv_code}.")

            component = f"PV_{pv_code}"
            entries[f"pv:{pv_code}"] = {
                "pv_code": pv_code,
                "component": component,
                "area_m2": float(pv_area),
                "lifetime": int(lifetime),
                "production_per_area": _coerce_float(pv_db.loc[pv_code, embodied_column]),
            }
        return entries
