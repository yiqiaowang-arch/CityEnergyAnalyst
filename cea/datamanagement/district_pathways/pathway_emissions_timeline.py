"""District pathway emissions timeline frontend built on the shared backend."""

from __future__ import annotations

import os
from bisect import bisect_right
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, cast

import geopandas as gpd
import numpy as np
import pandas as pd

from cea.analysis.lca.emission_timeline import (
    get_component_quantities,
    load_building_properties_for_emissions,
)
from cea.analysis.lca.emission_timeline_backend import (
    COMPONENT_TO_SRC_COMPONENT,
    TECHNICAL_SYSTEM_COMPONENTS,
    EmissionTimelineBackend,
    LifecycleChange,
    OperationalSegment,
    ResolvedComponent,
    TimelineFrontendBase,
    emission_timeline_columns,
    envelope_intensities_per_m2,
    feedstock_policies_from_config,
    finalise_emission_timeline_dataframe,
    resolve_emissions_end_year,
)
from cea.config import Configuration
from cea.constants import (
    EMISSIONS_EMBODIED_TECHNICAL_SYSTEMS,
    SERVICE_LIFE_OF_TECHNICAL_SYSTEMS,
)
from cea.datamanagement.database.envelope_lookup import EnvelopeLookup
from cea.datamanagement.district_pathways.pathway_years import (
    ensure_state_years_exist,
    get_building_construction_years,
    get_required_state_years,
)
from cea.inputlocator import InputLocator

__author__ = "Yiqiao Wang, Zhongming Shi"
__copyright__ = "Copyright 2026, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Yiqiao Wang", "Zhongming Shi"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Reynold Mok"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


MATERIAL_SRC_COMPONENTS: set[str] = {"wall", "roof", "base", "floor", "part"}
ASSEMBLY_SRC_COMPONENTS: set[str] = {"win"}
LAYERED_COMPONENT_TO_DB: dict[str, str] = {
    "wall": "wall",
    "roof": "roof",
    "base": "floor",
    "floor": "floor",
    "part": "wall",
}
COMP_AREA_MAP: list[tuple[str, str, str]] = [
    (component, src_component, f"A{component}")
    for component, src_component in COMPONENT_TO_SRC_COMPONENT.items()
    if src_component in MATERIAL_SRC_COMPONENTS
    or src_component in ASSEMBLY_SRC_COMPONENTS
    or src_component == "technical_systems"
]
CONSTRUCTION_TYPE_WIN_FIELD = "type_win"
SUPPLY_TYPE_FIELDS: tuple[str, ...] = (
    "supply_type_hs",
    "supply_type_cs",
    "supply_type_dhw",
    "supply_type_el",
)
TECH_SYSTEM_COMPONENT_BY_SUPPLY_FIELD: dict[str, str] = {
    "supply_type_hs": "technical_system_hs",
    "supply_type_cs": "technical_system_cs",
    "supply_type_dhw": "technical_system_dhw",
    "supply_type_el": "technical_system_el",
}


def _parse_code_change(value: Any) -> tuple[str, str] | None:
    if isinstance(value, tuple) and len(value) == 2:
        old_val, new_val = value
        return str(old_val), str(new_val)
    return None


@dataclass(frozen=True)
class MaterialLayer:
    """One raw pathway material-layer slot.

    This stays intentionally small:

    - `name` identifies the material entry in `MATERIALS.csv`
    - `thickness_m` stores the pathway-authored geometric thickness

    It does not cache embodied intensities. Those depend on the currently
    loaded materials database and its unit semantics, so pathway resolves them
    later when translating layers into backend-facing `ResolvedComponent`
    objects.
    """

    name: str
    thickness_m: float


def _empty_layers() -> list[MaterialLayer]:
    return [MaterialLayer(name="", thickness_m=0.0) for _ in range(3)]


def _coerce_3_layers(layers: list[MaterialLayer] | None) -> list[MaterialLayer]:
    src = list(layers) if layers else []
    out: list[MaterialLayer] = []
    for idx in range(3):
        if idx < len(src):
            name = str(src[idx].name or "").strip()
            thickness = float(src[idx].thickness_m)
            if thickness < 0.0:
                thickness = 0.0
            if not name:
                thickness = 0.0
            out.append(MaterialLayer(name=name, thickness_m=thickness))
        else:
            out.append(MaterialLayer(name="", thickness_m=0.0))
    return out


def _is_layer_active(layer: MaterialLayer) -> bool:
    return bool(str(layer.name).strip()) and float(layer.thickness_m) > 0.0


def _active_layers(layers: list[MaterialLayer]) -> list[MaterialLayer]:
    return [layer for layer in layers if _is_layer_active(layer)]


@dataclass(frozen=True)
class ArchetypeChangesTimeline:
    """Cumulative pathway snapshots plus year-local authored events.

    Pathway edits are stored as year-based deltas, but embodied calculations
    often need the effective construction state *at* a specific year. This
    object therefore keeps both:

    - snapshot-style `..._at_or_before(year)` accessors for effective state
    - `..._for_year(year)` accessors for authored events happening in that year
    """

    years_sorted: list[int]
    layers_at_year: dict[int, dict[str, dict[str, list[MaterialLayer]]]]
    layer_events_by_year: dict[int, dict[str, dict[str, list[tuple[str, MaterialLayer]]]]]
    construction_types_at_year: dict[int, dict[str, dict[str, str]]]
    construction_type_events_by_year: dict[int, dict[str, dict[str, tuple[str, str]]]]

    def layers_snapshot_at_or_before(self, year: int) -> dict[str, dict[str, list[MaterialLayer]]]:
        year_int = int(year)
        if year_int in self.layers_at_year:
            return dict(self.layers_at_year.get(year_int, {}) or {})
        prior_years = [value for value in self.years_sorted if value <= year_int]
        if not prior_years:
            return {}
        return dict(self.layers_at_year.get(prior_years[-1], {}) or {})

    def events_for_year(self, year: int) -> dict[str, dict[str, list[tuple[str, MaterialLayer]]]]:
        return dict(self.layer_events_by_year.get(int(year), {}) or {})

    def construction_types_snapshot_at_or_before(self, year: int) -> dict[str, dict[str, str]]:
        year_int = int(year)
        if year_int in self.construction_types_at_year:
            return dict(self.construction_types_at_year.get(year_int, {}) or {})
        prior_years = [value for value in self.years_sorted if value <= year_int]
        if not prior_years:
            return {}
        return dict(self.construction_types_at_year.get(prior_years[-1], {}) or {})

    def construction_type_events_for_year(self, year: int) -> dict[str, dict[str, tuple[str, str]]]:
        return dict(self.construction_type_events_by_year.get(int(year), {}) or {})


@dataclass(frozen=True)
class EffectiveBuildingState:
    """One cached effective embodied state for a building at a snapshot year.

    This is the pathway-side equivalent of a small pre-resolved building DB
    slice. It keeps the embodied-relevant state that the frontend repeatedly
    needs while assembling lifecycle changes:

    - effective window code
    - effective supply codes
    - effective 3-layer snapshots by source component
    - resolved window lifetime for that effective window code
    """

    window_code: str
    window_lifetime: int
    supply_codes_by_field: dict[str, str]
    layers_by_src_component: dict[str, list[MaterialLayer]]


class EffectiveBuildingStateCache:
    """Carry forward pathway snapshot state for one building construction type.

    The cache materialises effective state only for key snapshot years, then
    serves intermediate years by carrying forward the last known snapshot.
    This keeps pathway logic readable without repeatedly reading state DB files
    or scattering many independent `*_at(year)` lookups throughout the
    frontend.
    """

    def __init__(
        self,
        *,
        const_type: str,
        construction_year: int,
        years_sorted: Sequence[int],
        archetype_timeline: ArchetypeChangesTimeline,
        envelope_lookup: EnvelopeLookup,
        win_lifetime_default: int,
    ) -> None:
        self.const_type = str(const_type)
        self.archetype_timeline = archetype_timeline
        self.envelope_lookup = envelope_lookup
        self.win_lifetime_default = int(win_lifetime_default)
        self.snapshot_years = sorted({int(construction_year), *[int(year) for year in years_sorted]})
        self.states_by_snapshot_year = {
            year: self._build_state(year)
            for year in self.snapshot_years
        }
        self.has_embodied_state = self.const_type in self.archetype_timeline.layers_snapshot_at_or_before(
            int(construction_year)
        )

    def at(self, year: int) -> EffectiveBuildingState:
        year_int = int(year)
        idx = bisect_right(self.snapshot_years, year_int) - 1
        if idx < 0:
            raise ValueError(
                f"No effective pathway state is available at or before year {year_int} "
                f"for construction type '{self.const_type}'."
            )
        return self.states_by_snapshot_year[self.snapshot_years[idx]]

    def _build_state(self, year: int) -> EffectiveBuildingState:
        codes_snapshot = self.archetype_timeline.construction_types_snapshot_at_or_before(year)
        construction_codes = dict(codes_snapshot.get(self.const_type, {}) or {})
        window_code = str(construction_codes.get(CONSTRUCTION_TYPE_WIN_FIELD, ""))
        lifetime_any = self.envelope_lookup.get_item_value(window_code, "Service_Life") if window_code else None
        window_lifetime = int(lifetime_any) if lifetime_any is not None else self.win_lifetime_default
        if window_lifetime <= 0:
            raise ValueError(
                f"Invalid Service_Life={window_lifetime} for component 'win' "
                f"(const_type={self.const_type}, year={int(year)})."
            )

        layers_snapshot = self.archetype_timeline.layers_snapshot_at_or_before(year)
        raw_layers = dict(layers_snapshot.get(self.const_type, {}) or {})
        layers_by_src_component = {
            str(component): _coerce_3_layers(list(layers))
            for component, layers in raw_layers.items()
        }
        supply_codes_by_field = {
            field: str(construction_codes.get(field, ""))
            for field in SUPPLY_TYPE_FIELDS
        }
        return EffectiveBuildingState(
            window_code=window_code,
            window_lifetime=window_lifetime,
            supply_codes_by_field=supply_codes_by_field,
            layers_by_src_component=layers_by_src_component,
        )


def prepare_pathway_archetype_timeline(
    *,
    years_sorted: list[int],
    base_year: int | None,
    log_data: Mapping[int, dict[str, Any]],
    archetype_layers: dict[str, dict[str, list[MaterialLayer]]],
    archetype_construction_types: dict[str, dict[str, str]],
) -> ArchetypeChangesTimeline:
    layer_events_by_year: dict[int, dict[str, dict[str, list[tuple[str, MaterialLayer]]]]] = {}
    layers_at_year: dict[int, dict[str, dict[str, list[MaterialLayer]]]] = {}
    construction_at_year: dict[int, dict[str, dict[str, str]]] = {}
    construction_events_by_year: dict[int, dict[str, dict[str, tuple[str, str]]]] = {}

    if not years_sorted:
        return ArchetypeChangesTimeline(
            years_sorted=[],
            layers_at_year={},
            layer_events_by_year={},
            construction_types_at_year={},
            construction_type_events_by_year={},
        )

    effective_years = list(years_sorted)
    if base_year is not None and int(base_year) < int(min(years_sorted)):
        effective_years = sorted({int(base_year), *[int(year) for year in years_sorted]})
        base = int(base_year)
        layer_events_by_year[base] = {}
        layers_at_year[base] = {
            archetype: {component: layers[:] for component, layers in components.items()}
            for archetype, components in archetype_layers.items()
        }
        construction_events_by_year[base] = {}
        construction_at_year[base] = {
            archetype: dict(codes)
            for archetype, codes in archetype_construction_types.items()
        }

    for year in effective_years:
        entry = log_data.get(year, {}) or {}
        year_mods = entry.get("modifications", {}) or {}
        year_layer_events: dict[str, dict[str, list[tuple[str, MaterialLayer]]]] = {}
        year_construction_events: dict[str, dict[str, tuple[str, str]]] = {}

        for archetype, components in year_mods.items():
            archetype_name = str(archetype)
            if archetype_name not in archetype_layers:
                continue
            for component, patch in (components or {}).items():
                component_name = str(component)
                if component_name == "construction_type":
                    patch_dict = dict(patch or {})
                    current = archetype_construction_types.setdefault(archetype_name, {})
                    for key_raw, value_raw in patch_dict.items():
                        if value_raw is None:
                            continue
                        key = str(key_raw)
                        new_value = str(value_raw)
                        old_value = str(current.get(key, ""))
                        if new_value and new_value != old_value:
                            year_construction_events.setdefault(archetype_name, {})[key] = (old_value, new_value)
                            current[key] = new_value
                    continue

                if component_name == "supply_systems":
                    continue
                if component_name not in LAYERED_COMPONENT_TO_DB:
                    raise ValueError(
                        f"Unsupported modified component '{component_name}' in district pathway log. "
                        f"Supported components: {sorted(LAYERED_COMPONENT_TO_DB.keys())}"
                    )
                old_layers = archetype_layers[archetype_name].get(component_name)
                new_layers = _apply_layer_patch(old_layers, patch or {})
                events = _diff_layers(old_layers, new_layers)
                if events:
                    year_layer_events.setdefault(archetype_name, {})[component_name] = events
                archetype_layers[archetype_name][component_name] = new_layers

        layer_events_by_year[year] = year_layer_events
        layers_at_year[year] = {
            archetype: {component: layers[:] for component, layers in components.items()}
            for archetype, components in archetype_layers.items()
        }
        construction_events_by_year[year] = year_construction_events
        construction_at_year[year] = {
            archetype: dict(codes)
            for archetype, codes in archetype_construction_types.items()
        }

    return ArchetypeChangesTimeline(
        years_sorted=list(effective_years),
        layers_at_year=layers_at_year,
        layer_events_by_year=layer_events_by_year,
        construction_types_at_year=construction_at_year,
        construction_type_events_by_year=construction_events_by_year,
    )


def _apply_layer_patch(
    old_layers: list[MaterialLayer] | None,
    patch: dict[str, Any],
) -> list[MaterialLayer]:
    base_layers = _coerce_3_layers(old_layers)
    slots: list[tuple[str, float]] = [(layer.name, float(layer.thickness_m)) for layer in base_layers]

    for idx in (1, 2, 3):
        material_key = f"material_name_{idx}"
        thickness_key = f"thickness_{idx}_m"
        name, thickness = slots[idx - 1]
        if material_key in patch:
            raw = patch.get(material_key)
            if raw is not None and not (isinstance(raw, float) and np.isnan(raw)):
                name = str(raw).strip()
        if thickness_key in patch:
            raw = patch.get(thickness_key)
            if raw is not None and not (isinstance(raw, float) and np.isnan(raw)):
                thickness = float(raw)
        if thickness < 0.0:
            thickness = 0.0
        if not str(name).strip():
            name = ""
            thickness = 0.0
        slots[idx - 1] = (name, float(thickness))

    return [MaterialLayer(name=name, thickness_m=thickness) for name, thickness in slots]


def _parse_density(value: Any) -> float | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text or text == "-":
        return None
    if "-" in text and all(part.strip().replace(".", "", 1).isdigit() for part in text.split("-")):
        lower, upper = text.split("-", 1)
        return 0.5 * (float(lower) + float(upper))
    try:
        return float(text)
    except Exception:
        return None


def _to_float(value: Any, default: float = 0.0) -> float:
    if value is None:
        return float(default)
    if isinstance(value, float) and np.isnan(value):
        return float(default)
    if isinstance(value, pd.Series):
        if len(value) == 0:
            return float(default)
        return _to_float(value.iloc[0], default=default)
    try:
        return float(value)
    except Exception:
        return float(default)


def read_material_db(locator: InputLocator) -> pd.DataFrame:
    path = locator.get_database_components_materials()
    if not os.path.exists(path):
        raise FileNotFoundError(f"Materials database not found: {path}")
    df = pd.read_csv(path)
    if "name" not in df.columns:
        raise ValueError("Materials database missing required 'name' column")
    return df.set_index("name", drop=False)


def _layers_from_envelope_row(row: pd.Series) -> list[MaterialLayer]:
    layers: list[MaterialLayer] = []
    for idx in (1, 2, 3):
        material = row.get(f"material_name_{idx}")
        thickness = row.get(f"thickness_{idx}_m")

        name = ""
        if material is not None and not (isinstance(material, float) and np.isnan(material)):
            name = str(material).strip()

        thickness_m = 0.0
        if thickness is not None and not (isinstance(thickness, float) and np.isnan(thickness)):
            try:
                thickness_m = float(thickness)
            except Exception:
                thickness_m = 0.0
        if thickness_m < 0.0:
            thickness_m = 0.0
        if not name:
            thickness_m = 0.0

        layers.append(MaterialLayer(name=name, thickness_m=thickness_m))
    return layers


def get_component_layers(
    env: EnvelopeLookup,
    *,
    db_name: str,
    code: str,
) -> list[MaterialLayer]:
    df = getattr(env.envelope, db_name)
    if df is None:
        return _empty_layers()
    if code not in df.index:
        raise ValueError(f"Envelope code '{code}' not found in DB '{db_name}'.")
    row = df.loc[code]
    required = {"material_name_1", "thickness_1_m"}
    if not required.issubset(set(df.columns)):
        raise ValueError(
            f"Envelope DB '{db_name}' is missing required layer columns {sorted(required)}."
        )
    return _layers_from_envelope_row(row)


def _material_intensity_per_m2(
    materials: pd.DataFrame,
    layer: MaterialLayer,
) -> tuple[float, float, float]:
    """Resolve one raw layer into per-m2 embodied intensities.

    This is pathway-specific translation logic from `MaterialLayer` to the
    backend's `ResolvedComponent` fields. It reads `MATERIALS.csv` directly
    because pathway works on individual material layers, while the current
    assemblies helpers derive legacy per-assembly columns in-memory and do not
    expose a reusable public per-layer API.

    The returned tuple is:
    - production [database-native intensity times kg/m2]
    - demolition / recycling [database-native intensity times kg/m2]
    - biogenic [database-native intensity times kg/m2]

    Note:
    `biogenic_carbon_in_product` is forwarded as stored in the materials
    database. If a database provides this field in kg C/kg rather than a
    kgCO2e-compatible intensity, conversion must happen before or inside this
    helper; the backend itself only applies the sign convention.
    """
    if layer.name not in materials.index:
        raise ValueError(f"Material '{layer.name}' not found in MATERIALS.csv")
    rec = materials.loc[layer.name]
    unit = str(rec.get("unit", "")).strip().lower()
    if unit != "kg":
        raise ValueError(
            f"Material '{layer.name}' has unit '{rec.get('unit')}', expected 'kg' for layered constructions."
        )

    density = _parse_density(rec.get("density"))
    if density is None or density <= 0:
        raise ValueError(f"Material '{layer.name}' has invalid density: {rec.get('density')}")

    mass_per_m2 = density * float(layer.thickness_m)
    production = _to_float(rec.get("GHG_emission_production"), default=0.0) * mass_per_m2
    demolition = _to_float(rec.get("GHG_emission_recycling"), default=0.0) * mass_per_m2
    biogenic = _to_float(rec.get("biogenic_carbon_in_product"), default=0.0) * mass_per_m2
    return production, demolition, biogenic


def _diff_layers(
    old_layers: list[MaterialLayer] | None,
    new_layers: list[MaterialLayer] | None,
) -> list[tuple[str, MaterialLayer]]:
    old3 = _coerce_3_layers(old_layers)
    new3 = _coerce_3_layers(new_layers)
    events: list[tuple[str, MaterialLayer]] = []
    for idx in range(3):
        old_layer = old3[idx]
        new_layer = new3[idx]
        if old_layer.name != new_layer.name or abs(old_layer.thickness_m - new_layer.thickness_m) > 1e-9:
            events.append(("remove", old_layer))
            events.append(("add", new_layer))
    return events


def _load_building_const_types(locator: InputLocator) -> dict[str, str]:
    zone = gpd.read_file(locator.get_zone_geometry())
    if "name" not in zone.columns:
        raise ValueError("Zone geometry is missing required 'name' column.")
    if "const_type" not in zone.columns:
        raise ValueError("Zone geometry is missing required 'const_type' column.")
    out: dict[str, str] = {}
    for _, row in zone.iterrows():
        name = str(row["name"])
        const_type = row["const_type"]
        if const_type is None or (isinstance(const_type, float) and np.isnan(const_type)):
            continue
        out[name] = str(const_type)
    return out


def _building_demolition_years(log_data: dict[int, dict[str, Any]]) -> dict[str, int]:
    out: dict[str, int] = {}
    for year in sorted(int(value) for value in log_data.keys()):
        entry = log_data.get(year, {}) or {}
        events = entry.get("building_events", {}) or {}
        demolished = events.get("demolished_buildings", []) or []
        for building in demolished:
            out.setdefault(str(building), int(year))
    return out


class PathwayTimelineFrontend(TimelineFrontendBase):
    """Archetype-authored, building-applied frontend for pathway timelines.

    This frontend owns pathway-specific state reconstruction only:

    - interpret archetype/log snapshots and year-local events
    - resolve the effective code or layer stack at a given year
    - translate those states into backend-facing `LifecycleChange` objects

    Effective embodied state is cached per construction type via
    `EffectiveBuildingStateCache`, so the frontend can ask for
    `effective_state_cache.at(year)` instead of repeating many independent
    code/layer snapshot lookups.
    """

    def __init__(
        self,
        *,
        building_name: str,
        locator: InputLocator,
        const_type: str,
        area_dict: dict[str, float],
        materials: pd.DataFrame,
        years_sorted: list[int],
        start_year: int,
        end_year: int,
        construction_year: int,
        demolition_year: int | None,
        archetype_timeline: ArchetypeChangesTimeline,
        service_life_by_src_component: Mapping[str, int | None],
        operational_by_state_year: Mapping[int, pd.DataFrame | None],
        allow_missing_operational: bool,
    ) -> None:
        self.name = str(building_name)
        self.locator = locator
        self.const_type = str(const_type)
        self.area_dict = area_dict
        self.materials = materials
        self.years_sorted = list(years_sorted)
        self._start_year = int(start_year)
        self._end_year = int(end_year)
        self.construction_year = int(construction_year)
        self.demolition_year = int(demolition_year) if demolition_year is not None else None
        self.archetype_timeline = archetype_timeline
        self.service_life_by_src_component = service_life_by_src_component
        self.operational_by_state_year = operational_by_state_year
        self.allow_missing_operational = bool(allow_missing_operational)
        self.envelope_lookup = EnvelopeLookup.from_locator(locator)
        self.emission_per_tech = float(EMISSIONS_EMBODIED_TECHNICAL_SYSTEMS) / len(TECHNICAL_SYSTEM_COMPONENTS)
        win_lifetime_any = self.service_life_by_src_component.get("win")
        if win_lifetime_any is None:
            raise ValueError(
                f"Missing Service_Life for component 'win' (const_type={self.const_type}). "
                "Service life is mandatory for window replacement scheduling."
            )
        self.win_lifetime_default = int(win_lifetime_any)
        if self.win_lifetime_default <= 0:
            raise ValueError(
                f"Invalid Service_Life={self.win_lifetime_default} for component 'win' (const_type={self.const_type})."
            )
        tech_lifetime_any = self.service_life_by_src_component.get("technical_systems")
        self.tech_lifetime = int(tech_lifetime_any) if tech_lifetime_any is not None else int(SERVICE_LIFE_OF_TECHNICAL_SYSTEMS)
        self.effective_state_cache = EffectiveBuildingStateCache(
            const_type=self.const_type,
            construction_year=self.construction_year,
            years_sorted=self.years_sorted,
            archetype_timeline=self.archetype_timeline,
            envelope_lookup=self.envelope_lookup,
            win_lifetime_default=self.win_lifetime_default,
        )
        self.has_embodied_state = self.effective_state_cache.has_embodied_state

    def building_name(self) -> str:
        return self.name

    def start_year(self) -> int:
        return self._start_year

    def end_year(self) -> int:
        return self._end_year

    def initial_change(self) -> LifecycleChange | None:
        if not self.has_embodied_state:
            return LifecycleChange(
                year=self.construction_year,
                mode="add",
                components=(),
                note="Constructed | No envelope layer snapshot available for this construction type; skipping embodied emissions.",
            )

        return LifecycleChange(
            year=self.construction_year,
            mode="add",
            components=tuple(self._initial_components(self.construction_year)),
            note="Constructed",
        )

    def initial_due_years(self) -> dict[str, int]:
        if not self.has_embodied_state:
            return {}

        areas_by_src_component: dict[str, float] = {}
        for _, src_component, area_key in COMP_AREA_MAP:
            areas_by_src_component[src_component] = areas_by_src_component.get(src_component, 0.0) + float(
                self.area_dict.get(area_key, 0.0)
            )

        due_years: dict[str, int] = {}
        for src_component, total_area in areas_by_src_component.items():
            if total_area <= 0.0 or src_component in ("technical_systems", "win"):
                continue
            due_years[src_component] = self.construction_year + self._service_life_for_src_component(src_component)

        if float(self.area_dict.get("Awin_ag", 0.0)) > 0.0:
            due_years["win"] = self.construction_year + self.effective_state_cache.at(self.construction_year).window_lifetime

        if float(self.area_dict.get("Atechnical_systems", 0.0)) > 0.0:
            for component in TECHNICAL_SYSTEM_COMPONENTS:
                due_years[component] = self.construction_year + self.tech_lifetime

        return due_years

    def authored_changes_by_year(self) -> dict[int, tuple[LifecycleChange, ...]]:
        if not self.has_embodied_state:
            return {}

        changes_by_year: dict[int, tuple[LifecycleChange, ...]] = {}
        for year in self.years_sorted:
            if not self._exists_at(year):
                continue

            year_changes: list[LifecycleChange] = []
            year_mods = self.archetype_timeline.events_for_year(year).get(self.const_type, {}) or {}
            for src_component, events in year_mods.items():
                detail_parts = [f"-{layer.name} {layer.thickness_m:.3f}m" for action, layer in events if action == "remove"]
                detail_parts.extend(f"+{layer.name} {layer.thickness_m:.3f}m" for action, layer in events if action == "add")
                note = f"Modified {src_component}" + (f": {', '.join(detail_parts)}" if detail_parts else "")
                year_changes.extend(self._modification_changes(year=year, src_component=src_component, events=events, note=note))

            patch = self.archetype_timeline.construction_type_events_for_year(year).get(self.const_type, {}) or {}
            if CONSTRUCTION_TYPE_WIN_FIELD in patch:
                change = _parse_code_change(patch.get(CONSTRUCTION_TYPE_WIN_FIELD))
                if change is not None:
                    old_code, new_code = change
                    if new_code and new_code != old_code:
                        area_win = float(self.area_dict.get("Awin_ag", 0.0))
                        components: tuple[ResolvedComponent, ...] = ()
                        if area_win > 0.0:
                            production_new, _, biogenic_new = envelope_intensities_per_m2(
                                self.envelope_lookup,
                                code=new_code,
                            )
                            _, demolition_old, _ = envelope_intensities_per_m2(
                                self.envelope_lookup,
                                code=old_code,
                            )
                            components = (
                                ResolvedComponent(
                                    component="win_ag",
                                    area_m2=area_win,
                                    production_per_area=production_new,
                                    demolition_per_area=demolition_old,
                                    biogenic_per_area=biogenic_new,
                                ),
                            )
                        year_changes.append(
                            LifecycleChange(
                                year=year,
                                mode="replace",
                                components=components,
                                note=f"Window code changed: {old_code} -> {new_code}",
                                reset_keys=("win",),
                            )
                        )

            area_tech = float(self.area_dict.get("Atechnical_systems", 0.0))
            for supply_field, sys_component in TECH_SYSTEM_COMPONENT_BY_SUPPLY_FIELD.items():
                if supply_field not in patch:
                    continue
                change = _parse_code_change(patch.get(supply_field))
                if change is None:
                    continue
                old_code, new_code = change
                if not new_code or new_code == old_code:
                    continue
                label = supply_field.replace("supply_type_", "")
                year_changes.append(
                    LifecycleChange(
                        year=year,
                        mode="production_only",
                        components=tuple(self._resolved_technical_components([sys_component], area_m2=area_tech)),
                        note=f"Supply system changed ({label}): {old_code} -> {new_code}",
                        reset_keys=(sys_component,),
                    )
                )

            if year_changes:
                changes_by_year[int(year)] = tuple(year_changes)

        return changes_by_year

    def replacement_change(self, key: str, year: int) -> LifecycleChange | None:
        if not self.has_embodied_state or not self._exists_at(year):
            return None

        if key in TECHNICAL_SYSTEM_COMPONENTS:
            return LifecycleChange(
                year=year,
                mode="production_only",
                components=tuple(self._resolved_technical_components([key], area_m2=float(self.area_dict.get("Atechnical_systems", 0.0)))),
            )

        if key == "win":
            return LifecycleChange(
                year=year,
                mode="replace",
                components=tuple(self._full_replacement_components(year=year, src_component="win")),
                note="Service life reached: windows",
            )

        state = self.effective_state_cache.at(year)
        layers = state.layers_by_src_component.get(key, _empty_layers())
        layer_desc = ", ".join(f"{layer.name or '-'} {layer.thickness_m:.3f}m" for layer in layers)
        return LifecycleChange(
            year=year,
            mode="replace",
            components=tuple(self._full_replacement_components(year=year, src_component=key)),
            note=f"Service life reached: {key} (replace {layer_desc})",
        )

    def next_due_year(self, key: str, year: int) -> int | None:
        if not self.has_embodied_state:
            return None
        if key in TECHNICAL_SYSTEM_COMPONENTS:
            return int(year) + int(self.tech_lifetime)
        if key == "win":
            return int(year) + int(self.effective_state_cache.at(year).window_lifetime)
        if key in MATERIAL_SRC_COMPONENTS | {"win"}:
            return int(year) + int(self._service_life_for_src_component(key))
        return None

    def operational_segments(self) -> tuple[OperationalSegment, ...]:
        active_end_year = self._end_year if self.demolition_year is None else min(self._end_year, self.demolition_year - 1)
        if active_end_year < self.construction_year:
            return ()

        segments: list[OperationalSegment] = []
        last_known: pd.Series | None = None
        for idx, state_year in enumerate(self.years_sorted):
            if state_year > self._end_year:
                break
            df_state = self.operational_by_state_year.get(state_year)
            if df_state is not None and self.name in df_state.index:
                op_state_year = self._operational_row(df_state)
                last_known = op_state_year
            elif self.allow_missing_operational and last_known is not None:
                op_state_year = last_known
            else:
                op_state_year = pd.Series(dtype=float)

            next_state_year = self.years_sorted[idx + 1] if idx + 1 < len(self.years_sorted) else (self._end_year + 1)
            interval_start = max(int(state_year), self.construction_year)
            interval_end = min(int(next_state_year) - 1, active_end_year)
            if interval_start > interval_end or op_state_year.empty:
                continue

            segments.append(
                OperationalSegment(
                    start_year=interval_start,
                    end_year=interval_end,
                    column_values={str(column): float(value) for column, value in op_state_year.items()},
                )
            )

        return tuple(segments)

    def demolition_change(self) -> LifecycleChange | None:
        if not self.has_embodied_state or self.demolition_year is None:
            return None
        if not (self._start_year <= self.demolition_year <= self._end_year):
            return None
        return LifecycleChange(
            year=self.demolition_year,
            mode="demolish",
            components=tuple(self._demolition_components(self.demolition_year)),
            note="Demolished",
        )

    def _exists_at(self, year: int) -> bool:
        return int(year) >= self.construction_year and (
            self.demolition_year is None or int(year) < self.demolition_year
        )

    def _service_life_for_src_component(self, src_component: str) -> int:
        lifetime_any = self.service_life_by_src_component.get(src_component)
        if lifetime_any is None:
            raise ValueError(
                f"Missing Service_Life for component '{src_component}' (const_type={self.const_type}). "
                "Service life is mandatory for the district pathway emissions timeline replacement scheduling."
            )
        lifetime = int(lifetime_any)
        if lifetime <= 0:
            raise ValueError(
                f"Invalid Service_Life={lifetime} for component '{src_component}' (const_type={self.const_type})."
            )
        return lifetime

    def _initial_components(self, year: int) -> list[ResolvedComponent]:
        state = self.effective_state_cache.at(year)
        layers_snapshot = state.layers_by_src_component
        components: list[ResolvedComponent] = []
        window_code = state.window_code
        for component, src_component, area_key in COMP_AREA_MAP:
            area = float(self.area_dict.get(area_key, 0.0))
            if area <= 0.0:
                continue
            if src_component == "technical_systems":
                components.extend(self._resolved_technical_components(TECHNICAL_SYSTEM_COMPONENTS, area_m2=area))
                continue
            if src_component == "win":
                production, _, biogenic = envelope_intensities_per_m2(self.envelope_lookup, code=window_code)
                components.append(
                    ResolvedComponent(
                        component=component,
                        area_m2=area,
                        production_per_area=production,
                        biogenic_per_area=biogenic,
                    )
                )
                continue

            layers = layers_snapshot.get(src_component, _empty_layers())
            for action, layer in _diff_layers(_empty_layers(), layers):
                if action == "remove" or not _is_layer_active(layer):
                    continue
                production, _, biogenic = _material_intensity_per_m2(self.materials, layer)
                components.append(
                    ResolvedComponent(
                        component=component,
                        area_m2=area,
                        production_per_area=production,
                        biogenic_per_area=biogenic,
                    )
                )
        return components

    def _modification_changes(
        self,
        *,
        year: int,
        src_component: str,
        events: list[tuple[str, MaterialLayer]],
        note: str,
    ) -> list[LifecycleChange]:
        additions: list[ResolvedComponent] = []
        demolitions: list[ResolvedComponent] = []
        for component, comp_src, area_key in COMP_AREA_MAP:
            if comp_src != src_component:
                continue
            area = float(self.area_dict.get(area_key, 0.0))
            if area <= 0.0 or comp_src == "technical_systems":
                continue
            for action, layer in events:
                if not layer.name:
                    continue
                production, demolition, biogenic = _material_intensity_per_m2(self.materials, layer)
                if action == "add":
                    additions.append(
                        ResolvedComponent(
                            component=component,
                            area_m2=area,
                            production_per_area=production,
                            biogenic_per_area=biogenic,
                        )
                    )
                else:
                    demolitions.append(
                        ResolvedComponent(
                            component=component,
                            area_m2=area,
                            demolition_per_area=demolition,
                        )
                    )

        changes: list[LifecycleChange] = []
        if demolitions:
            changes.append(
                LifecycleChange(
                    year=year,
                    mode="demolish",
                    components=tuple(demolitions),
                    note=note,
                    reset_keys=(src_component,),
                )
            )
            note = ""
        if additions:
            changes.append(
                LifecycleChange(
                    year=year,
                    mode="add",
                    components=tuple(additions),
                    note=note or None,
                    reset_keys=(src_component,),
                )
            )
        return changes

    def _full_replacement_components(self, *, year: int, src_component: str) -> list[ResolvedComponent]:
        components: list[ResolvedComponent] = []
        if src_component == "technical_systems":
            return components
        state = self.effective_state_cache.at(year)

        for component, comp_src, area_key in COMP_AREA_MAP:
            if comp_src != src_component:
                continue
            area = float(self.area_dict.get(area_key, 0.0))
            if area <= 0.0:
                continue

            if comp_src == "win":
                production, demolition, biogenic = envelope_intensities_per_m2(
                    self.envelope_lookup,
                    code=state.window_code,
                )
                components.append(
                    ResolvedComponent(
                        component=component,
                        area_m2=area,
                        production_per_area=production,
                        demolition_per_area=demolition,
                        biogenic_per_area=biogenic,
                    )
                )
                continue

            layers = state.layers_by_src_component.get(src_component, _empty_layers())
            for layer in _active_layers(layers):
                production, demolition, biogenic = _material_intensity_per_m2(self.materials, layer)
                components.append(
                    ResolvedComponent(
                        component=component,
                        area_m2=area,
                        production_per_area=production,
                        demolition_per_area=demolition,
                        biogenic_per_area=biogenic,
                    )
                )
        return components

    def _demolition_components(self, year: int) -> list[ResolvedComponent]:
        components: list[ResolvedComponent] = []
        state = self.effective_state_cache.at(year)
        layers_snapshot = state.layers_by_src_component
        for component, src_component, area_key in COMP_AREA_MAP:
            area = float(self.area_dict.get(area_key, 0.0))
            if area <= 0.0 or src_component == "technical_systems":
                continue

            if src_component == "win":
                _, demolition, _ = envelope_intensities_per_m2(
                    self.envelope_lookup,
                    code=state.window_code,
                )
                components.append(
                    ResolvedComponent(
                        component=component,
                        area_m2=area,
                        demolition_per_area=demolition,
                    )
                )
                continue

            layers = layers_snapshot.get(src_component, _empty_layers())
            for action, layer in _diff_layers(_empty_layers(), layers):
                if action == "remove" or not _is_layer_active(layer):
                    continue
                _, demolition, _ = _material_intensity_per_m2(self.materials, layer)
                components.append(
                    ResolvedComponent(
                        component=component,
                        area_m2=area,
                        demolition_per_area=demolition,
                    )
                )
        return components

    def _resolved_technical_components(
        self,
        components: Sequence[str],
        *,
        area_m2: float,
    ) -> list[ResolvedComponent]:
        if area_m2 <= 0.0:
            return []
        return [
            ResolvedComponent(
                component=str(component),
                area_m2=area_m2,
                production_per_area=self.emission_per_tech,
            )
            for component in components
        ]

    def _operational_row(self, df_state: pd.DataFrame) -> pd.Series:
        cols = [
            column
            for column in df_state.columns
            if isinstance(column, str) and column.endswith("_kgCO2e")
        ]
        if not cols:
            return pd.Series(dtype=float)
        selected_row = df_state.loc[self.name]
        if isinstance(selected_row, pd.DataFrame):
            if selected_row.empty:
                return pd.Series(dtype=float)
            selected_row = selected_row.iloc[0]
        return selected_row.reindex(cols).astype(float)


def create_district_pathway_emissions_timeline(
    config: Configuration,
    *,
    pathway_name: str,
    allow_missing_operational: bool = False,
) -> pd.DataFrame:
    """Build the district pathway yearly emissions timeline.

    High-level flow in 7 steps:

    1. Resolve required pathway years and load / normalise the pathway log.
    2. Resolve the shared timeline horizon and backend configuration.
    3. Validate that each state year has the radiation and operational inputs
       needed to derive areas and operational emissions.
    4. Collect per-building construction metadata, component areas, and
       per-state operational totals.
    5. Load archetype, envelope, and materials databases, then reconstruct the
       cumulative archetype snapshot timeline used by the pathway frontend.
    6. Instantiate one `PathwayTimelineFrontend` per eligible building and let
       the shared backend assemble each yearly building timeline.
    7. Finalise, save, and aggregate per-building timelines into the district
       pathway timeline dataframe returned by this function.

    The function is intentionally orchestration-heavy: file I/O, pathway-state
    validation, archetype reconstruction, per-building frontend creation, and
    CSV persistence all happen here, while the shared backend still owns the
    generic timeline assembly rules.
    """
    main_locator = InputLocator(config.scenario)
    years = get_required_state_years(config, pathway_name)
    log_data = ensure_state_years_exist(
        config,
        years=years,
        update_yaml=True,
        pathway_name=pathway_name,
    )

    start_year = min(years)
    end_year = resolve_emissions_end_year(config.emissions)
    backend = EmissionTimelineBackend(
        feedstock_policies=feedstock_policies_from_config(config.emissions),
    )

    base_cols = emission_timeline_columns(
        include_legacy_aliases=True,
        include_phase_totals=True,
    ) + [f"operation_{demand}_kgCO2e" for demand in ("Qhs_sys", "Qww_sys", "Qcs_sys", "E_sys")]

    def _empty_timeline() -> pd.DataFrame:
        idx = [f"Y_{year}" for year in range(start_year, end_year + 1)]
        df = pd.DataFrame(0.0, index=idx, columns=base_cols)
        df.index.name = "period"
        df["Note"] = ""
        return df

    years_sorted = sorted(years)
    missing_radiation: dict[int, str] = {}
    missing_operational_years: list[int] = []
    for year in years_sorted:
        state_locator = InputLocator(main_locator.get_state_in_time_scenario_folder(pathway_name, year))
        buildings = list(state_locator.get_zone_building_names())
        if buildings:
            for building in buildings:
                rad_path = state_locator.get_radiation_building(building)
                if not os.path.exists(rad_path):
                    missing_radiation.setdefault(year, rad_path)
                    break
        op_path = state_locator.get_total_yearly_operational_building()
        if not os.path.exists(op_path):
            missing_operational_years.append(year)

    if missing_radiation:
        sample_lines = "\n".join(f"- {year}: missing {path}" for year, path in sorted(missing_radiation.items()))
        raise FileNotFoundError(
            "Some state years are missing solar-radiation outputs required for surface areas.\n"
            "Run the `pathway-simulations` script (includes Radiation) for all state years, then rerun this pathway emissions timeline.\n"
            "Missing examples:\n" + sample_lines
        )

    if missing_operational_years and not allow_missing_operational:
        years_str = ", ".join(map(str, sorted(set(missing_operational_years))))
        raise FileNotFoundError(
            "Missing operational-by-building results for some state years: "
            f"{years_str}.\n"
            "Run the `pathway-simulations` script (includes Emissions) to generate outputs for all state years, "
            "or call create_district_pathway_emissions_timeline(..., allow_missing_operational=True) for best-effort carry-forward."
        )

    building_construction_years = get_building_construction_years(main_locator)
    building_const_types = _load_building_const_types(main_locator)
    demolition_years = _building_demolition_years(log_data)
    areas_by_building: dict[str, dict[str, float]] = {}
    operational_by_state_year: dict[int, pd.DataFrame] = {}

    for year in years_sorted:
        state_locator = InputLocator(main_locator.get_state_in_time_scenario_folder(pathway_name, year))
        buildings_in_state = list(state_locator.get_zone_building_names())
        if not buildings_in_state:
            continue

        op_path = state_locator.get_total_yearly_operational_building()
        if os.path.exists(op_path):
            operational_by_state_year[year] = pd.read_csv(op_path, index_col="name")

        missing_area_buildings = [building for building in buildings_in_state if building not in areas_by_building]
        if missing_area_buildings:
            building_properties = load_building_properties_for_emissions(state_locator, missing_area_buildings)
            for building in missing_area_buildings:
                areas_by_building[building] = get_component_quantities(building_properties, building)

    archetype_df = pd.read_csv(main_locator.get_database_archetypes_construction_type(), index_col="const_type")
    env_lookup = EnvelopeLookup.from_locator(main_locator)
    materials = read_material_db(main_locator)

    archetypes_needed = sorted({value for value in building_const_types.values()})
    archetype_layers: dict[str, dict[str, list[MaterialLayer]]] = {}
    archetype_service_life: dict[str, dict[str, int | None]] = {}
    archetype_construction_types: dict[str, dict[str, str]] = {}
    for archetype in archetypes_needed:
        if archetype not in archetype_df.index:
            raise ValueError(f"Archetype '{archetype}' not found in construction types database.")
        row = archetype_df.loc[archetype]
        layered_codes = {component: str(row.get(f"type_{component}")) for component in LAYERED_COMPONENT_TO_DB}
        code_win = str(row.get("type_win"))
        archetype_construction_types[archetype] = {
            CONSTRUCTION_TYPE_WIN_FIELD: code_win,
            **{field: str(row.get(field)) for field in SUPPLY_TYPE_FIELDS},
        }
        archetype_layers[archetype] = {
            component: get_component_layers(
                env_lookup,
                db_name=LAYERED_COMPONENT_TO_DB[component],
                code=layered_codes[component],
            )
            for component in LAYERED_COMPONENT_TO_DB
        }
        archetype_service_life[archetype] = {
            **{
                component: cast(int | None, env_lookup.get_item_value(layered_codes[component], "Service_Life"))
                for component in LAYERED_COMPONENT_TO_DB
            },
            "win": cast(int | None, env_lookup.get_item_value(code_win, "Service_Life")),
            "technical_systems": int(SERVICE_LIFE_OF_TECHNICAL_SYSTEMS),
        }

    in_range_construction_years = [
        int(year)
        for year in building_construction_years.values()
        if year is not None and start_year <= int(year) <= end_year
    ]
    baseline_year = min([start_year, *in_range_construction_years]) if in_range_construction_years else start_year
    archetype_timeline = prepare_pathway_archetype_timeline(
        years_sorted=years_sorted,
        base_year=int(baseline_year),
        log_data=log_data,
        archetype_layers=archetype_layers,
        archetype_construction_types=archetype_construction_types,
    )

    building_timelines: dict[str, pd.DataFrame] = {}
    all_buildings = sorted(set(building_const_types.keys()) | set(building_construction_years.keys()))
    for building in all_buildings:
        construction_year = building_construction_years.get(building)
        if construction_year is None or construction_year < start_year or construction_year > end_year:
            continue

        const_type = building_const_types.get(building)
        area_dict = areas_by_building.get(building)
        if const_type is None or area_dict is None:
            continue

        frontend = PathwayTimelineFrontend(
            building_name=building,
            locator=main_locator,
            const_type=const_type,
            area_dict=area_dict,
            materials=materials,
            years_sorted=years_sorted,
            start_year=start_year,
            end_year=end_year,
            construction_year=int(construction_year),
            demolition_year=demolition_years.get(building),
            archetype_timeline=archetype_timeline,
            service_life_by_src_component=archetype_service_life[const_type],
            operational_by_state_year=operational_by_state_year,
            allow_missing_operational=allow_missing_operational,
        )
        building_timelines[building] = backend.build(frontend)

    per_building_folder = main_locator.get_building_pathway_emissions_timelines_folder(pathway_name)
    os.makedirs(per_building_folder, exist_ok=True)
    for building, timeline in building_timelines.items():
        save_path = main_locator.get_building_pathway_emissions_timeline_file(pathway_name, building_name=building)
        try:
            df_save = timeline.copy()
            df_save["name"] = building
            df_save.to_csv(save_path, float_format="%.2f")
        except PermissionError as exc:
            raise PermissionError(
                "Permission denied writing a per-building pathway emissions timeline CSV. "
                "This often happens on Windows when the CSV is open in Excel or locked by OneDrive sync. "
                f"Close '{save_path}' and rerun."
            ) from exc

    if building_timelines:
        out = backend.aggregate_by_index(list(building_timelines.values())).set_index("period")
        out.index.name = "period"
    else:
        out = _empty_timeline()

    out = finalise_emission_timeline_dataframe(out)
    save_path = main_locator.get_district_pathway_emissions_timeline_path(pathway_name)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    try:
        out.to_csv(save_path, float_format="%.2f")
    except PermissionError as exc:
        raise PermissionError(
            "Permission denied writing the district pathway emissions timeline output. "
            "This often happens on Windows when the CSV is open in Excel or locked by OneDrive sync. "
            f"Close '{save_path}' and rerun."
        ) from exc

    return out


def main(config: Configuration) -> None:
    pathway_name = config.pathway_simulations.existing_pathway_name
    if not pathway_name:
        raise ValueError(
            "No existing pathway name provided. "
            "Please provide an existing pathway name to create the district pathway emissions timeline."
        )
    df = create_district_pathway_emissions_timeline(config, pathway_name=pathway_name)
    print(f"District pathway emissions timeline saved with {len(df)} years.")


if __name__ == "__main__":
    main(Configuration())
