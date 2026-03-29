"""Shared backend for yearly emission timelines.

This module defines the canonical timeline contract used by both workflows:

- frontends resolve domain-specific state into generic lifecycle/operational inputs
- the backend turns those inputs into timeline rows, totals, and aliases

In practice this means:
- `ResolvedComponent` describes one fully resolved embodied component contribution
- `LifecycleChange` describes one event in one year
- `OperationalSegment` describes one constant yearly operational interval
- `TimelineFrontendBase` is the adapter contract for normal and pathway frontends
- `EmissionTimelineBackend` owns the generic timeline assembly logic
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, TypeAlias, cast

import numpy as np
import pandas as pd

from cea.analysis.lca.hourly_operational_emission import _tech_name_mapping

TimelineYear: TypeAlias = int | list[int] | str | list[str]

EMISSION_PHASES: tuple[str, ...] = ("production", "demolition", "biogenic")
TECHNICAL_SYSTEM_COMPONENTS: tuple[str, ...] = (
    "technical_system_hs",
    "technical_system_cs",
    "technical_system_dhw",
    "technical_system_el",
)
COMPONENT_TO_SRC_COMPONENT: dict[str, str] = {
    "wall_ag": "wall",
    "wall_bg": "base",
    "wall_part": "part",
    "win_ag": "win",
    "roof": "roof",
    "upperside": "roof",
    "underside": "base",
    "floor": "floor",
    "base": "base",
    "technical_systems": "technical_systems",
}
TIMELINE_COMPONENTS: tuple[str, ...] = tuple(COMPONENT_TO_SRC_COMPONENT.keys())
EMISSION_COMPONENTS: tuple[str, ...] = tuple(
    component for component in TIMELINE_COMPONENTS if component != "technical_systems"
) + TECHNICAL_SYSTEM_COMPONENTS
PHASE_TOTAL_COLUMNS: tuple[str, ...] = tuple(
    f"{phase}_kgCO2e" for phase in EMISSION_PHASES
)
OPERATION_TOTAL_COLUMNS: tuple[str, ...] = tuple(
    f"operation_{demand}_kgCO2e" for demand in _tech_name_mapping.keys()
)


def normalise_timeline_years(year: TimelineYear) -> list[str]:
    if isinstance(year, int):
        return [f"Y_{int(year)}"]
    if isinstance(year, str):
        return [year]
    if isinstance(year, list):
        if not year:
            return []
        if isinstance(year[0], int):
            return [f"Y_{int(value)}" for value in cast(list[int], year)]
        if isinstance(year[0], str):
            return list(cast(list[str], year))
    raise ValueError(f"Year must be int, str, list[int], or list[str]; got {type(year)}.")


def years_from_timeline_index(index: pd.Index) -> list[int]:
    years: list[int] = []
    for label in index:
        if isinstance(label, int):
            years.append(int(label))
            continue
        text = str(label).strip()
        if text.startswith("Y_"):
            text = text[2:]
        if not text or not text.lstrip("-").isdigit():
            raise ValueError(
                "Index must contain years only (int years or 'Y_YYYY' labels). "
                f"Got invalid label: {label!r}"
            )
        years.append(int(text))
    return years


def resolve_emissions_end_year(emissions_cfg: Any, default_year: int = 2100) -> int:
    raw = getattr(emissions_cfg, "year_end", None)
    if raw in (None, ""):
        return int(default_year)
    return int(raw)


def feedstock_policies_from_config(
    emissions_cfg: Any,
) -> dict[str, tuple[int, int, float]] | None:
    ref_year = getattr(emissions_cfg, "grid_decarbonise_reference_year", None)
    target_year = getattr(emissions_cfg, "grid_decarbonise_target_year", None)
    target_factor = getattr(emissions_cfg, "grid_decarbonise_target_emission_factor", None)

    if ref_year is not None and target_year is not None and target_factor is not None:
        return {"GRID": (int(ref_year), int(target_year), float(target_factor))}
    if ref_year is None and target_year is None and target_factor is None:
        return None
    raise ValueError(
        "If one of grid_decarbonise_reference_year, grid_decarbonise_target_year, "
        "or grid_decarbonise_target_emission_factor is set, all must be set."
    )


def emission_timeline_columns(
    *,
    include_legacy_aliases: bool = False,
    include_phase_totals: bool = False,
) -> list[str]:
    columns = [
        f"{phase}_{component}_kgCO2e"
        for phase in EMISSION_PHASES
        for component in EMISSION_COMPONENTS
    ]
    if include_legacy_aliases:
        columns.extend(f"{phase}_technical_systems_kgCO2e" for phase in EMISSION_PHASES)
    if include_phase_totals:
        columns.extend(PHASE_TOTAL_COLUMNS)
    return list(dict.fromkeys(columns))


def envelope_intensities_per_m2(env_lookup: Any, *, code: str) -> tuple[float, float, float]:
    """Return per-area embodied intensities from the envelope lookup.

    This is the shared reader for assembly-style envelope entries used by both
    frontends. It returns:

    - production intensity [kgCO2e/m2]
    - demolition / recycling intensity [kgCO2e/m2]
    - biogenic intensity from the legacy `GHG_biogenic_kgCO2m2` field

    The helper deliberately mirrors the current envelope lookup contract:
    if split production / recycling fields are missing, it falls back to the
    legacy aggregate `GHG_kgCO2m2` field and assumes zero demolition.

    Note:
    `GHG_biogenic_kgCO2m2` is forwarded as stored in the database. The backend
    later applies the sign convention for uptake (`biogenic_*` is written as a
    negative contribution for `add` / `replace` events), but it does not do a
    unit conversion here.
    """
    code_str = str(code)
    try:
        production_any = env_lookup.get_item_value(code_str, "GHG_production_kgCO2m2")
        demolition_any = env_lookup.get_item_value(code_str, "GHG_recycling_kgCO2m2")
    except KeyError:
        production_any = env_lookup.get_item_value(code_str, "GHG_kgCO2m2")
        demolition_any = 0.0
    biogenic_any = env_lookup.get_item_value(code_str, "GHG_biogenic_kgCO2m2")
    if production_any is None or demolition_any is None or biogenic_any is None:
        raise ValueError(
            f"Envelope database returned None for one of the required fields for item {code_str}."
        )
    return float(production_any), float(demolition_any), float(biogenic_any)


def _technical_alias_column(phase: str) -> str:
    return f"{phase}_technical_systems_kgCO2e"


def _technical_quarter_columns(phase: str) -> list[str]:
    return [f"{phase}_{component}_kgCO2e" for component in TECHNICAL_SYSTEM_COMPONENTS]


def finalise_emission_timeline_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for column in emission_timeline_columns(
        include_legacy_aliases=True,
        include_phase_totals=True,
    ):
        if column not in out.columns:
            out[column] = 0.0

    for phase in EMISSION_PHASES:
        alias_col = _technical_alias_column(phase)
        quarter_cols = _technical_quarter_columns(phase)
        out[alias_col] = out[quarter_cols].fillna(0.0).sum(axis=1)

        phase_total_col = f"{phase}_kgCO2e"
        phase_component_cols = [
            col
            for col in out.columns
            if isinstance(col, str)
            and col.startswith(f"{phase}_")
            and col.endswith("_kgCO2e")
            and col not in {alias_col, phase_total_col}
        ]
        out[phase_total_col] = (
            out[phase_component_cols].fillna(0.0).sum(axis=1)
            if phase_component_cols
            else 0.0
        )

    return out


def sum_timeline_dataframes_by_building(
    result_list: Sequence[tuple[str, pd.DataFrame]],
) -> pd.DataFrame:
    if not result_list:
        raise ValueError("result_list must be non-empty")

    sample_df = result_list[0][1].copy()
    sample_df = sample_df.drop(columns=[c for c in ("date", "period", "name") if c in sample_df.columns])
    numeric_cols = list(sample_df.select_dtypes(include="number").columns)
    summed_df = pd.DataFrame(
        data=0.0,
        index=[building for building, _ in result_list],
        columns=numeric_cols,
    )
    summed_df.index.rename("name", inplace=True)

    for building, df in result_list:
        df_copy = df.copy()
        df_copy = df_copy.drop(columns=[c for c in ("date", "period", "name") if c in df_copy.columns])
        df_numeric = df_copy.select_dtypes(include="number")
        summed_df.loc[building] = summed_df.loc[building].add(df_numeric.sum(axis=0), fill_value=0.0)

    return summed_df


def sum_timeline_dataframes_by_index(dfs: Sequence[pd.DataFrame]) -> pd.DataFrame:
    if not dfs:
        raise ValueError("dfs must be non-empty")

    sample_df = dfs[0]
    has_date_column = "date" in sample_df.columns
    has_year_index = isinstance(sample_df.index[0], str) and str(sample_df.index[0]).startswith("Y_")

    dfs_for_sum: list[pd.DataFrame] = []
    date_series: pd.Series | None = None
    for df in dfs:
        df_copy = df.copy()
        if "date" in df_copy.columns:
            if date_series is None:
                date_series = df_copy["date"].copy()
            df_copy = df_copy.drop(columns=["date"])
        if "period" in df_copy.columns:
            df_copy = df_copy.drop(columns=["period"])
        if "name" in df_copy.columns:
            df_copy = df_copy.drop(columns=["name"])
        if "Note" in df_copy.columns:
            df_copy = df_copy.drop(columns=["Note"])
        dfs_for_sum.append(df_copy.select_dtypes(include="number"))

    if has_year_index:
        min_year = min(min(years_from_timeline_index(df.index)) for df in dfs_for_sum)
        max_year = max(max(years_from_timeline_index(df.index)) for df in dfs_for_sum)
        reindex_range: Sequence[str] | pd.RangeIndex = [
            f"Y_{year}" for year in range(int(min_year), int(max_year) + 1)
        ]
    else:
        index_min = min(int(df.index.min()) for df in dfs_for_sum)
        index_max = max(int(df.index.max()) for df in dfs_for_sum)
        reindex_range = pd.RangeIndex(index_min, index_max + 1)

    out = (
        pd.concat(dfs_for_sum)
        .groupby(level=0, sort=True)
        .sum()
        .reindex(reindex_range, fill_value=0.0)
    )

    if has_date_column and date_series is not None:
        out_with_index = out.reset_index(drop=True)
        if len(date_series) >= len(out_with_index):
            out_with_index.insert(0, "date", date_series.iloc[:len(out_with_index)].to_list())
        else:
            full_dates = pd.concat([date_series] * (len(out_with_index) // len(date_series) + 1))
            out_with_index.insert(0, "date", full_dates.iloc[:len(out_with_index)].to_list())
        date_cols = ["date"]
        other_cols = [col for col in out_with_index.columns if col not in date_cols]
        return out_with_index[date_cols + other_cols]

    if has_year_index:
        out_with_index = out.reset_index(drop=True)
        out_with_index.insert(0, "period", [str(label) for label in out.index])
        year_cols = ["period"]
        other_cols = [col for col in out_with_index.columns if col not in year_cols]
        return out_with_index[year_cols + other_cols]

    return out.reset_index(drop=True)


def _discount_over_year_indexed(
    base: pd.Series,
    *,
    ref_year: int,
    target_year: int,
    target_fraction: float,
) -> pd.Series:
    if target_year <= ref_year:
        raise ValueError("Target year must be greater than reference year.")
    if target_fraction < 0:
        raise ValueError("Target fraction must be non-negative.")

    years = years_from_timeline_index(base.index)
    series = base.reindex(base.index).astype(float)

    years_arr = np.array(years, dtype=int)
    factors = np.ones(len(base.index), dtype=float)
    mask_linear = (years_arr >= int(ref_year)) & (years_arr <= int(target_year))
    if mask_linear.any():
        factors[mask_linear] = np.linspace(1.0, float(target_fraction), int(mask_linear.sum()))
    factors[years_arr > int(target_year)] = float(target_fraction)
    return series * factors


def _coerce_float(value: Any) -> float:
    return float(value)


@dataclass(frozen=True)
class ResolvedComponent:
    """One fully resolved embodied component contribution.

    This is the smallest embodied unit the backend knows how to write.
    By the time a frontend creates this object, all source-specific lookup is
    already done:

    - the canonical component name is chosen, e.g. `wall_ag` or `technical_system_hs`
    - the affected quantity is known as `area_m2`
    - the per-area emission intensities are already resolved

    The backend does not know whether this came from a normal building lookup,
    an archetype layer snapshot, or a pathway code change. It only multiplies
    these intensities by `area_m2` according to the lifecycle event mode.
    """

    component: str
    area_m2: float
    production_per_area: float = 0.0
    demolition_per_area: float = 0.0
    biogenic_per_area: float = 0.0


@dataclass(frozen=True)
class LifecycleChange:
    """One lifecycle event applied at one year.

    This is the main handoff object from a frontend to the backend.

    A `LifecycleChange` answers:
    - when the event happens: `year`
    - what generic lifecycle rule to apply: `mode`
    - which resolved embodied components are affected: `components`
    - what note should be attached to that year: `note`
    - which renewal clocks should be reset afterwards: `reset_keys`

    Typical examples:
    - initial construction of a building
    - scheduled replacement of a wall or window
    - pathway-authored code change of one supply subsystem
    - demolition at end of life
    """

    year: int
    mode: Literal["add", "replace", "demolish", "production_only"]
    components: tuple[ResolvedComponent, ...]
    note: str | None = None
    reset_keys: tuple[str, ...] = ()


@dataclass(frozen=True)
class OperationalSegment:
    """One constant operational-emission interval.

    Frontends use this to describe yearly operational emissions without
    writing dataframe rows directly. The backend expands the segment over the
    covered year range and writes the operational columns.

    Examples:
    - normal workflow: one carried-forward segment from construction to `year_end`
    - pathway workflow: one segment per state interval
    """

    start_year: int
    end_year: int
    column_values: Mapping[str, float]


@dataclass
class RenewalScheduler:
    """Small helper that tracks the next due year for recurring components.

    The scheduler itself is deliberately generic. It does not know what a wall,
    window, or technical system means. It only tracks keyed due years and helps
    the backend choose the next year that needs attention.
    """

    due_years: dict[str, int]
    end_year: int

    def next_year(self, authored_years: Sequence[int]) -> int | None:
        candidates = [
            int(year)
            for year in authored_years
            if int(year) <= int(self.end_year)
        ]
        candidates.extend(
            int(year)
            for year in self.due_years.values()
            if int(year) <= int(self.end_year)
        )
        return min(candidates) if candidates else None

    def due_keys(self, year: int) -> list[str]:
        return sorted(
            key for key, due_year in self.due_years.items() if int(due_year) == int(year)
        )

    def set_due(self, key: str, year: int | None) -> None:
        if year is None:
            self.due_years.pop(str(key), None)
            return
        self.due_years[str(key)] = int(year)


class TimelineFrontendBase(ABC):
    """Contract for workflow-specific timeline frontends.

    A frontend is responsible for resolving domain data into canonical inputs
    for the backend. In other words, the frontend decides *what exists* and
    *what changes*, while the backend decides *how those changes are written*.

    Current implementations:
    - `NormalTimelineFrontend`: building-resolved normal workflow
    - `PathwayTimelineFrontend`: archetype/log-driven pathway workflow
    """

    @abstractmethod
    def building_name(self) -> str:
        """Return the building identifier for this frontend instance."""
        raise NotImplementedError

    @abstractmethod
    def start_year(self) -> int:
        """Return the first year included in the timeline."""
        raise NotImplementedError

    @abstractmethod
    def end_year(self) -> int:
        """Return the last year included in the timeline."""
        raise NotImplementedError

    @abstractmethod
    def initial_change(self) -> LifecycleChange | None:
        """Return the baseline embodied event that makes the building exist.

        Conceptually this is not a "modification". It is the starting embodied
        state of the building entering the timeline:

        - normal frontend: the building-resolved envelope, window, technical
          system, and PV baseline at the construction year
        - pathway frontend: the archetype/log-resolved snapshot that exists at
          the building construction year

        After this initial event, the scheduler only needs renewals and other
        authored changes.
        """
        raise NotImplementedError

    @abstractmethod
    def initial_due_years(self) -> dict[str, int]:
        """Return the first scheduled renewal year for each recurring key."""
        raise NotImplementedError

    @abstractmethod
    def authored_changes_by_year(self) -> dict[int, tuple[LifecycleChange, ...]]:
        """Return non-recurring authored changes grouped by year."""
        raise NotImplementedError

    @abstractmethod
    def replacement_change(self, key: str, year: int) -> LifecycleChange | None:
        """Return the lifecycle event caused by a scheduled renewal key at `year`."""
        raise NotImplementedError

    @abstractmethod
    def next_due_year(self, key: str, year: int) -> int | None:
        """Return the next due year after a renewal key fires at `year`.

        This stays frontend-owned because the service-life rule may depend on
        source-specific state at that year, e.g. a pathway window code change
        can change the next window lifetime without changing backend logic.
        """
        raise NotImplementedError

    @abstractmethod
    def operational_segments(self) -> tuple[OperationalSegment, ...]:
        """Return operational intervals to be expanded over the timeline."""
        raise NotImplementedError

    @abstractmethod
    def demolition_change(self) -> LifecycleChange | None:
        """Return the demolition event, if the building is demolished in-horizon."""
        raise NotImplementedError


class EmissionTimelineBackend:
    """Assemble one canonical yearly timeline from a frontend description.

    The backend owns the generic mechanics:
    - create the dataframe and canonical columns
    - apply operational segments
    - apply initial/authored/scheduled/demolition lifecycle events
    - rebuild totals and legacy aliases

    It does not load scenario files or interpret archetypes/building metadata
    itself. That remains the frontend's responsibility.
    """

    def __init__(
        self,
        *,
        feedstock_policies: Mapping[str, tuple[int, int, float]] | None = None,
    ) -> None:
        self.feedstock_policies = dict(feedstock_policies or {})

    def build(self, frontend: TimelineFrontendBase) -> pd.DataFrame:
        start_year = int(frontend.start_year())
        end_year = int(frontend.end_year())
        if start_year > end_year:
            raise ValueError("Timeline start year must not be greater than end year.")

        timeline = self._empty_timeline(start_year=start_year, end_year=end_year)
        demolition_change = frontend.demolition_change()
        active_end_year = (
            min(end_year, int(demolition_change.year) - 1)
            if demolition_change is not None
            else end_year
        )

        self._apply_operational_segments(
            timeline=timeline,
            segments=frontend.operational_segments(),
        )

        initial_change = frontend.initial_change()
        self._apply_lifecycle_change(timeline, initial_change)

        authored_changes = dict(frontend.authored_changes_by_year())
        scheduler = RenewalScheduler(
            due_years=dict(frontend.initial_due_years()),
            end_year=active_end_year,
        )

        while True:
            year = scheduler.next_year(sorted(authored_changes.keys()))
            if year is None:
                break

            if year in authored_changes:
                for change in authored_changes.pop(year):
                    self._apply_lifecycle_change(timeline, change)
                    for key in change.reset_keys:
                        scheduler.set_due(key, frontend.next_due_year(key, year))

            for key in scheduler.due_keys(year):
                change = frontend.replacement_change(key, year)
                self._apply_lifecycle_change(timeline, change)
                scheduler.set_due(key, frontend.next_due_year(key, year))

        self._apply_lifecycle_change(timeline, demolition_change)
        return self.finalise(timeline)

    @staticmethod
    def aggregate_by_index(dfs: Sequence[pd.DataFrame]) -> pd.DataFrame:
        return sum_timeline_dataframes_by_index(dfs)

    @staticmethod
    def aggregate_by_building(named_dfs: Sequence[tuple[str, pd.DataFrame]]) -> pd.DataFrame:
        return sum_timeline_dataframes_by_building(named_dfs)

    @staticmethod
    def finalise(df: pd.DataFrame) -> pd.DataFrame:
        return finalise_emission_timeline_dataframe(df)

    @staticmethod
    def feedstock_policies_from_config(
        emissions_cfg: Any,
    ) -> dict[str, tuple[int, int, float]] | None:
        return feedstock_policies_from_config(emissions_cfg)

    @staticmethod
    def _empty_timeline(*, start_year: int, end_year: int) -> pd.DataFrame:
        idx = [f"Y_{year}" for year in range(int(start_year), int(end_year) + 1)]
        columns = emission_timeline_columns(
            include_legacy_aliases=True,
            include_phase_totals=True,
        ) + list(OPERATION_TOTAL_COLUMNS)
        df = pd.DataFrame(0.0, index=idx, columns=columns)
        df.index.name = "period"
        df["Note"] = ""
        return df

    def _apply_lifecycle_change(
        self,
        timeline: pd.DataFrame,
        change: LifecycleChange | None,
    ) -> None:
        if change is None:
            return
        label = f"Y_{int(change.year)}"
        if label not in timeline.index:
            return

        if change.note:
            self._append_note(timeline, year=int(change.year), message=str(change.note))

        include_production = change.mode in {"add", "replace", "production_only"}
        include_demolition = change.mode in {"replace", "demolish"}
        include_biogenic = change.mode in {"add", "replace"}

        for component in change.components:
            area = float(component.area_m2)
            if area <= 0.0:
                continue

            production = float(component.production_per_area) * area if include_production else 0.0
            demolition = float(component.demolition_per_area) * area if include_demolition else 0.0
            biogenic = (-float(component.biogenic_per_area) * area) if include_biogenic else 0.0

            if production:
                self._log_value(
                    timeline,
                    year=label,
                    column=f"production_{component.component}_kgCO2e",
                    value=production,
                )
            if demolition:
                self._log_value(
                    timeline,
                    year=label,
                    column=f"demolition_{component.component}_kgCO2e",
                    value=demolition,
                )
            if biogenic:
                self._log_value(
                    timeline,
                    year=label,
                    column=f"biogenic_{component.component}_kgCO2e",
                    value=biogenic,
                )

    @staticmethod
    def _log_value(
        timeline: pd.DataFrame,
        *,
        year: str,
        column: str,
        value: float,
    ) -> None:
        if column not in timeline.columns:
            timeline[column] = 0.0
        current = timeline.at[year, column]
        timeline.at[year, column] = _coerce_float(current) + _coerce_float(value)

    @staticmethod
    def _append_note(
        timeline: pd.DataFrame,
        *,
        year: int,
        message: str,
    ) -> None:
        label = f"Y_{int(year)}"
        if label not in timeline.index:
            return
        current = str(timeline.at[label, "Note"] or "").strip()
        text = str(message).strip()
        if not text:
            return
        if not current:
            timeline.at[label, "Note"] = text
        elif text not in current:
            timeline.at[label, "Note"] = f"{current} | {text}"

    def _apply_operational_segments(
        self,
        *,
        timeline: pd.DataFrame,
        segments: Sequence[OperationalSegment],
    ) -> None:
        if not segments:
            return

        operational = pd.DataFrame(index=timeline.index)
        for segment in segments:
            segment_start = max(int(segment.start_year), years_from_timeline_index(timeline.index)[0])
            segment_end = min(int(segment.end_year), years_from_timeline_index(timeline.index)[-1])
            if segment_start > segment_end:
                continue
            labels = [f"Y_{year}" for year in range(segment_start, segment_end + 1)]
            for column, value in segment.column_values.items():
                if column not in operational.columns:
                    operational[column] = np.nan
                operational.loc[labels, column] = float(value)

        operational = operational.fillna(0.0).astype(float)
        if not operational.empty and self.feedstock_policies:
            self._apply_feedstock_policies_inplace(operational)

        for demand in _tech_name_mapping.keys():
            cols_demand = [
                column
                for column in operational.columns
                if isinstance(column, str)
                and column.startswith(f"{demand}_")
                and column.endswith("_kgCO2e")
            ]
            if cols_demand:
                timeline.loc[:, f"operation_{demand}_kgCO2e"] = operational[cols_demand].sum(axis=1).to_numpy(dtype=float)

        pv_columns = [
            column
            for column in operational.columns
            if isinstance(column, str) and column.startswith("PV_") and column.endswith("_kgCO2e")
        ]
        for column in pv_columns:
            timeline[column] = operational[column].to_numpy(dtype=float)

    def _apply_feedstock_policies_inplace(self, operational: pd.DataFrame) -> None:
        for raw_key, raw_policy in self.feedstock_policies.items():
            ref_year, target_year, target_fraction = raw_policy
            feedstock_upper = str(raw_key).strip().upper()
            for demand in _tech_name_mapping.keys():
                for column in operational.columns:
                    if not isinstance(column, str):
                        continue
                    prefix = f"{demand}_"
                    suffix = "_kgCO2e"
                    if not (column.startswith(prefix) and column.endswith(suffix)):
                        continue
                    feedstock = column[len(prefix):-len(suffix)]
                    if feedstock.strip().upper() != feedstock_upper:
                        continue
                    operational[column] = _discount_over_year_indexed(
                        operational[column],
                        ref_year=int(ref_year),
                        target_year=int(target_year),
                        target_fraction=float(target_fraction),
                    )

            if feedstock_upper == "GRID":
                for column in operational.columns:
                    if not isinstance(column, str):
                        continue
                    if column.startswith("PV_") and column.endswith("_kgCO2e"):
                        operational[column] = _discount_over_year_indexed(
                            operational[column],
                            ref_year=int(ref_year),
                            target_year=int(target_year),
                            target_fraction=float(target_fraction),
                        )
