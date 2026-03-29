"""Shared emission-column selection helpers for lifecycle timeline plots."""

from __future__ import annotations

from typing import Any

SERVICE_TO_TECH: dict[str, str] = {
    "electricity": "E_sys",
    "space_heating": "Qhs_sys",
    "space_cooling": "Qcs_sys",
    "dhw": "Qww_sys",
}


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    if isinstance(value, tuple):
        return [str(v) for v in value]
    if isinstance(value, str):
        return [v.strip() for v in value.split(",") if v.strip()]
    return [str(value)]


def build_requested_emission_base_columns(plot_config: Any) -> list[str]:
    """Return lifecycle timeline base column names without the `_kgCO2e` suffix."""

    categories = _as_list(getattr(plot_config, "y_category_to_plot", []))
    operation_services = _as_list(getattr(plot_config, "operation_services", []))
    envelope_components = _as_list(getattr(plot_config, "envelope_components", []))

    pv_code_raw = getattr(plot_config, "pv_code", None)
    pv_code = str(pv_code_raw).strip() if pv_code_raw is not None else ""

    requested: list[str] = []

    if "operation" in categories:
        for service in operation_services:
            if service in SERVICE_TO_TECH:
                requested.append(f"operation_{SERVICE_TO_TECH[service]}")
            elif service == "pv_electricity_offset" and pv_code:
                requested.append(f"PV_{pv_code}_GRID_offset")
            elif service == "pv_electricity_export" and pv_code:
                requested.append(f"PV_{pv_code}_GRID_export")

    for phase in ("production", "demolition", "biogenic"):
        if phase not in categories:
            continue
        for component in envelope_components:
            if component == "pv" and pv_code:
                requested.append(f"{phase}_PV_{pv_code}")
            else:
                requested.append(f"{phase}_{component}")

    return list(dict.fromkeys(requested))
