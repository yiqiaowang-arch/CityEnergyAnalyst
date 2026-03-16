from __future__ import annotations

from typing import TypedDict

import pandas as pd

from cea.inputlocator import InputLocator

DISTRICT_SERVICES = ("DC", "DH")
SERVICE_LABELS = {"DC": "cooling", "DH": "heating"}
DH_SERVICE_ORDER = ("space_heating", "domestic_hot_water")
DH_SERVICE_DEMAND_FIELDS = {
    "space_heating": "Qhs_sys_MWhyr",
    "domestic_hot_water": "Qww_sys_MWhyr",
}
LEGACY_DH_DEMAND_FIELD = "QH_sys_MWhyr"


class ServiceBuildingSets(TypedDict):
    district_buildings: list[str]
    demand_buildings: list[str]
    eligible_buildings: list[str]


ServiceEligibility = dict[str, ServiceBuildingSets]


def build_supply_scale_mapping(supply_db: pd.DataFrame) -> dict[str, str]:
    """Convert a supply database table into an explicit `code -> scale` mapping."""
    return {
        str(code): str(scale)
        for code, scale in zip(supply_db["code"], supply_db["scale"])
    }


def get_state_service_eligibility(state_locator: InputLocator) -> ServiceEligibility:
    """Return district-service assignment, demand, and eligible-building sets for one state."""
    service_eligibility: ServiceEligibility = {}

    for network_type in DISTRICT_SERVICES:
        eligibility: ServiceBuildingSets = {
            "district_buildings": [],
            "demand_buildings": [],
            "eligible_buildings": [],
        }
        try:
            district_buildings = sorted(
                set(get_buildings_with_district_service(state_locator, network_type))
            )
            demand_buildings = sorted(
                set(get_buildings_with_demand_for_service(state_locator, network_type))
            )
            eligible_buildings = sorted(set(district_buildings) & set(demand_buildings))

            eligibility["district_buildings"] = district_buildings
            eligibility["demand_buildings"] = demand_buildings
            eligibility["eligible_buildings"] = eligible_buildings
        except FileNotFoundError as exc:
            print(f"Warning: Could not check {network_type} requirements - {exc}")
        except Exception as exc:
            print(f"Warning: Error checking {network_type} requirements - {exc}")

        service_eligibility[network_type] = eligibility

    return service_eligibility


def get_required_services(service_eligibility: ServiceEligibility) -> list[str]:
    """Return the district services whose eligible-building sets are non-empty."""
    return [
        network_type
        for network_type in DISTRICT_SERVICES
        if service_eligibility[network_type]["eligible_buildings"]
    ]


def report_service_requirements(year: int, service_eligibility: ServiceEligibility) -> None:
    """Print the per-service requirement summary for one state year."""
    for network_type in DISTRICT_SERVICES:
        district_buildings = service_eligibility[network_type]["district_buildings"]
        eligible_buildings = service_eligibility[network_type]["eligible_buildings"]
        service_name = SERVICE_LABELS[network_type]

        if eligible_buildings:
            print(
                f"State {year}: {network_type} network required for "
                f"{len(eligible_buildings)} building(s) with district {service_name} and positive demand."
            )
            continue

        if district_buildings:
            print(
                f"State {year}: Skipping {network_type} network because no buildings have both "
                f"district {service_name} supply and positive demand."
            )


def get_service_scale_mapping(locator: InputLocator, network_type: str) -> dict[str, str]:
    """Load the supply database scale mapping for the requested district service."""
    if network_type == "DH":
        heating_supply_db = pd.read_csv(locator.get_database_assemblies_supply_heating())
        hot_water_supply_db = pd.read_csv(locator.get_database_assemblies_supply_hot_water())
        scale_mapping = build_supply_scale_mapping(heating_supply_db)
        scale_mapping.update(build_supply_scale_mapping(hot_water_supply_db))
        return scale_mapping

    cooling_supply_db = pd.read_csv(locator.get_database_assemblies_supply_cooling())
    return build_supply_scale_mapping(cooling_supply_db)


def order_dh_services(services: set[str]) -> list[str]:
    """Return DH services in CEA's canonical priority order."""
    ordered_services = [service for service in DH_SERVICE_ORDER if service in services]
    ordered_services.extend(sorted(service for service in services if service not in DH_SERVICE_ORDER))
    return ordered_services


def get_buildings_with_district_dh_services(locator: InputLocator) -> dict[str, set[str]]:
    """Return the district DH services configured for each building in `supply.csv`."""
    supply_properties_df = pd.read_csv(locator.get_building_supply())
    scale_mapping = get_service_scale_mapping(locator, "DH")

    district_dh_services_by_building: dict[str, set[str]] = {}
    for _, row in supply_properties_df.iterrows():
        building_name = row["name"]
        building_services: set[str] = set()

        hs_code = row.get("supply_type_hs")
        if hs_code and not pd.isna(hs_code) and scale_mapping.get(hs_code) == "DISTRICT":
            building_services.add("space_heating")

        dhw_code = row.get("supply_type_dhw")
        if dhw_code and not pd.isna(dhw_code) and scale_mapping.get(dhw_code) == "DISTRICT":
            building_services.add("domestic_hot_water")

        if building_services:
            district_dh_services_by_building[building_name] = building_services

    return district_dh_services_by_building


def get_state_dh_network_services(locator: InputLocator) -> list[str]:
    """Return the ordered union of district DH services present in the state's `supply.csv`."""
    network_services: set[str] = set()
    for building_services in get_buildings_with_district_dh_services(locator).values():
        network_services.update(building_services)

    return order_dh_services(network_services)


def get_buildings_with_district_service(locator: InputLocator, network_type: str) -> list[str]:
    """Return buildings configured for district heating or district cooling in `supply.csv`."""
    if network_type == "DH":
        return list(get_buildings_with_district_dh_services(locator))

    supply_properties_df = pd.read_csv(locator.get_building_supply())
    scale_mapping = get_service_scale_mapping(locator, network_type)
    buildings_with_district = []
    for _, row in supply_properties_df.iterrows():
        building_name = row["name"]
        has_district_service = False

        if network_type == "DC" and "supply_type_cs" in supply_properties_df.columns:
            cs_code = row.get("supply_type_cs")
            if cs_code and not pd.isna(cs_code) and scale_mapping.get(cs_code) == "DISTRICT":
                has_district_service = True

        if has_district_service:
            buildings_with_district.append(building_name)

    return buildings_with_district


def get_buildings_with_positive_demand(
    total_demand: pd.DataFrame,
    demand_field: str,
) -> set[str]:
    """Return buildings with positive annual demand in one total-demand column."""
    if demand_field not in total_demand.columns:
        return set()

    demand_values = pd.to_numeric(total_demand[demand_field], errors="coerce").fillna(0.0)
    return set(total_demand.loc[demand_values > 0.0, "name"].tolist())


def get_buildings_with_demand_for_dh_services(locator: InputLocator, total_demand: pd.DataFrame) -> list[str]:
    """Return DH buildings whose assigned district services have positive annual demand."""
    district_dh_services_by_building = get_buildings_with_district_dh_services(locator)

    demand_buildings_by_service = {
        service: get_buildings_with_positive_demand(total_demand, demand_field)
        for service, demand_field in DH_SERVICE_DEMAND_FIELDS.items()
    }

    if LEGACY_DH_DEMAND_FIELD in total_demand.columns:
        legacy_demand_buildings = get_buildings_with_positive_demand(total_demand, LEGACY_DH_DEMAND_FIELD)
        for service in DH_SERVICE_DEMAND_FIELDS:
            if not demand_buildings_by_service[service]:
                demand_buildings_by_service[service] = legacy_demand_buildings

    demand_buildings = [
        building_name
        for building_name, building_services in district_dh_services_by_building.items()
        if any(
            building_name in demand_buildings_by_service[service]
            for service in building_services
        )
    ]
    return sorted(demand_buildings)


def get_buildings_with_demand_for_service(locator: InputLocator, network_type: str) -> list[str]:
    """Return buildings with positive annual demand for the requested district service."""
    total_demand = pd.read_csv(locator.get_total_demand())

    if network_type == "DH":
        return get_buildings_with_demand_for_dh_services(locator, total_demand)
    else:
        demand_field = "QC_sys_MWhyr"

    return sorted(get_buildings_with_positive_demand(total_demand, demand_field))
