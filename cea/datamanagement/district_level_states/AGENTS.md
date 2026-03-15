# District Level States

## Main API
- `DistrictEventTimeline(config: Configuration, timeline_name: str)` - Timeline manager bound to the main scenario.
- `DistrictEventTimeline.bake_states_from_log() -> None` - Rebuild all `state_{year}` folders from cumulative YAML modifications.
- `DistrictEventTimeline.simulate_states(simulation_mode: Literal["pending", "all"], workflow: list[dict[str, Any]], state_workflows: dict[int, list[dict[str, Any]]] | None = None) -> None` - Generic timeline helper for callers that still need explicit mode control; the Step 4 wrapper no longer exposes a config parameter for this.
- `state_simulation.main.simulate_all_states(config: Configuration, timeline_name: str) -> None` - Run Step 4 for every baked `state_{year}` folder.
- `state_simulation.workflow_assembly.prepare_workflow_for_state(config: Configuration, timeline_name: str, year: int, state_years: list[int]) -> list[dict[str, Any]]` - Build the full Step 4 workflow for one state year.
- `state_simulation.service_checks.get_state_service_eligibility(state_locator: InputLocator) -> ServiceEligibility` - Return district-service assignment, demand, and eligible-building sets for `DC` / `DH`.

## Architecture
- Timeline YAML is the source of truth. State folders are materialised projections plus `.district_timeline_signature.json` status metadata.
- State freshness is status-based only. Do not reintroduce hash or recipe-signature freshness checks for Step 4.
- Step 4 builds a per-state workflow in this order:
  `config -> [radiation | radiation-crax] -> occupancy -> demand -> [photovoltaic] -> [network-layout -> thermal-network] -> emissions`
- Runtime sequencing matters: execute the base workflow through `demand` first, then resolve the post-demand tail (`[network-layout -> thermal-network] -> emissions`) from the newly created `Total_demand.csv`.
- Step 4 now always reruns all state years. Network reuse makes the state sequence interdependent, so the dedicated Step 4 wrapper no longer exposes a pending-only mode.
- Step 4 uses `radiation-crax` only when every building in the state has `void_deck == 0`; otherwise it keeps Daysim. Emit a warning when CRAX is selected so the user can see the engine swap in the log.
- Step 4 skips `photovoltaic` when `config.emissions.include_pv` is `False`, because emissions ignore PV contributions in that configuration.
- Step 4 reads saved settings from the underlying script sections. Do not reintroduce Step 4-specific wrapper parameters for `radiation`, `demand`, `photovoltaic`, `network-layout`, `thermal-network`, or `emissions`.
- Step 4 removes a stale current-year `thermal_network_{year}` folder immediately before rerunning `network-layout`. This keeps interrupted or partial reruns from failing on the network-name collision validator while preserving copied previous-year baseline networks.
- Network steps are conditional. If no eligible district services exist for a state, skip both `network-layout` and `thermal-network`.
- Service eligibility is per service and per state:
  `district supply assignment ∩ positive demand`
- Step 4 always runs `thermal-network` in single-phase mode. Do not route Step 4 through thermal-network multi-phase mode.
- Existing-network reuse is additive. Copy the previous state's root layout plus `DC/DH/layout` subfolders locally, then run `network-layout` with `existing-network` plus `network-layout-mode = augment`.
- Step 4 now lets `network-layout` read connected buildings and DH service granularity from the state's `supply.csv` (`overwrite-supply-settings = False`). This preserves per-building DH service metadata for `thermal-network`.
- When a later state requests fewer connected buildings than the inherited network contains, keep the inherited layout. `network-layout` warns about the extra inherited buildings, while `thermal-network` still simulates only the currently active consumers for that year.

## File Map
- `create_new_timeline.py` - Step 0 entrypoint for selecting or creating the timeline container.
- `district_atomic_changes_define.py` - Step 1 entrypoint for saving reusable atomic change templates.
- `atomic_changes.py` - Atomic change persistence and merge helpers.
- `conflict_detector.py` - Detects overlapping archetype/component edits before they are applied.
- `district_events_apply_changes.py` - Step 2 entrypoint that writes selected changes into timeline YAML for a target year.
- `timeline_log.py` - Canonical YAML load/save helpers.
- `timeline_integrity.py` - Integrity gate between YAML, state folders, and cumulative modifications.
- `bake_states_from_log.py` - Step 3 wrapper around `DistrictEventTimeline.bake_states_from_log()`.
- `state_scenario.py` - State creation, bake, deletion, simulation dispatch, and status metadata handling.
- `state_simulation/` - Step 4 package containing the script entrypoint, workflow assembly, service checks, and network carry-forward logic.
- `district_emission_timeline.py` - District embodied + operational timeline aggregation from state outputs and timeline log.
- `timeline_years.py` - Helpers for required years driven by YAML plus building construction years.

## Key Patterns
### DO: Treat cumulative YAML modifications as source of truth
```python
cumulative = timeline.cumulative_by_year()
year_recipe = cumulative[year]
_apply_state_construction_changes(config, timeline_name, year, year_recipe)
```

### DO: Keep Step 4 emissions last
```python
workflow = prepare_base_workflow_for_state(config, timeline_name, year)
if required_services:
    workflow.append(build_network_layout_step(...))
    workflow.append(build_thermal_network_step(...))
workflow.append(EMISSIONS_STEP)
```

### DO: Let state supply settings drive Step 4 network membership and DH service mix
```python
{
    "overwrite-supply-settings": False,
    "itemised-dh-services": dh_network_services,
    "network-layout-mode": "augment",
    "existing-network": previous_network_name,
}
```

### DON'T: Hardcode DH services or rely on connected-building overrides in Step 4
```python
# Derive DH services from the state's supply settings instead.
# Mixed SH-only / SH+DHW networks need per-building metadata for thermal-network.
```

## Related Files
- `../AGENTS.md` - Project-wide coding and config rules.
- `../../technologies/network_layout/main.py` - Existing-network reconciliation (`validate`, `augment`, `filter`).
- `../../technologies/thermal_network/main.py` - Single-phase vs multi-phase entrypoint rules.
