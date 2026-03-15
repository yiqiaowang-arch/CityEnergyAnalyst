# State Simulation Package

## Main API
- `simulate_all_states(config: Configuration, timeline_name: str) -> None` - Run Step 4 for every baked `state_{year}` folder.
- `main(config: Configuration) -> None` - Script entrypoint for `state-simulations`.
- `prepare_base_workflow_for_state(config: Configuration, timeline_name: str, year: int) -> list[dict[str, Any]]` - Build the pre-network workflow for one state.
- `prepare_post_demand_workflow_for_state(config: Configuration, timeline_name: str, year: int, state_years: list[int]) -> list[dict[str, Any]]` - Build the network tail plus emissions after demand exists.
- `get_state_service_eligibility(state_locator: InputLocator) -> ServiceEligibility` - Read district-service assignment, demand, and eligible-building sets for `DC` / `DH`.
- `copy_previous_network_for_state(main_locator: InputLocator, state_locator: InputLocator, timeline_name: str, year: int, state_years: list[int]) -> str | None` - Copy the latest reusable network baseline into the current state.

## File Roles
- `main.py` - Human-readable Step 4 entrypoint: orchestration, two-pass state execution, timeline-log updates, and district material timeline generation.
- `workflow_assembly.py` - Build pre-network and post-demand workflow steps, including CRAX / PV decisions and emissions-last ordering.
- `service_checks.py` - DH/DC supply assignment checks, demand checks, and per-state requirement reporting.
- `network_handling.py` - Previous-network lookup, reusable layout copying, and stale current-year network cleanup.

## Settings Ownership
- `state-simulations` owns only timeline selection.
- Runtime settings for `radiation`, `demand`, `photovoltaic`, `network-layout`, `thermal-network`, and `emissions` must be read from their own config sections.
- PV execution is controlled only by `config.emissions.include_pv`.
- Existing-network reuse is additive under `network-layout-mode = augment`. If a later state requests fewer connected buildings than the inherited network contains, the reused layout is kept, `network-layout` warns about the extra inherited buildings, and `thermal-network` simulates only the active consumers in that state.
- Step 4 delegates district building selection and DH service granularity to the state's `supply.csv` by running `network-layout` with `overwrite-supply-settings = False`. This allows `network-layout` to write per-building DH service metadata for `thermal-network`.

## Key Patterns
### DO: Keep the two-pass execution order
```python
base_workflow = workflow_assembly.prepare_base_workflow_for_state(config, timeline_name, year)
state.simulate(config, workflow=base_workflow, mark_simulated=False)

post_demand_workflow = workflow_assembly.prepare_post_demand_workflow_for_state(
    config, timeline_name, year, state_years
)
state.simulate(config, workflow=post_demand_workflow, recorded_workflow=base_workflow + post_demand_workflow)
```

### DO: Resolve network steps only after demand exists
```python
service_eligibility = service_checks.get_state_service_eligibility(state_locator)
required_services = service_checks.get_required_services(service_eligibility)
```

### DO: Derive DH service order from state supply settings
```python
dh_network_services = service_checks.get_state_dh_network_services(state_locator)
```

### DO: Normalise pandas-derived mappings before exposing typed dicts
```python
scale_mapping = build_supply_scale_mapping(supply_db)
```

### DO: Reuse previous networks additively
```python
previous_network_name = network_handling.copy_previous_network_for_state(
    main_locator, state_locator, timeline_name, year, state_years
)
```

### DON'T: Reintroduce pending-mode behaviour in Step 4
```python
# Step 4 always reruns all state years because network reuse makes years interdependent.
```
