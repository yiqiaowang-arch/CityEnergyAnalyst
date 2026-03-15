# District Level States Plan

## What District Level States Is

District level states is the timeline workflow that turns one main scenario into a sequence of state-in-time scenarios.

It lets the user:
- define cumulative district changes by year
- materialise those changes into baked `state_{year}` folders
- simulate each baked state
- roll the results back up into district timeline outputs

## General Steps

1. Step 0: Select or create the timeline container
2. Step 1: Define reusable atomic changes
3. Step 2: Apply changes to years in the timeline log
4. Step 3: Bake `state_{year}` scenario folders
5. Step 4: Simulate all baked states and generate district material / emission timeline outputs

## Current Topic

Step 4 (`state-simulations`) now supports:
- two-pass per-state execution
- conditional `network-layout -> thermal-network`
- additive reuse of previous networks
- supply-driven DH service granularity for mixed SH-only / SH+DHW cases

The active discussion has moved from workflow ordering to DH service correctness inside reused district heating networks.

## Implemented So Far

### Step 4 architecture

- Step 4 entrypoint now lives in:
  - `cea/datamanagement/district_level_states/state_simulation/main.py`
- Supporting modules:
  - `workflow_assembly.py`
  - `service_checks.py`
  - `network_handling.py`

### Runtime sequence

Per state year:
- `config`
- `radiation` or `radiation-crax`
- `occupancy`
- `demand`
- optional `photovoltaic`
- optional `network-layout`
- optional `thermal-network`
- `emissions`

Emissions remain last.

### State execution model

- Step 4 always reruns all baked state years.
- `pending` is no longer part of the Step 4 wrapper.
- Base simulations run first.
- Network eligibility is resolved only after `demand` creates `Total_demand.csv`.

### Network reuse

- Reuse is still additive.
- Step 4 looks backwards for the latest earlier state year with a reusable network.
- Interrupted reruns pre-clean only the current-year `thermal_network_{year}` folder.

### PV and radiation shortcuts

- `photovoltaic` is skipped when `config.emissions.include_pv` is `False`
- `radiation-crax` is used when all buildings in the state have `void_deck == 0`

### DH service fix

The previous Step 4 path forced:
- `overwrite-supply-settings = True`
- `itemised-dh-services = ['space_heating', 'domestic_hot_water']`

That was too coarse and caused mixed SH-only / SH+DHW buildings to be misrepresented.

Current fix:
- Step 4 now runs `network-layout` with `overwrite-supply-settings = False`
- DH network services are derived from the state's `supply.csv`
- per-building DH service metadata is preserved for `thermal-network`

This should correctly support cases like:
- building A: district SH + district DHW
- building B: district SH only

## Files Touched In The Latest Phase

- `cea/datamanagement/district_level_states/state_simulation/service_checks.py`
- `cea/datamanagement/district_level_states/state_simulation/workflow_assembly.py`
- `cea/tests/test_state_simulations.py`
- `cea/scripts.yml`
- `cea/datamanagement/district_level_states/AGENTS.md`
- `cea/datamanagement/district_level_states/state_simulation/AGENTS.md`
- `cea/AGENTS.md`

## Verification

Ran:

```powershell
pixi run python -m py_compile cea\datamanagement\district_level_states\state_simulation\service_checks.py cea\datamanagement\district_level_states\state_simulation\workflow_assembly.py cea\tests\test_state_simulations.py
pixi run python -m pytest cea\tests\test_script_parameters.py cea\tests\test_state_simulations.py
```

Result:
- `15 passed`

Known non-failing noise:
- pytest cache warning on Windows
- `InputLocator` temp-directory cleanup warnings at process exit

## Current Open Issue

See `issue.md`.

Short version:
- mixed SH-only / SH+DHW is now addressed
- pure DHW-only district years are still likely limited by the underlying network stack

## Fast Reading Order

For the current issue, start with:

1. `state_simulation/workflow_assembly.py`
2. `state_simulation/service_checks.py`
3. `../../technologies/network_layout/main.py`
4. `../../technologies/thermal_network/main.py`
5. `../../technologies/thermal_network/simplified/model.py`
6. `../../technologies/substation.py`

## Recommended Next Checks

1. Run a live Step 4 scenario again after the latest DH-service fix.
2. Test these edge cases explicitly:
   - mixed SH-only / SH+DHW DH buildings
   - DH service disappears for one year, then returns later
   - partial building exit from an active reused DH network
   - pure DHW-only district year
3. If pure DHW-only still fails, inspect:
   - `cea/technologies/network_layout/main.py`
   - `cea/technologies/thermal_network/simplified/model.py`
   because both still use `QH_sys_MWhyr` as the DH demand gate.
