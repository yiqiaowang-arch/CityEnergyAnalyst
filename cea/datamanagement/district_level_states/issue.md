# Current Issue

## What District Level States Is

District level states is the timeline workflow for evolving one scenario across multiple state years.

Instead of simulating only one static scenario, it:
- stores cumulative district changes by year in a timeline log
- bakes a `state_{year}` scenario folder for each milestone year
- simulates each baked state
- aggregates the state results back into district-level timeline outputs

## General Steps

The workflow is organised as:

1. Step 0: Select or create the timeline container
2. Step 1: Define reusable atomic changes
3. Step 2: Apply selected changes to target years in the timeline log
4. Step 3: Bake `state_{year}` scenario folders from the cumulative log
5. Step 4: Simulate all baked state years and build the district material / emission timeline

The current issue is inside Step 4, specifically where district heating network services are assembled and passed into `network-layout` and `thermal-network`.

## Topic

District-state Step 4 needed finer DH service handling.

The original Step 4 network path treated district heating too coarsely:
- it always passed `itemised-dh-services = ['space_heating', 'domestic_hot_water']`
- it used `overwrite-supply-settings = True`
- that bypassed strict supply-based validation
- and it prevented `network-layout` from writing per-building DH service metadata

As a result, mixed cases were wrong:
- if one building used district SH only
- and another used district SH + DHW
- `thermal-network` could still treat both buildings as using both DH services

## What Was Changed

- Step 4 now lets `network-layout` read district building membership and DH service mix from the state's `supply.csv`
- Step 4 derives the ordered DH service union from `supply.csv`
- Step 4 no longer relies on Step 4-local connected-building overrides for `network-layout`

This fixes the mixed SH-only / SH+DHW case in principle because `thermal-network` can now receive per-building DH service metadata.

## What We Are Discussing Now

The remaining open concern is pure DHW-only district years.

Likely limitation:
- `state_simulation/service_checks.py` still treats DH demand in a coarse way
- `network-layout` still filters DH by `QH_sys_MWhyr`
- `thermal-network` still starts DH demand handling from `QH_sys_MWhyr`

So even though DH service granularity is now better, a state where buildings use district DHW without district space heating may still not run correctly.

## File Structure Needed To Understand This Issue

Read the current issue in three layers.

### 1. Step 4 wrapper layer

- `cea/datamanagement/district_level_states/state_simulation/main.py`
  - Step 4 entrypoint
  - runs the two-pass per-state workflow
- `cea/datamanagement/district_level_states/state_simulation/workflow_assembly.py`
  - builds the actual Step 4 workflow
  - decides what `network-layout` and `thermal-network` receive
  - this is where the latest DH-service fix was wired in
- `cea/datamanagement/district_level_states/state_simulation/service_checks.py`
  - reads the state `supply.csv` and `Total_demand.csv`
  - determines DH/DC eligibility and the state's DH service union
- `cea/datamanagement/district_level_states/state_simulation/network_handling.py`
  - handles additive reuse, stale output cleanup, and previous-network lookup
  - important for network carry-forward, but not the main DH-service bug source

### 2. Network-layout layer

- `cea/technologies/network_layout/main.py`
  - decides whether connected buildings come from `supply.csv` or override parameters
  - validates `itemised-dh-services`
  - writes `building_services.json` when using supply-driven mode
  - this file is the key bridge between Step 4 and thermal-network
- `cea/technologies/network_layout/plant_node_operations.py`
  - converts DH service lists into plant types like `PLANT_hs_ww`
  - affects how thermal-network interprets DH service priority

### 3. Thermal-network / substation layer

- `cea/technologies/thermal_network/main.py`
  - loads `building_services.json` and passes it into the model
- `cea/technologies/thermal_network/simplified/model.py`
  - runs the single-phase thermal-network simulation
  - still appears to gate DH demand with `QH_sys_MWhyr`, which is why DHW-only years remain suspicious
- `cea/technologies/substation.py`
  - applies the per-building DH service split
  - this is where mixed cases like SH-only vs SH+DHW are finally resolved per building

## Suggested Reading Order

If someone needs to debug the current issue quickly, open files in this order:

1. `cea/datamanagement/district_level_states/state_simulation/workflow_assembly.py`
2. `cea/datamanagement/district_level_states/state_simulation/service_checks.py`
3. `cea/technologies/network_layout/main.py`
4. `cea/technologies/thermal_network/main.py`
5. `cea/technologies/thermal_network/simplified/model.py`
6. `cea/technologies/substation.py`

## Next Step

Run live edge-case tests, especially:
- mixed SH-only / SH+DHW
- pure DHW-only district year

If DHW-only fails, the next implementation task is to trace DH demand gating in:
- `cea/technologies/network_layout/main.py`
- `cea/technologies/thermal_network/simplified/model.py`

and decide how district DHW-only networks should be represented end-to-end.
