# Tests

## Main API
- `cea/tests/test_state_simulations.py` - Regression coverage for district-level Step 4 workflow assembly and service checks.
- `cea/tests/test_script_parameters.py` - Script-parameter surface checks for `scripts.yml`.
- `cea/tests/test_config.py` - Configuration and parameter-behaviour regression tests.

## Key Patterns
### DO: Keep scenario tests isolated
```python
self.temp_dir = tempfile.mkdtemp(dir=os.getcwd())
self.config.scenario = self.temp_dir
```

### DO: Patch expensive or noisy filesystem cleanup in `InputLocator`
```python
self.inputlocator_cleanup = patch.object(
    InputLocator,
    "_cleanup_temp_directory",
    return_value=None,
)
```

### DO: Mock CSV inputs directly for Step 4 service checks
```python
with patch.object(service_checks.pd, "read_csv", side_effect=fake_read_csv):
    eligibility = service_checks.get_state_service_eligibility(locator)
```

### DO: Use temporary scenario folders or patched locators for network-reuse path tests
```python
state_folder = locator.get_state_in_time_scenario_folder("timeline", 2025)
os.makedirs(state_folder, exist_ok=True)
```

### DO: Keep network-layout output tests in memory on Windows
```python
monkeypatch.setattr(gpd.GeoDataFrame, "to_file", fake_to_file, raising=False)
monkeypatch.setattr("builtins.open", fake_open)
```

### DON'T: Put tests outside `cea/tests/`
```python
# New regression tests belong here, even for deeply nested modules.
```

## Related Files
- `../datamanagement/district_level_states/state_simulation/service_checks.py` - Step 4 district-service eligibility logic.
- `../datamanagement/district_level_states/state_simulation/workflow_assembly.py` - Step 4 workflow construction.
- `../technologies/network_layout/main.py` - Service-specific layout pruning and metadata refresh for reused networks.
