import os
import shutil
import tempfile
import unittest
from unittest.mock import Mock, patch

import pandas as pd

from cea.config import Configuration
from cea.datamanagement.district_level_states.state_simulation import main as state_simulation_main
from cea.datamanagement.district_level_states.state_simulation import network_handling
from cea.datamanagement.district_level_states.state_simulation import workflow_assembly
from cea.inputlocator import InputLocator


def make_service_eligibility(
    *,
    dc_eligible: list[str] | None = None,
    dh_eligible: list[str] | None = None,
    dc_district: list[str] | None = None,
    dh_district: list[str] | None = None,
) -> dict[str, dict[str, list[str]]]:
    dc_eligible = sorted(dc_eligible or [])
    dh_eligible = sorted(dh_eligible or [])
    return {
        "DC": {
            "district_buildings": sorted(dc_district if dc_district is not None else dc_eligible),
            "demand_buildings": list(dc_eligible),
            "eligible_buildings": list(dc_eligible),
        },
        "DH": {
            "district_buildings": sorted(dh_district if dh_district is not None else dh_eligible),
            "demand_buildings": list(dh_eligible),
            "eligible_buildings": list(dh_eligible),
        },
    }


def get_step_names(workflow: list[dict]) -> list[str]:
    return [step.get("script", "config") for step in workflow]


class TestStateSimulationsWorkflow(unittest.TestCase):
    def setUp(self):
        self.inputlocator_cleanup = patch.object(
            InputLocator,
            "_cleanup_temp_directory",
            return_value=None,
        )
        self.inputlocator_cleanup.start()
        self.temp_dir = tempfile.mkdtemp(dir=os.getcwd())
        self.config = Configuration()
        self.config.scenario = self.temp_dir
        self.config.emissions.include_pv = False

    def tearDown(self):
        self.inputlocator_cleanup.stop()
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_prepare_workflow_for_state_without_services_skips_network_steps(self):
        with patch.object(
            workflow_assembly.service_checks,
            "get_state_service_eligibility",
            return_value=make_service_eligibility(),
        ), patch.object(
            workflow_assembly,
            "should_use_crax_radiation",
            return_value=False,
        ):
            workflow = workflow_assembly.prepare_workflow_for_state(
                self.config,
                timeline_name="timeline",
                year=2030,
                state_years=[2030],
            )

        self.assertEqual(
            get_step_names(workflow),
            ["config", "radiation", "occupancy", "demand", "emissions"],
        )
        self.assertEqual(workflow[-1]["script"], "emissions")

    def test_prepare_workflow_for_state_uses_crax_and_warns_when_no_void_decks_exist(self):
        with patch.object(
            workflow_assembly.service_checks,
            "get_state_service_eligibility",
            return_value=make_service_eligibility(),
        ), patch.object(
            workflow_assembly,
            "should_use_crax_radiation",
            return_value=True,
        ), patch("builtins.print") as print_mock:
            workflow = workflow_assembly.prepare_workflow_for_state(
                self.config,
                timeline_name="timeline",
                year=2030,
                state_years=[2030],
            )

        self.assertEqual(
            get_step_names(workflow),
            ["config", "radiation-crax", "occupancy", "demand", "emissions"],
        )
        printed_messages = "\n".join(
            " ".join(str(arg) for arg in call.args) for call in print_mock.call_args_list
        )
        self.assertIn("Using radiation-crax instead of radiation (DAYSIM)", printed_messages)

    def test_prepare_workflow_for_state_keeps_photovoltaic_when_emissions_include_pv(self):
        self.config.emissions.include_pv = True

        with patch.object(
            workflow_assembly.service_checks,
            "get_state_service_eligibility",
            return_value=make_service_eligibility(),
        ), patch.object(
            workflow_assembly,
            "should_use_crax_radiation",
            return_value=False,
        ):
            workflow = workflow_assembly.prepare_workflow_for_state(
                self.config,
                timeline_name="timeline",
                year=2030,
                state_years=[2030],
            )

        self.assertEqual(
            get_step_names(workflow),
            ["config", "radiation", "occupancy", "demand", "photovoltaic", "emissions"],
        )

    def test_prepare_workflow_for_state_service_cases(self):
        cases = [
            (
                "dh_only",
                make_service_eligibility(dh_eligible=["B2"]),
                ["DH"],
                [],
                ["B2"],
            ),
            (
                "dc_only",
                make_service_eligibility(dc_eligible=["B1"]),
                ["DC"],
                ["B1"],
                [],
            ),
            (
                "both",
                make_service_eligibility(dc_eligible=["B1"], dh_eligible=["B2"]),
                ["DC", "DH"],
                ["B1"],
                ["B2"],
            ),
        ]

        for case_name, service_eligibility, expected_services, expected_dc, expected_dh in cases:
            with self.subTest(case_name=case_name), patch.object(
                workflow_assembly.service_checks,
                "get_state_service_eligibility",
                return_value=service_eligibility,
            ), patch.object(
                workflow_assembly.service_checks,
                "get_state_dh_network_services",
                return_value=["space_heating", "domestic_hot_water"],
            ), patch.object(
                workflow_assembly,
                "should_use_crax_radiation",
                return_value=False,
            ), patch.object(
                workflow_assembly.network_handling,
                "copy_previous_network_for_state",
                return_value=None,
            ):
                workflow = workflow_assembly.prepare_workflow_for_state(
                    self.config,
                    timeline_name="timeline",
                    year=2035,
                    state_years=[2035],
                )

            self.assertEqual(
                get_step_names(workflow),
                [
                    "config",
                    "radiation",
                    "occupancy",
                    "demand",
                    "network-layout",
                    "thermal-network",
                    "emissions",
                ],
            )

            network_layout_step = workflow[-3]
            thermal_network_step = workflow[-2]

            self.assertEqual(network_layout_step["script"], "network-layout")
            self.assertEqual(network_layout_step["parameters"]["include-services"], expected_services)
            self.assertFalse(network_layout_step["parameters"]["overwrite-supply-settings"])
            self.assertEqual(network_layout_step["parameters"]["network-layout-mode"], "augment")
            self.assertTrue(network_layout_step["parameters"]["auto-modify-network"])
            self.assertIsNone(network_layout_step["parameters"]["existing-network"])
            self.assertNotIn("cooling-connected-buildings", network_layout_step["parameters"])
            self.assertNotIn("heating-connected-buildings", network_layout_step["parameters"])
            if "DH" in expected_services:
                self.assertEqual(
                    network_layout_step["parameters"]["itemised-dh-services"],
                    ["space_heating", "domestic_hot_water"],
                )
            else:
                self.assertNotIn("itemised-dh-services", network_layout_step["parameters"])

            self.assertEqual(thermal_network_step["script"], "thermal-network")
            self.assertEqual(
                thermal_network_step["parameters"]["network-name"],
                ["thermal_network_2035"],
            )
            self.assertEqual(
                thermal_network_step["parameters"]["network-type"],
                expected_services,
            )
            self.assertFalse(thermal_network_step["parameters"]["multi-phase-mode"])
            self.assertEqual(workflow[-1]["script"], "emissions")

    def test_prepare_workflow_for_state_reuses_previous_network_with_augment_parameters(self):
        timeline_name = "timeline"
        state_years = [2025, 2040]
        expected_previous_network = "thermal_network_2025"

        with patch.object(
            workflow_assembly.service_checks,
            "get_state_service_eligibility",
            return_value=make_service_eligibility(dc_eligible=["B1"], dh_eligible=["B2"]),
        ), patch.object(
            workflow_assembly.service_checks,
            "get_state_dh_network_services",
            return_value=["space_heating"],
        ), patch.object(
            workflow_assembly,
            "should_use_crax_radiation",
            return_value=False,
        ), patch.object(
            workflow_assembly.network_handling,
            "copy_previous_network_for_state",
            return_value=expected_previous_network,
        ):
            workflow = workflow_assembly.prepare_workflow_for_state(
                self.config,
                timeline_name=timeline_name,
                year=2040,
                state_years=state_years,
            )

        network_layout_step = workflow[-3]
        thermal_network_step = workflow[-2]

        self.assertEqual(
            network_layout_step["parameters"]["existing-network"],
            expected_previous_network,
        )
        self.assertEqual(network_layout_step["parameters"]["network-layout-mode"], "augment")
        self.assertFalse(network_layout_step["parameters"]["overwrite-supply-settings"])
        self.assertEqual(
            network_layout_step["parameters"]["itemised-dh-services"],
            ["space_heating"],
        )
        self.assertNotIn("heating-connected-buildings", network_layout_step["parameters"])
        self.assertEqual(
            network_layout_step["parameters"]["include-services"],
            ["DC", "DH"],
        )
        self.assertEqual(
            thermal_network_step["parameters"]["network-name"],
            ["thermal_network_2040"],
        )
        self.assertEqual(
            thermal_network_step["parameters"]["network-type"],
            ["DC", "DH"],
        )

    def test_get_state_dh_network_services_orders_union_from_supply(self):
        locator = Mock()
        locator.get_building_supply.return_value = "supply.csv"

        with patch.object(
            workflow_assembly.service_checks.pd,
            "read_csv",
            return_value=pd.DataFrame(
                [
                    {"name": "B1", "supply_type_hs": "HS_DISTRICT", "supply_type_dhw": "DHW_BUILDING"},
                    {"name": "B2", "supply_type_hs": "HS_DISTRICT", "supply_type_dhw": "DHW_DISTRICT"},
                    {"name": "B3", "supply_type_hs": "HS_BUILDING", "supply_type_dhw": "DHW_BUILDING"},
                ]
            ),
        ), patch.object(
            workflow_assembly.service_checks,
            "get_service_scale_mapping",
            return_value={
                "HS_DISTRICT": "DISTRICT",
                "HS_BUILDING": "BUILDING",
                "DHW_DISTRICT": "DISTRICT",
                "DHW_BUILDING": "BUILDING",
            },
        ):
            dh_services = workflow_assembly.service_checks.get_state_dh_network_services(locator)

        self.assertEqual(dh_services, ["space_heating", "domestic_hot_water"])

    def test_get_state_dh_network_services_returns_space_heating_only_when_no_district_dhw_exists(self):
        locator = Mock()
        locator.get_building_supply.return_value = "supply.csv"

        with patch.object(
            workflow_assembly.service_checks.pd,
            "read_csv",
            return_value=pd.DataFrame(
                [
                    {"name": "B1", "supply_type_hs": "HS_DISTRICT", "supply_type_dhw": "DHW_BUILDING"},
                    {"name": "B2", "supply_type_hs": "HS_DISTRICT", "supply_type_dhw": "DHW_BUILDING"},
                ]
            ),
        ), patch.object(
            workflow_assembly.service_checks,
            "get_service_scale_mapping",
            return_value={
                "HS_DISTRICT": "DISTRICT",
                "DHW_BUILDING": "BUILDING",
            },
        ):
            dh_services = workflow_assembly.service_checks.get_state_dh_network_services(locator)

        self.assertEqual(dh_services, ["space_heating"])

    def test_copy_previous_network_for_state_uses_expected_source_and_target_paths(self):
        timeline_name = "timeline"
        state_years = [2025, 2040]
        expected_previous_network = "thermal_network_2025"
        locator = InputLocator(self.config.scenario)
        state_locator = InputLocator(
            locator.get_state_in_time_scenario_folder(timeline_name, 2040)
        )
        expected_source = InputLocator(
            locator.get_state_in_time_scenario_folder(timeline_name, 2025)
        ).get_thermal_network_folder_network_name_folder(expected_previous_network)
        expected_target = state_locator.get_thermal_network_folder_network_name_folder(
            expected_previous_network
        )

        with patch.object(
            network_handling,
            "find_previous_network",
            return_value=(expected_previous_network, 2025),
        ), patch.object(
            network_handling,
            "copy_previous_network_layout",
        ) as copy_mock:
            previous_network_name = network_handling.copy_previous_network_for_state(
                main_locator=locator,
                state_locator=state_locator,
                timeline_name=timeline_name,
                year=2040,
                state_years=state_years,
            )

        self.assertEqual(previous_network_name, expected_previous_network)
        copy_mock.assert_called_once_with(expected_source, expected_target)

    def test_copy_previous_network_layout_keeps_service_layouts_but_skips_other_service_outputs(self):
        source_network_folder = "C:/scenario/source_network"
        target_network_folder = "C:/scenario/target_network"

        def fake_isdir(path):
            return path in {
                os.path.join(source_network_folder, "DH"),
                os.path.join(source_network_folder, "DC"),
                os.path.join(source_network_folder, "DH", "layout"),
                os.path.join(source_network_folder, "DC", "layout"),
            }

        def fake_isfile(path):
            return path == os.path.join(source_network_folder, "layout.shp")

        with patch.object(
            network_handling.os,
            "listdir",
            return_value=["layout.shp", "DH", "DC"],
        ), patch.object(
            network_handling.os.path,
            "isdir",
            side_effect=fake_isdir,
        ), patch.object(
            network_handling.os.path,
            "isfile",
            side_effect=fake_isfile,
        ), patch.object(
            network_handling.os,
            "makedirs",
        ), patch.object(
            network_handling.shutil,
            "copy2",
        ) as copy2_mock, patch.object(
            network_handling.shutil,
            "copytree",
        ) as copytree_mock:
            network_handling.copy_previous_network_layout(
                source_network_folder=source_network_folder,
                target_network_folder=target_network_folder,
            )

        copy2_mock.assert_called_once_with(
            os.path.join(source_network_folder, "layout.shp"),
            os.path.join(target_network_folder, "layout.shp"),
        )
        self.assertEqual(
            copytree_mock.call_args_list,
            [
                unittest.mock.call(
                    os.path.join(source_network_folder, "DH", "layout"),
                    os.path.join(target_network_folder, "DH", "layout"),
                    dirs_exist_ok=True,
                ),
                unittest.mock.call(
                    os.path.join(source_network_folder, "DC", "layout"),
                    os.path.join(target_network_folder, "DC", "layout"),
                    dirs_exist_ok=True,
                ),
            ],
        )

    def test_cleanup_current_network_outputs_removes_stale_folder(self):
        locator = Mock()
        stale_network_folder = "C:/scenario/state_2030/outputs/data/thermal-network/thermal_network_2030"
        locator.get_thermal_network_folder_network_name_folder.return_value = stale_network_folder

        with patch.object(
            network_handling.os.path,
            "exists",
            return_value=True,
        ), patch.object(
            network_handling.shutil,
            "rmtree",
        ) as rmtree_mock, patch("builtins.print") as print_mock:
            network_handling.cleanup_current_network_outputs(locator, 2030)

        locator.get_thermal_network_folder_network_name_folder.assert_called_once_with(
            "thermal_network_2030"
        )
        rmtree_mock.assert_called_once_with(stale_network_folder)
        printed_messages = "\n".join(
            " ".join(str(arg) for arg in call.args) for call in print_mock.call_args_list
        )
        self.assertIn("Removing the stale current-year folder", printed_messages)

    def test_prepare_workflow_for_state_cleans_stale_outputs_before_network_steps(self):
        with patch.object(
            workflow_assembly.service_checks,
            "get_state_service_eligibility",
            return_value=make_service_eligibility(dh_eligible=["B2"]),
        ), patch.object(
            workflow_assembly.service_checks,
            "get_state_dh_network_services",
            return_value=["space_heating"],
        ), patch.object(
            workflow_assembly,
            "should_use_crax_radiation",
            return_value=False,
        ), patch.object(
            workflow_assembly.network_handling,
            "copy_previous_network_for_state",
            return_value=None,
        ), patch.object(
            workflow_assembly.network_handling,
            "cleanup_current_network_outputs",
        ) as cleanup_mock:
            workflow = workflow_assembly.prepare_workflow_for_state(
                self.config,
                timeline_name="timeline",
                year=2030,
                state_years=[2030],
            )

        cleanup_mock.assert_called_once()
        self.assertIn("network-layout", get_step_names(workflow))

    def test_simulate_all_states_builds_post_demand_workflow_after_base_workflow(self):
        simulate_calls = []
        base_workflow = [
            {"config": "."},
            {"script": "radiation-crax"},
            {"script": "occupancy"},
            {"script": "demand"},
        ]
        post_demand_workflow = [{"script": "emissions", "parameters": {}}]

        class FakeTimeline:
            def __init__(self):
                self.timeline_name = "timeline"
                self.main_locator = object()
                self.log_data = {2030: {}}

            def list_state_years_on_disk(self):
                return [2030]

            def save(self):
                return None

        class FakeState:
            def __init__(self, *args, **kwargs):
                pass

            def simulate(self, config, *, workflow, mark_simulated=True, recorded_workflow=None):
                simulate_calls.append(
                    {
                        "workflow": workflow,
                        "mark_simulated": mark_simulated,
                        "recorded_workflow": recorded_workflow,
                    }
                )

        def fake_prepare_post_demand_workflow_for_state(*args, **kwargs):
            self.assertEqual(len(simulate_calls), 1)
            return post_demand_workflow

        fake_timeline = FakeTimeline()

        with patch.object(
            state_simulation_main,
            "DistrictEventTimeline",
            return_value=fake_timeline,
        ), patch.object(
            state_simulation_main,
            "DistrictStateYear",
            FakeState,
        ), patch.object(
            state_simulation_main,
            "check_district_timeline_log_yaml_integrity",
        ), patch.object(
            state_simulation_main.workflow_assembly,
            "prepare_base_workflow_for_state",
            return_value=base_workflow,
        ), patch.object(
            state_simulation_main.workflow_assembly,
            "prepare_post_demand_workflow_for_state",
            side_effect=fake_prepare_post_demand_workflow_for_state,
        ):
            state_simulation_main.simulate_all_states(self.config, timeline_name="timeline")

        self.assertEqual(len(simulate_calls), 2)
        self.assertEqual(simulate_calls[0]["workflow"], base_workflow)
        self.assertFalse(simulate_calls[0]["mark_simulated"])
        self.assertIsNone(simulate_calls[0]["recorded_workflow"])
        self.assertEqual(simulate_calls[1]["workflow"], post_demand_workflow)
        self.assertTrue(simulate_calls[1]["mark_simulated"])
        self.assertEqual(
            simulate_calls[1]["recorded_workflow"],
            base_workflow + post_demand_workflow,
        )
        self.assertEqual(
            fake_timeline.log_data[2030]["simulation_workflow"],
            base_workflow + post_demand_workflow,
        )

    def test_simulate_all_states_always_runs_all_years(self):
        simulate_calls = []

        class FakeTimeline:
            def __init__(self):
                self.timeline_name = "timeline"
                self.main_locator = object()
                self.log_data = {2030: {}}

            def list_state_years_on_disk(self):
                return [2030]

            def save(self):
                return None

        class FakeState:
            def __init__(self, *args, **kwargs):
                pass

            def simulate(self, config, *, workflow, mark_simulated=True, recorded_workflow=None):
                simulate_calls.append(
                    {
                        "workflow": workflow,
                        "mark_simulated": mark_simulated,
                        "recorded_workflow": recorded_workflow,
                    }
                )

        with patch.object(
            state_simulation_main,
            "DistrictEventTimeline",
            return_value=FakeTimeline(),
        ), patch.object(
            state_simulation_main,
            "DistrictStateYear",
            FakeState,
        ), patch.object(
            state_simulation_main,
            "check_district_timeline_log_yaml_integrity",
            return_value=None,
        ), patch.object(
            state_simulation_main.workflow_assembly,
            "prepare_base_workflow_for_state",
            return_value=[{"config": "."}, {"script": "demand"}],
        ), patch.object(
            state_simulation_main.workflow_assembly,
            "prepare_post_demand_workflow_for_state",
            return_value=[{"script": "emissions", "parameters": {}}],
        ):
            state_simulation_main.simulate_all_states(self.config, timeline_name="timeline")

        self.assertEqual(len(simulate_calls), 2)


if __name__ == "__main__":
    unittest.main()
