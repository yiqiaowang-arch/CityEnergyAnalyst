"""
===========================
Analytical energy demand model algorithm
===========================

"""
from __future__ import division
import multiprocessing as mp
import pandas as pd
import functions as f
import globalvar
import inputlocator
import maker as m
from cea.utils import epwreader
from cea import thermal_loads
from cea.thermal_loads import BuildingProperties

__author__ = "Jimeno A. Fonseca"
__copyright__ = "Copyright 2015, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Jimeno A. Fonseca", "Daren Thomas", "Gabriel Happle"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "thomas@arch.ethz.ch"
__status__ = "Production"


reload(f)
reload(globalvar)


def demand_calculation(locator, weather_path, gv):
    """
    Algorithm to calculate the hourly demand of energy services in buildings
    using the integrated model of Fonseca et al. 2015. Applied energy.
    (http://dx.doi.org/10.1016/j.apenergy.2014.12.068)

    PARAMETERS
    ----------
    :param locator: An InputLocator to locate input files
    :type locator: inputlocator.InputLocator

    :param weather_path: A path to the EnergyPlus weather data file (.epw)
    :type weather_path: str

    :param gv: A GlobalVariable (context) instance
    :type gv: globalvar.GlobalVariable


    RETURNS
    -------

    :returns: None
    :rtype: NoneType


    INPUT / OUTPUT FILES
    --------------------

    - get_demand_results_folder: C:\reference-case\baseline\outputs\data\demand
    - get_temporary_folder: c:\users\darthoma\appdata\local\temp
    - get_temporary_file: c:\users\darthoma\appdata\local\temp\B153767T.csv (* for each building)
    - get_total_demand: C:\reference-case\baseline\outputs\data\demand\Total_demand.csv


    SIDE EFFECTS
    ------------

    Produces a demand file per building and a total demand file for the whole zone of interest.

    B153767T.csv: csv file for every building with hourly demand data
    Total_demand.csv: csv file of yearly demand data per buidling.
    """

    # starting date
    date = pd.date_range(gv.date_start, periods=8760, freq='H')

    # weather model
    weather_data = epwreader.epw_reader(weather_path)[['drybulb_C', 'relhum_percent', 'windspd_ms']]

    # building properties model
    building_properties = BuildingProperties(locator, gv)

    # schedules model
    list_uses = list(building_properties._prop_occupancy.drop('PFloor', axis=1).columns)
    schedules = m.schedule_maker(date, locator, list_uses)
    schedules_dict = {'list_uses': list_uses, 'schedules': schedules}

    # demand model
    num_buildings = len(building_properties)
    if gv.multiprocessing:
        thermal_loads_all_buildings_multiprocessing(building_properties, date, gv, locator, num_buildings,
                                                    schedules_dict,
                                                    weather_data)
    else:
        thermal_loads_all_buildings(building_properties, date, gv, locator, num_buildings, schedules_dict,
                                    weather_data)
    write_totals_csv(building_properties, locator)
    gv.log('done')


def write_totals_csv(building_properties, locator):
    """read in the temporary results files and append them to the Totals.csv file."""
    counter = 0
    for name in building_properties.list_building_names():
        temporary_file = locator.get_temporary_file('%(name)sT.csv' % locals())
        if counter == 0:
            df = pd.read_csv(temporary_file)
            counter += 1
        else:
            df2 = pd.read_csv(temporary_file)
            df = df.append(df2, ignore_index=True)
    df.to_csv(locator.get_total_demand(), index=False, float_format='%.2f')


def thermal_loads_all_buildings(building_properties, date, gv, locator, num_buildings, usage_schedules,
                                weather_data):
    for i, building in enumerate(building_properties.list_building_names()):
        bpr = building_properties[building]
        thermal_loads.calc_thermal_loads_new_ventilation(
            building, bpr, weather_data, usage_schedules, date, gv,
            locator.get_demand_results_folder(),
            locator.get_temporary_folder())
        gv.log('Building No. %(bno)i completed out of %(num_buildings)i', bno=i + 1, num_buildings=num_buildings)


def thermal_loads_all_buildings_multiprocessing(building_properties, date, gv, locator, num_buildings, usage_schedules,
                                                weather_data):
    pool = mp.Pool()
    gv.log(mp.cpu_count())
    joblist = []
    for building in building_properties.list_building_names():
        bpr = building_properties[building]
        job = pool.apply_async(thermal_loads.calc_thermal_loads_new_ventilation,
                               [building, bpr, weather_data, usage_schedules, date, gv,
                                locator.get_demand_results_folder(),
                                locator.get_temporary_folder()])
        joblist.append(job)
    for i, job in enumerate(joblist):
        job.get(60)
        gv.log('Building No. %(bno)i completed out of %(num_buildings)i', bno=i + 1, num_buildings=num_buildings)


def test_demand():
        locator = inputlocator.InputLocator(scenario_path=r'C:\reference-case\baseline')
        # for the interface, the user should pick a file out of of those in ...DB/Weather/...
        weather_path = locator.get_default_weather()
        gv = globalvar.GlobalVariables()
        demand_calculation(locator=locator, weather_path=weather_path, gv=gv)
        print "test_demand() succeeded"


if __name__ == '__main__':
    test_demand()

