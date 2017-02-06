"""
===========================
Primary energy and CO2 emissions model algorithm for mobility
===========================

M. Mosteiro Romero  script development          31.08.16

"""

from __future__ import division

import os

import pandas as pd
from geopandas import GeoDataFrame as gpdf

from cea import inputlocator

reload(inputlocator)

__author__ = "Martin Mosteiro Romero"
__copyright__ = "Copyright 2016, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Martin Mosteiro Romero"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"

def lca_mobility(locator):
    """
    Algorithm to calculate the primary energy and CO2 emissions for mobility
    in the area in order to compare with the 2000 Watt society benchmark
    based on the present day values calculated for the 2000 Watt society target.

    The current values for the Swiss case for each type of occupancy are based on the following sources:
        -   Swiss Society of Engineers and Architects (SIA), "SIA Efficiency Path 2040" (2011):
            'MULTI_RES', 'SINGLE_RES', 'SCHOOL', 'OFFICE'
        -   Kellenberger, D. et al. "Arealentwicklung fur die 2000-Watt gesellschaft: Leitfaden und Fallbeispiele" (2012):
            'HOTEL', 'RETAIL', 'FOODSTORE', 'RESTAURANT'; 'HOSPITAL' assumed same as for hotel
        -   Assumed from Fonseca, J. et al. "Assessing the environmental impact of future urban developments at neighborhood scale" (2015):
            'INDUSTRY'
    Due to lack of sources, the following target values are, as of yet, undefined:
        'GYM', 'SWIMMING', 'SERVERROOM', 'COOLROOM'
    
    Parameters
    ----------
    :param InputLocator locator: an InputLocator instance set to the scenario to work on

    Side effects
    ------------
    The following file is created by this script:
    - total_LCA_mobility:.csv
        csv file of yearly nonr-renewable primary energy demand and emissions due to mobility for each building
    """

    # local files
    demand = pd.read_csv(locator.get_total_demand())
    prop_occupancy = gpdf.from_file(locator.get_building_occupancy()).drop('geometry', axis=1)#.set_index('Name')
    factors_mobility = pd.read_excel(locator.get_data_benchmark_today(), sheetname='mobility').drop('Description', axis=1)

    # calculate total_LCA_mobility:.csv
    occupancy_type = factors_mobility['code']
    non_renewable_energy = factors_mobility['PEN']
    emissions = factors_mobility['CO2']

    mobility = prop_occupancy.merge(demand,on='Name')
    fields_to_plot = ['Name', 'GFA_m2', 'M_nre_pen_GJ', 'M_nre_pen_MJm2', 'M_ghg_ton', 'M_ghg_kgm2']
    mobility[fields_to_plot[3]] = 0
    mobility[fields_to_plot[5]] = 0
    for i in range(len(occupancy_type)):
        mobility[fields_to_plot[3]] += mobility[occupancy_type[i]] * non_renewable_energy[i]
        mobility[fields_to_plot[5]] += mobility[occupancy_type[i]] * emissions[i]
    mobility[fields_to_plot[2]] = mobility['GFA_m2'] * mobility[fields_to_plot[3]] / 1000
    mobility[fields_to_plot[4]] = mobility['GFA_m2'] * mobility[fields_to_plot[5]] / 1000

    mobility[fields_to_plot].to_csv(locator.get_lca_mobility(), index=False, float_format='%.2f')

def run_as_script(scenario_path=None):
    import cea.globalvar
    gv = cea.globalvar.GlobalVariables()
    if not scenario_path:
        scenario_path = gv.scenario_reference
    locator = cea.inputlocator.InputLocator(scenario_path=scenario_path)
    lca_mobility(locator=locator)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--scenario', help='Path to the scenario folder')
    args = parser.parse_args()
    run_as_script(scenario_path=args.scenario)
