"""
This file cleans and fills in the gap of metered data based on weather conditions
"""
from __future__ import division

import os
import cea.config
import cea.inputlocator

__author__ = "Jimeno A. Fonseca"
__copyright__ = "Copyright 2017, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Clayton Miller, Jimeno A. Fonseca"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


def data_cleaner(locator, scenario_path):

    pass


def main(config):
    assert os.path.exists(config.scenario), 'Scenario not found: %s' % config.scenario
    locator = cea.inputlocator.InputLocator(config.scenario)

    # print out all configuration variables used by this script
    print("Running calculation with scenario = %s" % config.scenario)
    print("Running calculation with archetypes = %s" % config.data_helper.archetypes)
    print("Running calculation with the next buildings = %s" % config.data_cleaner.buildings)
    print("Running calculation between the next period = %s - %s" % config.data_cleaner.start_date,config.data_cleaner.end_date)

    data_cleaner(locator = locator, scenario_path = config.scenario, buildings = config.datacleaner.buildings)

if __name__ == '__main__':
    main(cea.config.Configuration())
