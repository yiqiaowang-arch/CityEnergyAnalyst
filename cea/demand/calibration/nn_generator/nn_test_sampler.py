# coding=utf-8
"""
'nn_test_sampler.py' script generates one random sample for the entire case-study,
to be called sequentially if necessary
"""

__author__ = "Fazel Khayatian"
__copyright__ = "Copyright 2017, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Fazel Khayatian"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"

from cea.demand.calibration.latin_sampler import latin_sampler
from cea.demand.demand_main import properties_and_schedule
from cea.demand.calibration.bayesian_calibrator.calibration_sampling import apply_sample_parameters
from cea.demand import demand_main
import pickle
import cea
import json
#import h5py
import os
import numpy as np
import pandas as pd
from cea.demand.calibration.nn_generator.nn_settings import number_samples, random_variables,\
    target_parameters, boolean_vars
from cea.demand.calibration.nn_generator.input_prepare import input_prepare_main

def sampling_single(locator, random_variables, target_parameters, list_building_names, weather_path, gv):
    size_city = np.shape(list_building_names)
    size_city=size_city[0]

    bld_counter = 0
    # create list of samples with a LHC sampler and save to disk (*.csv)
    samples, pdf_list = latin_sampler(locator, size_city, random_variables)
    for building_name in (list_building_names):
        np.save(locator.get_calibration_samples(building_name), samples)
        problem = {'variables': random_variables,
                   'building_load': target_parameters, 'probabiltiy_vars': pdf_list}
        pickle.dump(problem, file(locator.get_calibration_problem(building_name), 'w'))
        sample = np.asarray(zip(random_variables, samples[bld_counter, :]))
        apply_sample_parameters(locator, sample)
        bld_counter = bld_counter + 1
    # read the saved *.csv file and replace Boolean with logical (True/False)
    overwritten = pd.read_csv(locator.get_building_overrides())
    bld_counter = 0
    for building_name in (list_building_names):
        sample = np.asarray(zip(random_variables, samples[bld_counter, :]))
        for boolean_mask in (boolean_vars):
            indices = np.where(sample == boolean_mask)

            if sample[indices[0], 1] == '0.0':
                sample[indices[0], 1] = 'False'
            else:
                sample[indices[0], 1] = 'True'

        overwritten.loc[overwritten.Name == building_name, random_variables] = sample[:, 1]
        bld_counter = bld_counter + 1

    # write to csv format
    overwritten.to_csv(locator.get_building_overrides())

    #   run cea demand
    demand_main.demand_calculation(locator, weather_path, gv)
    #   prepare the inputs for feeding into the neural network
    urban_input_matrix, urban_taget_matrix = input_prepare_main(list_building_names, locator, target_parameters, gv)
    #   drop half the inputs and targets to avoid overfitting and save RAM / Disk space
    urban_input_matrix, urban_taget_matrix = input_dropout(urban_input_matrix, urban_taget_matrix)
    #   get the pathfor saving the files
    nn_inout_path = locator.get_nn_inout_folder()
    #   save inputs with sequential naming
    file_path_inputs = os.path.join(nn_inout_path, "input%(i)s.csv" % locals())
    data_file_inputs = pd.DataFrame(urban_input_matrix)
    data_file_inputs.to_csv(file_path_inputs, header=False, index=False)
    #   save inputs with sequential naming
    file_path_targets = os.path.join(nn_inout_path, "target%(i)s.csv" % locals())
    data_file_targets = pd.DataFrame(urban_taget_matrix)
    data_file_targets.to_csv(file_path_targets, header=False, index=False)

    return urban_input_matrix, urban_taget_matrix


def run_as_script():
    gv = cea.globalvar.GlobalVariables()
    scenario_path = gv.scenario_reference
    locator = cea.inputlocator.InputLocator(scenario_path=scenario_path)
    weather_path = locator.get_default_weather()
    building_properties, schedules_dict, date = properties_and_schedule(gv, locator)
    list_building_names = building_properties.list_building_names()
    urban_input_matrix, urban_taget_matrix=sampling_single(locator, random_variables, target_parameters, list_building_names, weather_path, gv)




if __name__ == '__main__':
    run_as_script()
