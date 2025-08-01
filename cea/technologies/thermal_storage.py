"""
thermal storage
"""




import pandas as pd
from math import log
from cea.analysis.costs.equations import calc_capex_annualized
__author__ = "Thuy-An Nguyen"
__copyright__ = "Copyright 2015, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Thuy-An Nguyen", "Tim Vollrath", "Jimeno A. Fonseca"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "cea@arch.ethz.ch"
__status__ = "Production"


# investment and maintenance costs

def calc_Cinv_storage(V_tank_m3, locator, technology_type):
    """
    calculate the annualized investment cost of a thermal storage tank

    :param V_tank_m3: storage tank volume
    :type V_tank_m3: float

    :returns InvCa:

    """
    if V_tank_m3 > 0:
        storage_cost_data = pd.read_csv(locator.get_db4_components_conversion_conversion_technology_csv("THERMAL_ENERGY_STORAGES"))
        storage_cost_data = storage_cost_data[storage_cost_data['code'] == technology_type]

        # if the Q_design is below the lowest capacity available for the technology, then it is replaced by the least
        # capacity for the corresponding technology from the database
        if V_tank_m3 < storage_cost_data.iloc[0]['cap_min']:
            V_tank_m3 = storage_cost_data.iloc[0]['cap_min']

        storage_cost_data = storage_cost_data[
            (storage_cost_data['cap_min'] <= V_tank_m3) & (storage_cost_data['cap_max'] > V_tank_m3)]

        Inv_a = storage_cost_data.iloc[0]['a']
        Inv_b = storage_cost_data.iloc[0]['b']
        Inv_c = storage_cost_data.iloc[0]['c']
        Inv_d = storage_cost_data.iloc[0]['d']
        Inv_e = storage_cost_data.iloc[0]['e']
        Inv_IR = storage_cost_data.iloc[0]['IR_%']
        Inv_LT = storage_cost_data.iloc[0]['LT_yr']
        Inv_OM = storage_cost_data.iloc[0]['O&M_%'] / 100

        Capex_total_USD = Inv_a + Inv_b * (V_tank_m3) ** Inv_c + (Inv_d + Inv_e * V_tank_m3) * log(V_tank_m3)
        Capex_a_storage_USD = calc_capex_annualized(Capex_total_USD, Inv_IR, Inv_LT)
        Opex_fixed_storage_USD = Capex_total_USD * Inv_OM
    else:
        Capex_a_storage_USD = 0.0
        Opex_fixed_storage_USD = 0.0
        Capex_total_USD = 0.0

    return Capex_a_storage_USD, Opex_fixed_storage_USD, Capex_total_USD
