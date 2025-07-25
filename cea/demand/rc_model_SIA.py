# -*- coding: utf-8 -*-





import numpy as np
from cea.demand import constants

__author__ = "Gabriel Happle"
__copyright__ = "Copyright 2016, Architecture and Building Systems - ETH Zurich"
__credits__ = ["Gabriel Happle"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Daren Thomas"
__email__ = "thomas@arch.ethz.ch"
__status__ = "Production"

'''
RC model calculations according to sia 2044

Merkblatt 2044 Kilimatisierte Gebauede - Standard-Berechnungsverfahren fuer den Leistungs-und Energiebedarf
'''

# import constants
T_WARNING_LOW = constants.T_WARNING_LOW
T_WARNING_HIGH = constants.T_WARNING_HIGH

# TODO: documentation

# SIA 2044 constants
# section 2.1.3 in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011
h_cv_i = 2.5  # (W/m2K) heat transfer coefficient for convection on internal area
h_r_i = 5.5  # (W/m2K) heat trasnfer coefficient for long wave radiation
h_ic = 9.1  # (W/m2K) total heat transfer coefficient on internal area, = h_cv_i + 1.2*h_r_i
# section 2.1.4 in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011
f_sa = 0.1  # (-) share of solar radiation that turns to convective heat
f_r_l = 0.7  # (-)
f_r_p = 0.5  # (-)
f_r_a = 0.2  # (-)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1.3
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def calc_h_mc(a_m):
    """
    :param a_m: see ``bpr.rc_model['Am']``
    :return:
    """

    # get properties from bpr # TODO: to be addressed in issue #443
    # a_m = bpr.rc_model['Am']

    # (7) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    h_mc = h_ic * a_m

    return h_mc


def calc_h_ac(a_t: float) -> float:
    """
    Source: (10) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    :param a_t: area of all surfaces facing the building zone in [m2], equivalent to `bpr.rc_model['Atot']`. 
                See documentation for `calc_prop_rc_model`.
    :type a_t: float
    :return: `h_ac`: heat loss factor from air node to central node, used in 5R1C model.
    :rtype: float
    """

    # get properties from bpr # TODO: to be addressed in issue #443

    # 

    h_ac = a_t / (1 / h_cv_i - 1 / h_ic)

    return h_ac


def calc_h_op_m(Htr_op):

    # work around # TODO: to be addressed in issue #443
    # get h_op from ISO model (with basement factor)
    h_op_m = Htr_op
    # TODO: This formula should be adjusted to be compatible with SIA2044

    # (9) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # h_op_m = a_j_m * u_j  # summation
    # This formula in the future should take specific properties of the location of the building into account.

    return h_op_m


def calc_h_em(h_op_m: float, h_mc: float) -> float:
    """
    Source: (11) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    calculates `h_em` for the 5R1C model.

    :param h_op_m: total heat loss factor: `U_op * A_op`, including conduction, convection and radiation.
    :type h_op_m: float
    :param h_mc: heat loss factor from thermal mass to internal air due to convection and radiation.
    :type h_mc: float
    :return: heat loss factor due to only conduction between mass node and outside.
    :rtype: float
    """
    if h_op_m > 0:
        h_em = 1.0 / (1.0 / h_op_m - 1.0 / h_mc)
    else:
        # h_op_m = 0, no heat transfer from mass to outside air. 
        # Therefore h_em (part of heat transfer from mass to outside air) should also be 0.
        h_em = 0

    return h_em


def calc_h_j_em():

    # TODO: to be addressed in issue #443

    # (11) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    #h_j_em = (h_em * a_j_m * u_j) / h_op_m

    # This formula in the future should take specific properties of the location of the building into account.

    return None


def calc_h_ec(Htr_w: float) -> float:
    """
    Source: (12) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011.\\
    Calculates transmission heat loss factor of all components with negligible thermal mass. 
    This category can include window, linear thermal bridge and point thermal bridge.\\
    
    Current simplified formula: `h_ec = Htr_w`, where Htr_w is the only input; \\
    Or conceptually, `h_ec = a_j_l * u_j`, 
    where a_j_l is the area of j-th light component, u_j is its U-value.

    :param Htr_w: transmission heat loss factor of windows.
    :type Htr_w: float
    :return: transmission heat loss factor of light component `h_ec`, used in the 5R1C model.
    :rtype: float
    """
    # This formula in the future should take specific properties of the location of the building into account.
    # TODO: can incorporate point or linear thermal bridges

    h_ec = Htr_w  # h_ec is Htr_w of ISO13790 RC model

    return h_ec


def calc_h_ea(m_ve_mech: float, m_ve_window: float, m_ve_inf_simple: float) -> float:
    """
    source: (13) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011.\\
    Adapted for mass flows instead of volume flows

    :param m_ve_mech: mass flow rate for mechanical ventilation in [kg/s]
    :type m_ve_mech: float
    :param m_ve_window: mass flow rate for window ventilation in [kg/s]
    :type m_ve_window: float
    :param m_ve_inf_simple: mass flow rate for infiltration in [kg/s]
    :type m_ve_inf_simple: float
    :return: `h_ea`: heat loss factor from outside to the air node, used in 5R1C model.
    :rtype: float
    """
    cp = 1.005 / 3.6  # (Wh/kg/K)
    # TODO: check units of air flow

    # get values
    m_v_sys = m_ve_mech * 3600  # mass flow rate mechanical ventilation
    m_v_w = m_ve_window * 3600  # mass flow rate window ventilation
    m_v_inf = m_ve_inf_simple * 3600  # mass flow rate infiltration

    h_ea = (m_v_sys + m_v_w + m_v_inf) * cp

    return h_ea


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1.4
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def calc_phi_a(phi_hc_cv, phi_i_l, phi_i_a, phi_i_p, I_sol):

    # (14) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # Gabriel Happle 01.12.2016

    # get internal loads
    phi_i_l = phi_i_l
    phi_i_a = phi_i_a
    phi_i_p = phi_i_p

    # solar gains
    phi_s = I_sol  # solar gains

    # standard assumptions
    #f_sa = 0.1
    #f_r_l = 0.7
    #f_r_p = 0.5
    #f_r_a = 0.2

    phi_a = f_sa * phi_s + (1 - f_r_l) * phi_i_l + (1 - f_r_p) * phi_i_p +(1 - f_r_a) * phi_i_a + phi_hc_cv

    return phi_a


def calc_phi_c(phi_hc_r, phi_i_l, phi_i_a, phi_i_p, I_sol, f_ic, f_sc):

    # (15) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # Gabriel Happle 01.12.2016

    # get internal loads
    phi_i_l = phi_i_l
    phi_i_a = phi_i_a
    phi_i_p = phi_i_p

    # solar gains
    phi_s = I_sol  # solar gains

    # call functions for factor
    f_ic = f_ic
    f_sc = f_sc

    # standard assumptions
    #f_sa = 0.1
    #f_r_l = 0.7
    #f_r_p = 0.5
    #f_r_a = 0.2

    phi_c = f_ic * (f_r_l * phi_i_l + f_r_p * phi_i_p + f_r_a * phi_i_a + phi_hc_r) + (1 - f_sa) * f_sc * phi_s

    return phi_c


def calc_phi_i_p(Qs: float) -> float:  # _Wp, people):
    """internal gain due to occupancy.

    :param Qs: total internal load from people in [W]
    :type Qs: float
    :return: total internal load from people in [W]
    :rtype: float
    """

    # # internal gains from people
    # phi_i_p = people * Qs_Wp
    return Qs # phi_i_p


def calc_phi_i_a(Eaf: float, Epro: float) -> float:
    """internal gains from appliances and industrial processes.

    :param Eaf: appliance load in W
    :type Eaf: float
    :param Epro: industrial process load in W
    :type Epro: float
    :return: `phi_i_a`, total internal gain due to appliances and industrial processes in [W]
    :rtype: float
    """
    # internal gains from appliances, factor of 0.9 taken from old method calc_Qgain_sen()
    # TODO make function and dynamic, check factor
    phi_i_a = 0.9 * (Eaf + Epro)
    return phi_i_a


def calc_phi_i_l(Elf: float) -> float:
    """internal gains from lighting.

    :param Elf: lighting load in W
    :type Elf: float
    :return: `phi_i_l`, total internal gain due to lighting load in [W].
    :rtype: float
    """
    # factor of 0.9 taken from old method calc_Qgain_sen()
    # TODO make function and dynamic, check factor
    phi_i_l = 0.9 * Elf
    return phi_i_l


def calc_phi_m(phi_hc_r, phi_i_l: float, phi_i_a, phi_i_p, I_sol, f_im, f_sm):

    # (16) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # Gabriel Happle 01.12.2016

    # get internal loads
    phi_i_l = phi_i_l
    phi_i_a = phi_i_a
    phi_i_p = phi_i_p

    # solar gains
    phi_s = I_sol  # solar gains

    # call functions for factors
    f_im = f_im
    f_sm = f_sm

    # standard assumption
    #f_sa = 0.1
    #f_r_l = 0.7
    #f_r_p = 0.5
    #f_r_a = 0.2
    phi_m = f_im * (f_r_l * phi_i_l + f_r_p * phi_i_p + f_r_a * phi_i_a + phi_hc_r) + (1 - f_sa) * f_sm * phi_s

    return phi_m


def calc_f_ic(a_t: float, a_m: float, h_ec: float) -> float:
    """
    Source: (17) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011\\
    Gabriel Happle 01.12.2016

    :param a_t: area of all surfaces facing the building zone in [m2], equivalent to `bpr.rc_model['Atot']`. 
                See documentation for `calc_prop_rc_model`.
    :type a_t: float
    :param a_m: effective mass area in [m2], see `bpr.rc_model['Am']` and documentation for `calc_prop_rc_model`.
    :type a_m: float
    :param h_ec: transmission heat loss factor of light component, see function `calc_h_ec`.
    :type h_ec: float
    :return: fraction of internal radiant heat gain (occupants, equipments, etc.) distributed to central node.
    :rtype: float
    """
    f_ic = (a_t - a_m - h_ec / h_ic) / a_t

    return f_ic


def calc_f_sc(a_t: float, a_m: float, a_w: float, h_ec: float) -> float:
    """
    Source: (18) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011\\
    Gabriel Happle 01.12.2016

    :param a_t: area of all surfaces facing the building zone in [m2], equivalent to `bpr.rc_model['Atot']`. 
                See documentation for `calc_prop_rc_model`.
    :type a_t: float
    :param a_m: effective mass area in [m2], see `bpr.rc_model['Am']` and documentation for `calc_prop_rc_model`.
    :type a_m: float
    :param a_w: area of all windows. See `bpr.rc_model['Aw']` and documentation for `calc_prop_rc_model`.
    :type a_w: float
    :param h_ec: transmission heat loss factor of light component, see function `calc_h_ec`.
    :type h_ec: float
    :return: fraction of solar gain that feeds to central node in 5R1C model.
    :rtype: float
    """
    f_sc = (a_t-a_m-a_w-h_ec/h_ic) / (a_t - a_w)

    return f_sc


def calc_f_im(a_t, a_m):
    """

    :param a_t: see ``bpr.rc_model['Atot']``
    :param a_m: see ``bpr.rc_model['Am']``
    :return:
    """

    # (19) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # Gabriel Happle 01.12.2016

    f_im = a_m / a_t

    return f_im


def calc_f_sm(a_t, a_m, a_w):
    """
    :param a_t: bpr.rc_model['Atot']
    :param a_m: bpr.rc_model['Am']
    :param a_w: bpr.rc_model['Aw']
    :return:
    """

    # (20) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # Gabriel Happle 01.12.2016

    f_sm = a_m / (a_t - a_w)

    return f_sm

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1.5
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def calc_theta_ea(m_ve_mech, m_ve_window, m_ve_inf_simple, theta_ve_mech, T_ext):

    # get values
    m_v_sys = m_ve_mech  # mass flow rate mechanical ventilation
    m_v_w = m_ve_window  # mass flow rate window ventilation
    m_v_inf = m_ve_inf_simple  # mass flow rate infiltration
    theta_v_sys = theta_ve_mech  # supply air temperature of mechanical ventilation (i.e. after HEX)
    theta_e = T_ext  # outdoor air temperature

    # (21) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # adjusted for mass flows instead of volume flows and simplified by (rho*cp)/(rho*cp) = 1
    # Gabriel Happle 01.12.2016

    theta_ea = (m_v_sys * theta_v_sys + (m_v_w + m_v_inf) * theta_e) / (m_v_sys + m_v_w + m_v_inf)

    return theta_ea


def calc_theta_ec(T_ext):

    # WORKAROUND
    theta_ec = T_ext  # TODO: adjust to actual calculation

    # (22) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    # theta_ec = a_j_l * u_j * theta_e_j / h_ec

    # This formula in the future should take specific properties of the location of the building into account.

    # TODO: theta_e_j is dependent on adjacent space to surface (outdoor, adiabatic, ground, etc.)

    return theta_ec


def calc_theta_em(T_ext):

    # WORKAROUND
    theta_em = T_ext  # TODO: adjust to actual calculation

    # (23) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    # theta_em = h_j_em * theta_e_j / h_em

    # This formula in the future should take specific properties of the location of the building into account.

    # TODO: theta_e_j is dependent on adjacent space to surface (outdoor, adiabatic, ground, etc.)

    return theta_em


def calc_theta_e_star():

    # TODO: To be addressed in issue #446
    pass

    # (24) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    #

    # standard values for calculation
    # h_e_load_and_temp = 13.5  # (W/m2K) for load and temperature calculation
    # h_e_energy = 23  # (W/m2K) for energy demand calculation
    #
    # f_r_roof = 1  # (-)
    # f_r_wall = 0.5  # (-)
    # h_r = 5.5  # (-)
    # delta_t_er = 11  # (K)

    # if is_roof(surface):
        #f_r = f_r_roof
    #elif is_wall(surface):
        #f_r = f_r_wall
    #else:
        #raise()

    #if is_energy_calculation(calculation):
       # h_e = h_e_energy
    #elif is_load_or_temp_calculation(calculation):
       # h_e = h_e_load_and_temp
    #else:
        #raise()


    #theta_e_star = theta_e + (alpha_s * i_s_i) / h_e - (f_r * h_r * epsilon_0 * delta_t_er) / h_e

    #return theta_e_star

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.1.6
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def calc_theta_m_t(phi_m_tot, theta_m_t_1, h_em, h_3, c_m):
    # (25) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    theta_m_t = (theta_m_t_1 * (c_m - 0.5 * (h_3 + h_em)) + phi_m_tot) / (c_m + 0.5 * (h_3 + h_em))

    return theta_m_t


def calc_h_1(h_ea: float, h_ac: float) -> float:
    """
    source: (26) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    Calculates the series heat transmittance from central node to outdoor air via convection.

    :param h_ea: heat loss factor between outdoor air and air node.
    :type h_ea: float
    :param h_ac: heat loss factor between air node and central node (e.g., via surface convection). 
    :type h_ac: float
    :return: augmented heat loss factor from central node to outdoor air.
    :rtype: float
    """

    # get values
    h_ea = h_ea
    h_ac = h_ac

    h_1 = 1 / (1 / h_ea + 1 / h_ac)

    return h_1


def calc_h_2(h_1: float, h_ec: float) -> float:
    """
    Source: (27) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    Calculates total heat loss factor from central node to outdoor air via both convection and conduction.

    :param h_1: heat loss factor from central node via air node to outside air (e.g., through convection).
    :type h_1: float
    :param h_ec: heat loss factor from central node to outside air (e.g., through conduction).
    :type h_ec: float
    :return: total heat loss factor from central node to outdoor air via both paths.
    :rtype: float
    """

    h_2 = h_1 + h_ec

    return h_2


def calc_h_3(h_2: float, h_mc: float) -> float:
    """
    Source: (28) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    Calculates the heat loss factor from mass node via central node (via conduction) 
    to the outside (via conduction and convection).

    :param h_2: total heat loss factor from central node to outside air via convection and conduction
    :type h_2: float
    :param h_mc: heat loss factor from mass node to central node via conduction.
    :type h_mc: float
    :return: augmented heat loss factor from mass node via central node to outside air.
    :rtype: float
    """

    h_3 = 1.0 / (1.0 / h_2 + 1.0 / h_mc)
    return h_3


def calc_phi_m_tot(phi_m, phi_a, phi_c, theta_ea, theta_em, theta_ec, h_1, h_2, h_3, h_ec, h_ea, h_em):
    # (29) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    phi_m_tot = phi_m + h_em * theta_em + (h_3 * (phi_c + h_ec * theta_ec + h_1 * (phi_a / h_ea + theta_ea))) / h_2
    return phi_m_tot


def calc_theta_m(theta_m_t, theta_m_t_1):
    # (30) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    theta_m = (theta_m_t + theta_m_t_1) / 2
    return theta_m


def calc_theta_c(phi_a, phi_c, theta_ea, theta_ec, theta_m, h_1, h_mc, h_ec, h_ea):

    # get values
    h_mc = h_mc
    h_ec = h_ec
    h_ea = h_ea

    # (31) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    theta_c = (h_mc * theta_m + phi_c + h_ec * theta_ec + h_1 * (phi_a / h_ea + theta_ea)) / (h_mc + h_ec + h_1)

    return theta_c


def calc_T_int(phi_a, theta_ea, theta_c, h_ac, h_ea):
    # (32) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    T_int = (h_ac * theta_c + h_ea * theta_ea + phi_a) / (h_ac + h_ea)
    return T_int


def calc_theta_o(T_int, theta_c):
    # (33) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    theta_o = T_int * 0.31 + theta_c * 0.69
    return theta_o

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.2.7
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def calc_phi_hc_cv(phi_hc, f_hc_cv):

    # (58) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    phi_hc_cv = f_hc_cv * phi_hc

    return phi_hc_cv


def calc_phi_hc_r(phi_hc, f_hc_cv):

    # (59) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    phi_hc_r = (1 - f_hc_cv) * phi_hc

    return phi_hc_r


def calc_theta_tabs_su():

    # TODO: to be addressed in issue #444

    # (60) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # theta_tabs_su = theta_tabs_su_max - (theta_tabs_su_max - theta_tabs_su_min) * (theta_e -  theta_e_min)/(theta_e_max - theta_e_min)

    return None


def calc_phi_tabs():

    # TODO: to be addressed in issue #444

    # (61) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011
    # phi_tabs = h_tabs * (theta_tabs_su - theta_m_t_1)

    return None


def calc_h_tabs():

    # TODO: to be addressed in issue #444

    # (62) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    # typical values
    # a_tabs = 0.8 * a_ngf
    # r_tabs = 0.08  # (m2K/W)

    # h_tabs = a_tabs / r_tabs

    return None


def calc_phi_m_tot_tabs():

    # TODO: to be addressed in issue #444

    # (63) in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    # phi_m_tot = calc_phi_m_tot() + phi_tabs

    return None


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2.3.2
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def calc_rc_model_temperatures_no_heating_cooling(bpr, tsd, t, config):
    """
    Calculates R-C-Model temperatures are calculated with zero heating/cooling power according to SIA 2044 procedure.

    :py:func: `cea.demand.rc_model_SIA.calc_rc_model_temperatures_no_heating_cooling`

    Author: Gabriel Happle
    Date: FEB 2017

    :param bpr: Building Properties
    :type bpr: BuildingPropertiesRow
    :param tsd: Time series data of building
    :type tsd: dict
    :param t: time step / hour of the year
    :type t: int
    :return: R-C-Model node temperatures
    :rtype: dict
    """

    # no heating or cooling
    phi_hc_cv = 0.0
    phi_hc_r = 0.0

    # calculate r-c-model node temperatures
    temperatures = calc_rc_model_temperatures(phi_hc_cv, phi_hc_r, bpr, tsd, t, config)
    rc_model_temp = temperatures

    return rc_model_temp


def calc_rc_model_temperatures(phi_hc_cv, phi_hc_r, bpr, tsd, t, config):
    # calculate node temperatures of RC model
    theta_m_t_1 = tsd['theta_m'][t - 1]
    if np.isnan(theta_m_t_1):
        theta_m_t_1 = tsd['T_ext'][t - 1]

    # copy data required for calculation from `tsd` for this timestep
    m_ve_mech = tsd['m_ve_mech'][t]
    m_ve_window = tsd['m_ve_window'][t]
    m_ve_inf = tsd['m_ve_inf'][t]
    El = tsd['El'][t] * min(bpr.rc_model['Af']/bpr.rc_model['Aef'], 1.0) # account for a proportion of internal gains
    Ea = tsd['Ea'][t] * min(bpr.rc_model['Af']/bpr.rc_model['Aef'], 1.0)  # account for a proportion of internal gains
    Epro = tsd['Epro'][t]
    # account for a proportion of solar gains. This is very simplified for now.
    I_sol = tsd['I_sol_and_I_rad'][t] * np.sqrt(bpr.architecture.Hs_ag)
    T_ext = tsd['T_ext'][t]
    theta_ve_mech = tsd['theta_ve_mech'][t]

    # copy data from `bpr`
    Htr_op = bpr.rc_model['Htr_op']
    Htr_w = bpr.rc_model['Htr_w']
    Qs = tsd['Qs'][t]
    a_t = bpr.rc_model['Atot']
    a_m = bpr.rc_model['Am']
    a_w = bpr.rc_model['Awin_ag']
    c_m = bpr.rc_model['Cm'] / 3600  # (Wh/K) SIA 2044 unit is Wh/K, ISO unit is J/K

    T_int, theta_c, theta_m, theta_o, theta_ea, theta_ec, theta_em, h_ea, h_ec, h_em, h_op_m \
        = _calc_rc_model_temperatures(Ea, El, Epro, Htr_op, Htr_w, I_sol, Qs, T_ext, a_m, a_t, a_w, c_m, m_ve_inf,
                                                                     m_ve_mech, m_ve_window, phi_hc_cv,
                                                                     phi_hc_r, theta_m_t_1, theta_ve_mech)

    if T_WARNING_LOW > T_int or T_WARNING_LOW > theta_c or T_WARNING_LOW > theta_m \
            or T_int > T_WARNING_HIGH or theta_c > T_WARNING_HIGH or theta_m > T_WARNING_HIGH:
        if config.demand.overheating_warning:
            raise Exception("Temperature in RC-Model of building {} out of bounds! First occurred at timestep = {}. "
                            "The results were Tint = {}, theta_c = {}, theta_m = {}.\n"
                            "If it is an expected behavior, consider turning off over-heating warning in the "
                            "advanced parameters to continue the simulation.\n"
                            "If it is not expected, check building geometry and internal loads.\n" 
                            "Building might be too small in size or architecture parameter Hs_ag = {} might be too "
                            "small for this geometry. Current bounds of range for RC-model temperatures are "
                            "between {} and {}.".format(bpr.name, t, round(T_int, 2), round(theta_c, 2),  round(theta_m, 2),
                                                         bpr.architecture.Hs_ag, T_WARNING_LOW, T_WARNING_HIGH))

    rc_model_temp = {'theta_m': theta_m, 'theta_c': theta_c, 'T_int': T_int, 'theta_o': theta_o, 'theta_ea': theta_ea,
                     'theta_ec': theta_ec, 'theta_em': theta_em, 'h_ea': h_ea, 'h_ec': h_ec, 'h_em': h_em,
                     'h_op_m': h_op_m}
    return rc_model_temp


def _calc_rc_model_temperatures(Eaf, Elf, Epro, Htr_op, Htr_w, I_sol, Qs, T_ext, a_m, a_t, a_w, c_m,
                                m_ve_inf_simple, m_ve_mech, m_ve_window, phi_hc_cv, phi_hc_r, theta_m_t_1,
                                theta_ve_mech):
    # numba_cc compatible calculation
    h_ec = calc_h_ec(Htr_w=Htr_w)
    h_ac = calc_h_ac(a_t)
    h_ea = calc_h_ea(m_ve_mech, m_ve_window, m_ve_inf_simple)
    f_sc = calc_f_sc(a_t, a_m, a_w, h_ec)
    f_ic = calc_f_ic(a_t, a_m, h_ec)
    h_op_m = calc_h_op_m(Htr_op=Htr_op)
    h_mc = calc_h_mc(a_m=a_m)
    h_em = calc_h_em(h_op_m, h_mc)
    f_im = calc_f_im(a_t=a_t, a_m=a_m)
    f_sm = calc_f_sm(a_t=a_t, a_m=a_m, a_w=a_w)
    phi_i_l = calc_phi_i_l(Elf=Elf)
    phi_i_a = calc_phi_i_a(Eaf=Eaf, Epro=Epro) # include processes
    phi_i_p = calc_phi_i_p(Qs=Qs)
    h_1 = calc_h_1(h_ea=h_ea, h_ac=h_ac)
    phi_a = calc_phi_a(phi_hc_cv, phi_i_l, phi_i_a, phi_i_p, I_sol)
    phi_m = calc_phi_m(phi_hc_r, phi_i_l, phi_i_a, phi_i_p, I_sol, f_im, f_sm)
    phi_c = calc_phi_c(phi_hc_r, phi_i_l, phi_i_a, phi_i_p, I_sol, f_ic, f_sc)
    theta_ea = calc_theta_ea(m_ve_mech, m_ve_window, m_ve_inf_simple, theta_ve_mech, T_ext)
    theta_em = calc_theta_em(T_ext=T_ext)
    theta_ec = calc_theta_ec(T_ext=T_ext)
    h_2 = calc_h_2(h_1=h_1, h_ec=h_ec)
    h_3 = calc_h_3(h_2, h_mc)
    phi_m_tot = calc_phi_m_tot(phi_m, phi_a, phi_c, theta_ea, theta_em, theta_ec, h_1, h_2, h_3, h_ec, h_ea, h_em)
    theta_m_t = calc_theta_m_t(phi_m_tot, theta_m_t_1, h_em, h_3, c_m)
    theta_m = calc_theta_m(theta_m_t, theta_m_t_1)
    theta_c = calc_theta_c(phi_a, phi_c, theta_ea, theta_ec, theta_m, h_1, h_mc, h_ec, h_ea)
    T_int = calc_T_int(phi_a=phi_a, theta_ea=theta_ea, theta_c=theta_c, h_ac=h_ac, h_ea=h_ea)
    theta_o = calc_theta_o(T_int=T_int, theta_c=theta_c)
    return T_int, theta_c, theta_m, theta_o, theta_ea, theta_ec, theta_em, h_ea, h_ec, h_em, h_op_m


def calc_rc_model_temperatures_heating(phi_hc, bpr, tsd, t, config):
    """
    This function executes the equations of SIA 2044 R-C-Building-Model to calculate the node temperatures for a given
    heating energy demand

    .. py:function:: `cea.demand.rc_model_SIA.lookup_f_hc_cv_heating`
    .. py:function:: `cea.demand.rc_model_SIA.calc_phi_hc_cv`
    .. py:function:: `cea.demand.rc_model_SIA.calc_phi_hc_r`
    .. py:function:: `cea.demand.rc_model_SIA.calc_rc_model_temperatures`

    Author: Gabriel Happle
    Date: FEB 2017

    :param phi_hc: Heating or cooling energy demand of building
    :type phi_hc: float
    :param bpr: Building Properties
    :type bpr: BuildingPropertiesRow
    :param tsd: Time series data of building
    :type tsd: dict
    :param t: time step / hour of the year
    :type t: int
    :return: R-C-Building-Model node temperatures
    :rtype: dict
    """

    # lookup convection factor for heating/cooling system
    f_hc_cv = lookup_f_hc_cv_heating(bpr)

    # convective and radiative fractions of heating system
    phi_hc_cv = calc_phi_hc_cv(phi_hc, f_hc_cv)
    phi_hc_r = calc_phi_hc_r(phi_hc, f_hc_cv)

    # calculating R-C-Model node temperatures
    rc_model_temp = calc_rc_model_temperatures(phi_hc_cv, phi_hc_r, bpr, tsd, t, config)

    return rc_model_temp


def calc_rc_model_temperatures_cooling(phi_hc, bpr, tsd, t, config):
    """
    This function executes the equations of SIA 2044 R-C-Building-Model to calculate the node temperatures for a given
    cooling energy demand

    :py:func: `cea.demand.rc_model_SIA.lookup_f_hc_cv_cooling`
    :py:func: `cea.demand.rc_model_SIA.calc_phi_hc_cv`
    :py:func: `cea.demand.rc_model_SIA.calc_phi_hc_r`
    :py:func: `cea.demand.rc_model_SIA.calc_rc_model_temperatures`

    Author: Gabriel Happle
    Date: FEB 2017

    :param phi_hc: Heating or cooling energy demand of building
    :type phi_hc: float
    :param bpr: Building Properties
    :type bpr: BuildingPropertiesRow
    :param tsd: Time series data of building
    :type tsd: dict
    :param t: time step / hour of the year
    :type t: int
    :return: R-C-Building-Model node temperatures
    :rtype: dict
    """

    # lookup convection factor for heating/cooling system
    f_hc_cv = lookup_f_hc_cv_cooling(bpr)

    # convective and radiative fractions of heating system
    phi_hc_cv = calc_phi_hc_cv(phi_hc, f_hc_cv)
    phi_hc_r = calc_phi_hc_r(phi_hc, f_hc_cv)

    # calculating R-C-Model node temperatures
    rc_model_temp = calc_rc_model_temperatures(phi_hc_cv, phi_hc_r, bpr, tsd, t, config)

    return rc_model_temp


def has_heating_demand(bpr, tsd, t, config):
    """
    This function checks whether the building R-C-Model has a heating demand according to the procedure in SIA 2044.
    R-C-Model temperatures are calculated with zero heating power and checked versus the set-point temperature.
    Function includes a temperature tolerance according to the precision of the result reporting.

    :py:func: `cea.demand.rc_model_SIA.calc_rc_model_temperatures_no_heating_cooling`

    Author: Gabriel Happle
    Date: FEB 2017

    :param bpr: Building Properties
    :type bpr: BuildingPropertiesRow
    :param tsd: Time series data of building
    :type tsd: dict
    :param t: time step / hour of the year
    :type t: int
    :return: True or False
    :rtype: bool
    """

    temp_tolerance = 0.001  # temperature tolerance of temperature sensor (°C),
    #  i.e. heating is turned on if temperature is temp_tolerance below the set point
    # tolerance is consistent with maximum temperature difference that can be reported with the current precision
    # of the demand *.csv file

    ta_hs_set = tsd['ta_hs_set'][t]
    if np.isnan(ta_hs_set):
        # no set point = system off
        return False

    # calculate temperatures
    rc_model_temp = calc_rc_model_temperatures_no_heating_cooling(bpr, tsd, t, config)

    # True, if T_int < ta_hs_set, False, if T_int >= ta_hs_set
    return rc_model_temp['T_int'] < ta_hs_set - temp_tolerance


def has_cooling_demand(bpr, tsd, t, config):
    """
    This function checks whether the building R-C-Model has a cooling demand according to the procedure in SIA 2044.
    R-C-Model temperatures are calculated with zero cooling power and checked versus the set-point temperature.
    Function includes a temperature tolerance according to the precision of the result reporting.

    :py:func: `cea.demand.rc_model_SIA.calc_rc_model_temperatures_no_heating_cooling`

    Author: Gabriel Happle
    Date: FEB 2017

    :param bpr: Building Properties
    :type bpr: BuildingPropertiesRow
    :param tsd: Time series data of building
    :type tsd: dict
    :param t: time step / hour of the year
    :type t: int
    :return: True or False
    :rtype: bool
    """

    temp_tolerance = 0.001  # temperature tolerance of temperature sensor (°C),
    #  i.e. heating is turned on if temperature is temp_tolerance below the set point
    # tolerance is consistent with maximum temperature difference that can be reported with the current precision
    # of the demand *.csv file

    ta_cs_set = tsd['ta_cs_set'][t]
    if np.isnan(ta_cs_set):
        # no set point = system off
        return False

    # calculate temperatures
    rc_model_temp = calc_rc_model_temperatures_no_heating_cooling(bpr, tsd, t, config)

    # True, if temperature w/o conditioning is higher than cooling set point temperature, else False
    return rc_model_temp['T_int'] > ta_cs_set + temp_tolerance


def has_sensible_cooling_demand(t_int_0, tsd, t):
    temp_tolerance = 0.001  # temperature tolerance of temperature sensor (°C),
    #  i.e. heating is turned on if temperature is temp_tolerance below the set point
    # tolerance is consistent with maximum temperature difference that can be reported with the current precision
    # of the demand *.csv file

    ta_cs_set = tsd['ta_cs_set'][t]
    if np.isnan(ta_cs_set):
        # no set point = system off
        return False

    # True, if temperature w/o conditioning is higher than cooling set point temperature, else False
    return t_int_0 > ta_cs_set + temp_tolerance


def has_sensible_heating_demand(t_int_0, tsd, t):
    temp_tolerance = 0.001  # temperature tolerance of temperature sensor (°C),
    #  i.e. heating is turned on if temperature is temp_tolerance below the set point
    # tolerance is consistent with maximum temperature difference that can be reported with the current precision
    # of the demand *.csv file

    ta_hs_set = tsd['ta_hs_set'][t]
    if np.isnan(ta_hs_set):
        # no set point = system off
        return False

    # True, if T_int < ta_hs_set, False, if T_int >= ta_hs_set
    return t_int_0 < ta_hs_set - temp_tolerance





# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3.8.1
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def lookup_f_hc_cv_heating(bpr):

    # 3.1.8.1 in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    # look up factor
    f_hc_cv = bpr.hvac['convection_hs']

    return f_hc_cv


def lookup_f_hc_cv_cooling(bpr):

    # 3.1.8.1 in SIA 2044 / Korrigenda C1 zum Merkblatt SIA 2044:2011 / Korrigenda C2 zum Mekblatt SIA 2044:2011

    # look up factor
    f_hc_cv = bpr.hvac['convection_cs']

    return f_hc_cv


# use the optimized (numba_cc) versions of the functions in this module if available
try:
    # import Numba AOT versions of the functions above, overwriting them
    from .rc_model_sia_cc import (calc_phi_m, calc_phi_c, calc_theta_c, calc_phi_m_tot, calc_phi_a, calc_theta_m,
                                 calc_h_ea, calc_theta_m_t, calc_theta_ea, calc_h_em, calc_h_3)
except ImportError:
    # fall back to using the python version
    # print('failed to import from rc_model_sia_cc.pyd, falling back to pure python functions')
    pass
