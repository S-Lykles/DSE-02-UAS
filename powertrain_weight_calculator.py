import numpy as np
import pandas as pd

#Constants needed from constants file
#-------------------------------------------------------------------------------------------
eff_GB  = 0.96 #Gearbox efficiency from literature
eff_GEN = 0.96 #Generator efficiency from literature
eff_EM  = 0.96 #Electro motor efficiency from literature
eff_PM  = 0.99 #Power management module efficiency from literature
eff_FC  = 0.76 #Hydrogen fuel cell efficiency from literature

spec_tank_W = 11.5 #unit weight of hydrogen fuel tank per unit weight of liquid hydrogen from literature
spec_P_fuelcell = 0.3#kW/kg specific power of a fuell cell, basically how much power can a fuell cell deliver per kg of fuel cell weight from literature
spec_P_fuelcell_fut = 8#kW/kg from literature
E_rho_bat = 1.44 #MJ/kg batteries of 1.75 MJ/kg were found (https://sionpower.com/files/Company-Brochure-21B.pdf) so I went a bit lower
E_rho_H = 119.93 #MJ/kg #Liquid hydrogen specific energy for Low Heating Value CHECK IF LHV OR HHV IS NEEDED

#Flight characteristics from power curves
#-------------------------------------------------------------------------------------------
P_cruise = 10
t_cruise = 185000/40
P_loiter = 5
t_loiter = 10*3600
P_max = 30
t_atPmax = 250


def W_engine(P_cruise, P_max, eff_GB, eff_GEN, eff_EM, eff_PM, spec_P_fuelcell, spec_P_fuelcell_fut):
    '''
    Calculation of the engine weight for different powertrain configurations (conventional, turbo-electric, parallel and serial) and engine types

    Parameters
    ----------
    P_cruise: float
        power required in kW for cruise computed from power curves
    P_max: float
        maximum power in kW during operations computed from power curves
    eff_GB: float
        efficiency of gearbox from literature (DOI: 10.2514/6.2018-4228)
    eff_GEN: float
        efficiency of generator from literature (DOI: 10.2514/6.2018-4228)
    eff_EM: float
        efficiency of electro motor from literature (DOI: 10.2514/6.2018-4228)
    eff_PM: float
        efficiency of power management module from literature (DOI: 10.2514/6.2018-4228)

    Outputs
    -------
    W_all_configs: numpy array
        an array which contains the weights for each engine type and each configuration
    '''


    eff_con = eff_GB
    eff_turbel = eff_GB*eff_GEN*eff_PM*eff_EM
    eff_ser = eff_GB*eff_GEN*eff_PM*eff_EM
    eff_par = eff_GB

    P_con = P_max/eff_con
    P_turbel = P_max/eff_turbel
    P_ser = P_cruise/eff_ser
    P_par = P_cruise/eff_par

    P_req_powertrains = np.array([P_con, P_turbel, P_ser, P_par])

    W_2stroke_arr = np.array([])
    W_4stroke_arr = np.array([])
    W_rotary_arr = np.array([])
    W_turbshft_arr = np.array([])
    W_fuelcell_arr = np.array([])
    W_fuelcell_fut_arr = np.array([])

    for P in P_req_powertrains:
        W_2stroke = 8.804e-7*P**4 - 1.577e-4*P**3 + 8.233e-3*P**2 + 0.504*P
        W_2stroke_arr = np.append(W_2stroke_arr, W_2stroke)
        W_4stroke = 8.733e-7*P**4 - 2.858e-4*P**3 + 2.363e-2*P**2 + 0.460*P
        W_4stroke_arr = np.append(W_4stroke_arr, W_4stroke)
        W_rotary  = 9.331e-4*P**2 + 0.625*P
        W_rotary_arr = np.append(W_rotary_arr, W_rotary)
        W_turbshft = -0.00001*P**2 + 0.2781*P + 5.0058
        W_turbshft_arr = np.append(W_turbshft_arr, W_turbshft)
        #Hydrogen computationt, added on Sieds' request
        W_fuelcell = P / spec_P_fuelcell
        W_fuelcell_arr = np.append(W_fuelcell_arr, W_fuelcell)
        W_fuelcell_fut = P / spec_P_fuelcell_fut
        W_fuelcell_fut_arr = np.append(W_fuelcell_fut_arr, W_fuelcell_fut)

    W_all_configs = np.array([W_2stroke_arr, W_4stroke_arr, W_rotary_arr, W_turbshft_arr, W_fuelcell_arr, W_fuelcell_fut_arr])

    return W_all_configs


def W_battery_hybrid(P_cruise, P_max, t_atPmax, E_rho_bat, eff_GB, eff_EM, eff_PM):
    '''
    Calculation of the battery weight for different hybrd powertrain configurations (parallel and serial). It is assumed that the battery generates the extra electricity needed for VTOL compared to cruise power (which is provided by the combustion engine)

    Parameters
    ----------
    P_cruise: float
        power required in kW for cruise computed from power curves
    P_max: float
        maximum power in kW during operations computed from power curves
    t_atPmax: float
        time in s for how long the maximum power is needed to perform ONE take-off
    E_rho_bat
        Energy density of battery in MJ/kg
    eff_GB: float
        efficiency of gearbox from literature (DOI: 10.2514/6.2018-4228)
    eff_EM: float
        efficiency of electro motor from literature (DOI: 10.2514/6.2018-4228)
    eff_PM: float
        efficiency of power management module from literature (DOI: 10.2514/6.2018-4228)

    Outputs
    -------
    W_bat: numpy array
        an array which contains the weights for the batteries, on the first index is the serial hybrid battery weight, on the second index the parallel hybrid battery weight
    '''
    eff_ser_bat = eff_PM*eff_EM
    eff_par_bat = eff_PM*eff_EM*eff_GB

    P_ser_bat = (P_max - P_cruise)/eff_ser_bat #kW
    P_par_bat = (P_max - P_cruise)/eff_par_bat #kW

    P_bat = np.array([P_ser_bat, P_par_bat]) #kW
    E_bat = P_bat/1000*t_atPmax #MJ
    W_bat = E_bat/E_rho_bat

    return W_bat


def W_electric(P_cruise, t_cruise, P_loiter, t_loiter, P_max, t_atPmax, E_rho_bat, E_rho_H, eff_EM, eff_PM, eff_FC, spec_tank_W, spec_P_fuelcell):
    """
    Calculates battery weights and liquid hydrogen system weight for different missions if only electric propulsion is used.

    :param P_cruise: Cruise power in kW
    :param t_cruise: Cruise time for travelling 185km two times in s
    :param P_loiter: Loiter power in kW
    :param t_loiter: Cruise time for loitering in s
    :param P_max: Maximum power needed (during VTOL) in kW
    :param t_atPmax: time that P_max is needed for ONE take-off
    :param E_rho_bat: battery energy density
    :param E_rho_H: liquid hydrogen energy density
    :param eff_EM: electro motor efficiency
    :param eff_PM: power management module efficiency
    :param eff_FC: fuel cell conversion efficiency
    :param spec_tank_W: fuel tank weight per unit weight of liquid hydrogen
    :return:
        W_bat_supply: battery weight for supply delivery mission in the case of fully electric propulsion using batteries
        W_bat_endurance: battery weight for communication relay mission in the case of fully electric propulsion using batteries
        W_H_supply: weight of liquid hydrogen propulsion system (weight of liquid hydrogen and weight of storage tank) for supply delivery mission in case of electric propulsion from hydrogen fuel cell
        W_H_endurance: weight of liquid hydrogen propulsion system (weight of liquid hydrogen and weight of storage tank) for communication relay mission in case of electric propulsion from hydrogen fuel cell
    """
    #batteries
    eff = eff_PM*eff_EM
    E_needed_supply = (P_cruise*t_cruise + 4*P_max*t_atPmax)/eff #kJ
    E_needed_endurance = (P_cruise*t_cruise + 2*P_max*t_atPmax + P_loiter*t_loiter)/eff #kJ
    W_bat_supply = E_needed_supply/1000 / E_rho_bat #kg
    W_bat_endurance = E_needed_endurance/1000 / E_rho_bat #kg

    #hydrogen
    eff_H = eff_PM*eff_EM #The fuel cell efficiency doesn't belong here, right? (Onni)
    P_peak = P_max/eff_H
    W_fuelcell = P_peak / spec_P_fuelcell
    W_fuelcell_fut = P_peak / spec_P_fuelcell_fut
    E_needed_supply_H = (P_cruise * t_cruise + 4 * P_max * t_atPmax) / eff_H / eff_FC  # kJ
    E_needed_endurance_H = (P_cruise * t_cruise + 2 * P_max * t_atPmax + P_loiter * t_loiter) * eff_H / eff_FC # kJ

    W_H_supply = E_needed_supply_H/1000 / E_rho_H #kg
    W_H_supply += (spec_tank_W*W_H_supply + W_fuelcell) #kg adds the fueltank weight and fuel cell weight to hydrogen fuel weight
    W_H_endurance = E_needed_endurance_H/1000 / E_rho_H #kg
    W_H_endurance += (spec_tank_W*W_H_endurance + W_fuelcell) #kg

    W_H_fut_supply = E_needed_supply_H/1000 / E_rho_H #kg
    W_H_fut_supply += (spec_tank_W*W_H_fut_supply + W_fuelcell_fut) #kg adds the fueltank weight and fuel cell weight to hydrogen fuel weight
    W_fut_H_endurance = E_needed_endurance_H/1000 / E_rho_H #kg
    W_fut_H_endurance += (spec_tank_W*W_fut_H_endurance + W_fuelcell_fut) #kg


    return W_bat_supply, W_bat_endurance, W_H_supply, W_H_endurance, W_H_fut_supply, W_fut_H_endurance



def table_hybrid_propulsion_weights():
    """
    Produces a pandas dataframe (a fancy table) which presents all the weights for conventional and hybrid propulsion (so no full electric propulsion!) and engine types.
    :return:
    """
    #Import all the weights for the conventional, turbo-electric, serial hybrid and parallel hybrid powertrains with different engine types
    W_engines_hybrid = W_engine(P_cruise, P_max, eff_GB, eff_GEN, eff_EM, eff_PM, spec_P_fuelcell, spec_P_fuelcell_fut)
    W_batteries_hybrid = W_battery_hybrid(P_cruise, P_max, t_atPmax, E_rho_bat, eff_GB, eff_EM, eff_PM)

    #create array of zeros to add battery weight to engine weight of serial and parallel hybrid
    zero_array = np.zeros_like(W_engines_hybrid)
    zero_array[:, 2:] = W_batteries_hybrid
    result_array = W_engines_hybrid + zero_array
    result_array = np.round(result_array, 1)
    result_dataframe = pd.DataFrame(data = result_array, columns = ['Conventional', 'Turbo-electric', 'Serial Hybrid', 'Parallel Hybrid'], index = ['2-Stroke', '4-Stroke', 'Rotary', 'Turboshaft', 'Liquid Hydrogen Current Tech', 'Liquid Hydrogen Future Tech'])

    return result_dataframe

def table_electric_propulsion_weights():
    """
    Produces a pandas dataframe (a fancy table) which presents all the weights for conventional and hybrid propulsion (so no full electric propulsion!) and engine types.
    :return:
    """
    W_bat_supply,  W_bat_endurance, W_H_supply, W_H_endurance, W_fut_H_supply, W_fut_H_endurance = W_electric(P_cruise, t_cruise, P_loiter, t_loiter, P_max, t_atPmax, E_rho_bat, E_rho_H, eff_EM, eff_PM, eff_FC, spec_tank_W, spec_P_fuelcell)
    result_array = np.array([[W_bat_supply, W_bat_endurance],[W_H_supply, W_H_endurance], [W_fut_H_supply, W_fut_H_endurance]])
    result_array = np.round(result_array, 1)
    result_dataframe = pd.DataFrame(data = result_array, columns = ['Supply Delivery', 'Communication Relay'], index = ['Battery Powered', 'Liquid Hydrogen Current Tech', 'Liquid Hydrogen Future Tech'])
    return result_dataframe


print(table_hybrid_propulsion_weights())
print(table_electric_propulsion_weights())

#LATEX BOOKTABS CODE
#print(table_hybrid_propulsion_weights().to_latex(index=False, formatters={"name": str.upper}, float_format="{:.1f}".format,))