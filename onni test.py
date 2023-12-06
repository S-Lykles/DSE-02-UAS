import numpy as np

eff_GB  = 0.96 #Gearbox efficiency from literature
eff_GEN = 0.96 #Generator efficiency from literature
eff_EM  = 0.96 #Electro motor efficiency from literature
eff_PM  = 0.99 #Power management module efficiency from literature
eff_FC  = 0.76 #Hydrogen fuel cell efficiency from literature

def W_engine(P_cruise, P_max, eff_GB, eff_GEN, eff_EM, eff_PM):
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

    for P in P_req_powertrains:
        W_2stroke = 8.804e-7*P**4 - 1.577e-4*P**3 + 8.233e-3*P**2 + 0.504*P
        W_2stroke_arr = np.append(W_2stroke_arr, W_2stroke)
        W_4stroke = 8.733e-7*P**4 - 2.858e-4*P**3 + 2.363e-2*P**2 + 0.460*P
        W_4stroke_arr = np.append(W_4stroke_arr, W_4stroke)
        W_rotary  = 9.331e-4*P**2 + 0.625*P
        W_rotary_arr = np.append(W_rotary_arr, W_rotary)

    W_all_configs = np.array([W_2stroke_arr, W_4stroke_arr, W_rotary_arr])

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


def W_battery_electric(P_cruise, t_cruise, P_loiter, t_loiter, P_max, t_atPmax, E_rho_bat, E_rho_H, eff_EM, eff_PM, eff_FC):
    '''
    Calculation of the battery weight for fully electric propulsion.
    Parameters
    ----------
    P_cruise: float
        power required in kW for cruise computed from power curves
    t_cruise: float
        time needed in s to travel 185km and back at cruise speed
    P_loiter: float
        power required in kW for loitering. Computed from power curves
    t_loiter: float
        time needed to loiter for the endurance mission communication relay
    P_max: float
        maximum power in kW during operations computed from power curves
    t_atPmax: float
        time in s for how long the maximum power is needed to perform ONE take-off
    E_rho_bat: float
        energy density in MJ/kg
    eff_EM: float
        efficiency of electro motor from literature (DOI: 10.2514/6.2018-4228)
    eff_PM: float
        efficiency of power management module from literature (DOI: 10.2514/6.2018-4228)

    Outputs
    -------
    W_bat_supply: float
        Weight needed for batteries in kg for the supply delivery mission
    W_bat_endurance: float
        Weight needed for batteries in kg for the communication relay mission         
    '''
    #batteries
    eff = eff_PM*eff_EM
    E_needed_supply = (P_cruise*t_cruise + 4*P_max*t_atPmax)/eff #kJ
    E_needed_endurance = (P_cruise*t_cruise + 2*P_max*t_atPmax + P_loiter*t_loiter)*eff #kJ
    W_bat_supply = E_needed_supply/1000 / E_rho_bat #kg
    W_bat_endurance = E_needed_endurance/1000 / E_rho_bat #kg

    #hydrogen
    eff_H = eff_FC*eff_PM*eff_EM
    E_needed_supply_H = (P_cruise * t_cruise + 4 * P_max * t_atPmax) / eff_H  # kJ
    E_needed_endurance_H = (P_cruise * t_cruise + 2 * P_max * t_atPmax + P_loiter * t_loiter) * eff_H  # kJ
    W_H_supply = E_needed_supply_H/1000 / E_rho_H
    W_H_tank_supply = W_H_supply * 11.5
    W_H_endurance = E_needed_endurance_H/1000 / E_rho_H
    W_H_tank_endurance = W_H_supply



    return W_bat_supply, W_bat_endurance, W_H_supply, W_H_endurance









P_cruise = 10
t_cruise = 185000/40
P_loiter = 5
t_loiter = 10*3600
P_max = 30
t_atPmax = 250
E_rho_bat = 1.06 #MJ/kg
E_rho_H = 119.93 #MJ/kg LHV
a, b, c, d = W_battery_electric(P_cruise, t_cruise, P_loiter, t_loiter, P_max, t_atPmax, E_rho_bat, E_rho_H, eff_EM, eff_PM, eff_FC)
e = W_engine(P_cruise, P_max, eff_GB, eff_GEN, eff_EM, eff_PM)
f = W_battery_hybrid(P_cruise, P_max, t_atPmax, E_rho_bat, eff_GB, eff_EM, eff_PM)
print(a, b, c, d)
print(e)
print(f)