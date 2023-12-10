import numpy as np

# Unit factors
kw_to_hp = 1.341022
lb_to_kg = 2.20462262
m_to_ft = 0.3048

#Parameters listed for the weight estimation
#General parameters:
MTOW = 160 *2.20462262 #lbs
payload_com = 20 *2.20462262 #lbs
payload_sup = 50 *2.20462262 #lbs
W_fuel = 20 *2.20462262 #lbs
W_L = MTOW - W_fuel - payload_sup #Landing mass


# Dimensional parameters
l = 2.5 * 0.3048 # Length fuselage in ft
d = 0.8  * 0.3048 # Max diameter fuselage in ft
perimeter = d*np.pi # Max perimeter fuselage in ft
l_sm = 1  # Shock strut length for main gear [ft]
l_sn = 1  # Shock strut length for nose gear [ft]


# Wing parameters
b = 6 *3.28084 #ft
S = 3.763 *3.28084**2 #Wing surface main wing in ft^2
AR = b**2/S #Aspect ratio main wing


#Empennage parameters
# S_h = #Surface area horizontal tailwing in ft^2
# A_h = #Aspect ratio horizontal tailwing
# S_v = #Surface area vertical tailwing in ft^2
# A_v = #Aspect ratio vertical tailwing
# t_rh = #Maximum root thickness of the horizontal tailwing in ft
# t_rv = #Maximum root thickness of the vertical tailwing in ft
# chord_sweep_angle = #Sweep angle of the quarter chord vertical wing in radians


# Specs
endurance = 10 #hrs
range = 500 #km
V_cruise = 110 #km/h
V_max = 130 #km/h
h_max = 5000 #m
n_ult = 1.5
P_hov_max = 35 #kW Max power required during hover/take-off phase
P_cruise_max = 15 #kW Max power required during cruise


# # Propulsion system parameters dual phase
# N_electric = 4 #number of rotors electrically driven in particular configuration
# N_gas = 1 #number of rotors gasoline driven
# W_electro_motor = 0.1836*(P_hov_max/N_electric)+ 2.7076 #Dependent on power required
# W_gas_motor = 10#Dependent on power required
# W_rotor_gas = 3#Dependent on rotor design
# W_rotor_electric = 1#Dependent on rotor design
# W_fuel_sys = W_fuel/9 #Literature research
# W_battery = 0 #kg, Dependent on propulsion configuration (battery weight for electric VTOL)
# W_generator = 12#Dependent on configuration


# Mission and avionics system parameters
W_missioncomputer = 0.027*MTOW
W_nav_sys = 0.013*MTOW
W_flt_ctrl = 0.024*MTOW


# W_com =
# W_sensors =
# W_payload_sys =
# W_hydraulics =
# W_electric_circuit =


# Compound helicopter parameters
R = 1.4738579917347323 / m_to_ft#radius main rotor in [ft]
t_c = 0.12#average tickness over chord of a blade
N_blades = 4 #number of blades on rotor
V_tip = 140*(2*R*m_to_ft)**0.171 / m_to_ft #Blade tip speed for cruise conditions in [ft/s]
DL = MTOW/np.pi/R**2 # Disk loading
sigma = 0.06659266968898687
b_ch = 5 * m_to_ft# Wing span of compound helicopter [ft]
S_ch = 3.763 * m_to_ft**2 # Wing area of compound helicopter [ft^2]
AR_ch = b_ch**2 / S_ch # Aspect ratio wing compound helicopter [-]

# Propulsion system compound helicopter
N_electric = 2 #number of rotors electrically driven in particular configuration
N_gas = 1 #number of rotors gasoline driven
W_rotor_electric = 1 #Dependent on rotor design [kg]
W_fuel_sys = W_fuel/9 #Literature research [kg]
W_generator = 15#Dependent on configuration

