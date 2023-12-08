#Parameters listed for the weight estimation
#General parameters:
MTOW = 160 *2.20462262 #lbs
payload_com = 20 *2.20462262 #lbs
payload_sup = 50 *2.20462262 #lbs
W_fuel = 15 *2.20462262 #lbs
W_L = MTOW - W_fuel - payload_sup #Landing mass


# Dimensional parameters
l = 3.35 #m
d = 0.8 #Diameter fuselage m
l_sm = 1  # Shock strut length for main gear [ft]
l_sn = 1  # Shock strut length for nose gear [ft]


# Wing parameters
b = 6 *3.28084 #ft
S = 3.763 *3.28084**2 #Wing surface main wing in ft^2
A = b**2/S #Aspect ratio main wing


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


# Propulsion system parameters
N_electric = 4 #number of rotors electrically driven in particular configuration
N_gas = 1 #number of rotors gasoline driven
W_electro_motor = 0.1836*(P_hov_max/N_electric)+ 2.7076 #Dependent on power required
W_gas_motor = 10#Dependent on power required
W_rotor_gas = 3#Dependent on rotor design
W_rotor_electric = 1#Dependent on rotor design
W_fuel_sys = W_fuel/9 #Literature research
W_battery = 10 #kg, Dependent on propulsion configuration (battery weight for electric VTOL)
W_generator = 15#Dependent on configuration


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

#rpm =  #expected rpm of the rotors
#R = #radius main rotor in ft
# C = #average chord of a blade in ft
# N = #number of blades on rotor
# V_tip = #Blade tip speed for cruise conditions in ft/s
# S_ch = #Total blade area of one rotor in ft^2
# A_ch = #Disk are of one rotor in ft