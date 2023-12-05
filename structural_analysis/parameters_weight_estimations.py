#Parameters listed for the weight estimation
#General parameters:

MTOW = 160 #kg
payload_com = 20 #kg
payload_sup = 50 #kg
W_fuel = 20 #kg
W_L = MTOW - W_fuel - payload_sup #Landing mass


# Dual phase based on the FD180P(most similar specs in database with comparable main conf)
# Cessna method used for this particular main configuration

# Dimensional parameters
l = 3.35 #m
d = 0.8 #Diameter fuselage m
l_sm =  #Shock strut length fro main gear
l_sn = #Shock strut length for nose gear

# Wing parameters
b = 6 #m
S = 3.763 #Wing surface main wing in squared meters
A = b**2/S #Aspect ratio main wing

#Empennage parameters
# S_h = #Surface area horizontal tailwing
# A_h = #Aspect ratio horizontal tailwing
# S_v = #Surface area vertical tailwing
# A_v = #Aspect ratio vertical tailwing
# t_rh = #Maximum root thickness of the horizontal tailwing
# t_rv = #Maximum root thickness of the vertical tailwing
# chord_sweep_angle = #Sweep angle of the quarter chord vertical wing


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
#W_gas_motor = #Dependent on power required
#W_rotor = #Dependent on rotor design
W_fuel_sys = W_fuel/9 #Literature research
W_battery = 24 #kg, Dependent on propulsion configuration (battery weight for electric VTOL)
#W_generator = #Dependent on configuration


# Add up all the weights of the propulsion to find an estimation for the weight of this subsystem
W_prop = N_electric*(W_electro_motor + W_rotor) + N_gas*(W_gas_motor + W_rotor) + W_fuel_sys + W_battery + W_generator

# Mission and avionics system parameters
W_missioncomputer = 0.027*MTOW
W_nav_sys = 0.013*MTOW
W_flt_ctrl = 0.024*MTOW

#Add all masses to find avionics subsystem weight
W_avionics = W_missioncomputer + W_nav_sys + W_flt_ctrl
# W_com =
# W_sensors =
# W_payload_sys =
# W_hydraulics =
# W_electric_circuit =


# Compound helicopter parameters\

#rpm =  #expected rpm of the rotors
#R = #radius main rotor
# C = #average chord of a blade
# N = #number of blades on rotor
# V_tip = #Blade tip speed for cruise conditions
# S_ch = #Total blade area of one rotor
# A_ch = #Disk are of one rotor