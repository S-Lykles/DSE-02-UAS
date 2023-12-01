#Parameters listed for the weight estimation
#General parameters:

MTOW = 160 #kg
payload_com = 20 #kg
payload_sup = 50 #kg

# m_fuel =


# Dual phase based on the FD180P(most similar specs in database with comparable main conf)
# Cessna method used for this particular main configuration

# Dimensional parameters
l = 3.35 #m
wing_span = 6.5 #m
c_mean = 0.7 #m
S = l * c_mean #Wing surface in meters
A = l**2/S



# Specs FD180P
endurance = 10 #hrs
range = 500 #km
V_cruise = 110 #km/h
V_max = 130 #km/h
h_max = 5000 #m
n_ult =


# Propulsion system parameters
N_electric = 4 #number of rotors electrically driven in particular configuration
N_gas = 1 #number of rotors gasoline driven
W_electro_motor =
W_gas_motor =
W_rotor =
W_fuel_sys =
W_battery = 24 #kg
W_generator =

# Add up all the weights of the propulsion to find an estimation for the weight of this subsystem
W_prop = N_electric*(W_electro_motor + W_rotor) + N_gas*(W_gas_motor + W_rotor) + W_fuel_sys + W_battery + W_generator

# Mission and avionics system parameters
W_fuel = 20
W_com =
W_sensors =
W_payload_sys =
W_hydraulics =
W_boardcomputer =
W_electric_circuit =


