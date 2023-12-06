# Function to calculate engine power needed based on efficiency
def calculate_required_engine_power(max_output_power, efficiencies):
    overall_efficiency = 1.0
    for efficiency in efficiencies.values():
        overall_efficiency *= efficiency
    return max_output_power / overall_efficiency

# Function to calculate engine weight based on the provided formula
def calculate_engine_weight(max_engine_power, engine_formula):
    return engine_formula(max_engine_power)

# Example formulas for engine weight calculations
def weight_2stroke(max_engine_power):
    return max_engine_power * 0.25

def engine_formula_2(max_engine_power):
    return max_engine_power * 0.3


eff_GB = 0.96 #Gearbox efficiency from literature
eff_GEN = 0.96 #Generator efficiency from literature
eff_EM = 0.96 #Electro motor efficiency from literature
eff_PM = 0.99 #Power management module efficiency from literature

# Define efficiencies for different configurations
efficiencies_conventional = {
    "gearbox": eff_GB
}

efficiencies_turbo_electric = {
    "gearbox": eff_GB
}

efficiencies_serial = {
    "gearbox": eff_GB,
    "generator": eff_GEN,
    "power management module": eff_PM,
    "electro motor": eff_EM,
}

efficiencies_parallel = {
    "gearbox": eff_GB,
}



# Calculate required engine power for Configuration 1
max_output_power = 1000  # Input max output power
required_power_conventional = calculate_required_engine_power(max_output_power, efficiencies_conventional)

# Choose an engine and calculate its weight based on the calculated engine power
chosen_engine = weight_2stroke  # Choose the desired engine formula
engine_weight = calculate_engine_weight(required_power_conventional, chosen_engine)

# Output the results
print(f"Required engine power for Configuration 1: {required_power_conventional}")
print(f"Engine weight for Configuration 1: {engine_weight}")
