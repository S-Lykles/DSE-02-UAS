import numpy as np
import pandas as pd
from DSE.Rotor_Sizing import *
from DSE import const
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

file_dir = Path(__file__).parent

def data_file():
    data = pd.read_csv(file_dir / r"clarky.csv")

    x_axis = 'Alpha'
    y_axis = 'Cl'

    data.plot(x=x_axis, y=y_axis, marker='', linestyle='-')

    plt.title('Title')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.grid()
    plt.show()

    return

data_file()

def lift_coefficient(r):
    # Replace this with your actual lift coefficient function
    # For illustration, let's assume a simple linear function: 1.0 + 0.2 * r
    return 1.0 + 0.2 * r

# Define the chord length as a function of r
def chord_length(r):
    return 0.03 - (0.06 * r)

def lift_blade(r):
    Cl = 0.7370
    rho = 1.225
    omega = (5000*2*np.pi) / 60
    c = 0.03
    return 0.5 * Cl * rho * (omega * r)**2 * c

def drag_blade(r):
    Cd = 0
    rho = 1.225
    omega = (5000 * 2 * np.pi) / 60
    c = 0.03
    return 0.5 * Cd * rho * (omega * r) ** 2 * c

lower_limit = 0.0
upper_limit = 0.5

result, error = quad(lift_blade, lower_limit, upper_limit)

print(f"Result of integration: {result}")
print(f"Estimate of error: {error}")

r_values = np.linspace(lower_limit, upper_limit, 1000)
integrand_values = lift_blade(r_values)

plt.plot(r_values, integrand_values, label='Lift')
plt.fill_between(r_values, 0, integrand_values, alpha=0.2)

plt.title('Lift distribution over blade')
plt.xlabel('Radius [r]')
plt.ylabel('Lift [N]')
plt.grid()
plt.legend()
plt.show()