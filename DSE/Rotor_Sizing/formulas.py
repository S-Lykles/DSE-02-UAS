import numpy as np
import pandas as pd
import csv
from DSE.Rotor_Sizing import *
from DSE import *

csv_file_path = 'DSE/Rotor_Sizing/clarky.csv'
data = pd.read_csv(csv_file_path)

print(data)

def lift_blade(cl, cd,  rho, omega, r, c):
    L = (cl * rho * (omega*r)**2 * c) / 2
    D = (cd * rho * (omega*r)**2 * c) / 2
    RPM = (2 * np.pi) / omega
    return