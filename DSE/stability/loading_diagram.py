from DSE.aero import aero_constants
from DSE.Locations import locations
import matplotlib.pyplot as plt
from DSE.structures.center_of_gravity import cg_per_mission
import numpy as np

x_lemac = 0.5
def bar_cal(x, x_lemac = 0.5):

    # x_lemac = locations()[8]
    mac = aero_constants.c_bar
    x_bar = (x-x_lemac)/mac

    return x_bar

def centre_of_gravity(OEM, fuel_pd, fuel_le, payload_pd, payload_le):
    mass = OEM+fuel_pd+fuel_le+payload_pd+payload_le
    cg =  OEM*1.8+fuel_pd*1.5+fuel_le*1.5+payload_pd*1.8+payload_le*1.8
    return cg, mass

def loading_diagram_extremes():
    l_f = 6
    xlemac_lf = []
    x_lemac= np.arange(1,4,0.01)

    x_cg_bar_min =[]
    x_cg_bar_max =[]
    for i in range(len(x_lemac)):
        xlemac = x_lemac[i]
        xlemac_lf.append(xlemac/l_f)

        mass_plus_fuel_pd, cg_fuel_pd, _ = cg_per_mission(True, True, False, False, False)
        x_cg_fuel_pd_bar = bar_cal(cg_fuel_pd[0], xlemac)

        mass_plus_payload_pd, cg_payload_pd, _  = cg_per_mission(True, True, False, True, False)
        x_cg_payload_pd_bar = bar_cal(cg_payload_pd[0], xlemac)

        mass_plus_fuel_le, cg_fuel_le, _ = cg_per_mission(True, False, True, False, False)
        x_cg_fuel_le_bar = bar_cal(cg_fuel_le[0], xlemac)

        mass_plus_payload_le, cg_payload_le, _  = cg_per_mission(True, False, True, False, True)
        x_cg_payload_le_bar = bar_cal(cg_payload_le[0], xlemac)

        x_cg_bar = [x_cg_fuel_le_bar, x_cg_fuel_pd_bar, x_cg_payload_pd_bar, x_cg_payload_le_bar]
        x_cg_bar_min.append(np.min(x_cg_bar))
        x_cg_bar_max.append(np.max(x_cg_bar))


    slope = (x_cg_bar_min[0]-x_cg_bar_min[-1])/(xlemac_lf[0]-xlemac_lf[-1])
    return x_cg_bar_min, x_cg_bar_max, xlemac_lf, slope

loading_diagram_extremes()