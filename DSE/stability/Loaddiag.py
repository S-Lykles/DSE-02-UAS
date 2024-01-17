import numpy as np
from DSE import const
from DSE.structures.center_of_gravity import class_two_cg_estimation
import matplotlib.pyplot as plt
# from DSE import plot_setting



def load_diagram_plot(empty, fuel, payload_supply, plots = False):
    """UPDATE IF MULTIPLE FUEL TANKS OR PAYLOAD LOCATIONS ARE DEFINED

    Create potato plot based on Cg calculation

    needs OEW
    needs fuel weight (fuel location estimation based on class II estimation)
    needs payload weight (weight location based on class II estimation, assumed to be concentrated in c_g of payload)"""

    weight = []
    c_g = []

    #c_g variation for fuel first, no payload
    for i in range(101):
        interim = (class_two_cg_estimation(empty, i/100 * fuel)[0:1, 0])
        weight.append(interim[0])
        c_g.append(interim[1])

    #c_g variation for fuel second, payload included\
    for i in range(101):
        interim = (class_two_cg_estimation(empty, i/100 * fuel, payload_supply)[0:1, 0])
        weight.append(interim[0])
        c_g.append(interim[1])

    #c_g variation for payload, no fuel
    for i in range(101):
        interim = (class_two_cg_estimation(empty, 0.0, i/100 * payload_supply)[0:1, 0])
        weight.append(interim[0])
        c_g.append(interim[1])

    #c_g variation for payload, fuel included
    for i in range(101):
        interim = (class_two_cg_estimation(empty, fuel, i/100 * payload_supply)[0:1, 0])
        weight.append(interim[0])
        c_g.append(interim[1])

    if plots == True:
        plt.plot(c_g, weight)
        plt.show()

    return weight, c_g
