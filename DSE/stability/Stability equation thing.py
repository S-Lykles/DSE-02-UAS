import numpy as np
from DSE import const

# this doc, will contain static and dynamic stability and controllabilty equations, with bases in Roskamp Airplane Design
# Just a temp doc all functions if used will be moved to the corresponding document


def elevator_sizing(c_bar=0.619,Cm_0=-0.111,Cm_alpha=-0.029,alpha=0,alpha_0=0,CL_alpha_h= ,Sh_S,l_h,vtrans=34,uh=31):
    # speed range ( Stall <-> Max + safety margin)

    #possible import ?  Cm_0 = Cm_ac - CL_alpha_h*(alhpa_0 - i_h)* x_h/c* (S_h/s)* (u_h/u)**2
    #possilbe import ?  Cm_alpha = d_Cm / d_CL * CL_alpha

    delta = 25  # Elevator deflection range ( -25 <-> 25 degrees)
    u = vtrans
    Cm_delta_el = -1*(Cm_0 + Cm_alpha*(alpha - alpha_0)) / (delta)
    Tau_el = -1*Cm_delta_el / CL_alpha_h * (bh/be) * 1/Sh_S * c_bar/l_h * (u/uh)**2

    return Tau_el



