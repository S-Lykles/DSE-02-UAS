import numpy as np
import matplotlib.pyplot as plt
from Inputs_Preq_rotorcraft import *


v = np.arange(5, 50, 0.01)
P_profile_drag_lst = []
P_induced_1_lst = []
P_parasite_lst = []
P_loss_lst = []
P_tot_req_level_lst = []

for i in v:
    #Profile drag
    advanced_ratio = i / (omega*R)

    C_t = W / (rho * np.pi * R**2 * (omega*R)**2)
    Cl_bar = 6.6*(C_t / sig_max)
    alpha_m =Cl_bar / Cl_alpha

    CD_p1 = 0.0087 - 0.0216*alpha_m + (0.4 * alpha_m **2)
    CD_p2 = 0.011 + 0.4*(alpha_m**2)
    CD_p3 = 0.009 + 0.73*(alpha_m**2)
    C_D_p = (CD_p1 + CD_p3 + CD_p2) / 3

    P_hov = (1/8)*sig_max*C_D_p*rho*(omega*R)**3*np.pi*(R**2)

    P_profile_drag = P_hov * (1 + 4.65*(advanced_ratio**2))
    P_profile_drag_lst.append(P_profile_drag)

    #induced drag power
    v_ih = np.sqrt(W / (2*rho*np.pi*R**2))
    v_ibar = 1 / i
    v_i = v_ibar * v_ih
    P_induced_1 = k * T_level * v_i
    P_induced_1_lst.append(P_induced_1)

    #parsite power
    P_parasite = 0.5 * rho * A_eq * (i**3)
    P_parasite_lst.append(P_parasite)

    #power loss
    P_loss = 0.04*(P_profile_drag + P_induced_1 + P_parasite)
    P_loss_lst.append(P_loss)

    #total power
    P_tot_req_level = P_profile_drag + P_induced_1 + P_parasite + P_loss
    P_tot_req_level_lst.append(P_tot_req_level)

plt.figure(dpi=600)
plt.plot(v, P_profile_drag_lst, label = 'profile drag')
plt.plot(v, P_induced_1_lst, label = 'induced')
plt.plot(v, P_parasite_lst, label = 'parasite')
plt.plot(v, P_loss_lst, label = 'loss')
plt.plot(v, P_tot_req_level_lst, label = 'total')
#plt.xlim(20)
plt.grid()
plt.legend()
plt.show()
#Graph
print(np.shape(P_parasite_lst))
print(P_parasite_lst)
