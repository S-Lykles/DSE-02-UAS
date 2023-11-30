from Inputs_Preq_rotorcraft import *
from Bucket_tool import generate_Preq
from rotor_sizing import rotor_sizing_tool
import matplotlib.pyplot as plt


DL = 200
N = 1
A_eq = 0.09

v = np.arange(5, 60, 0.01)

R, D_v, omega, T_level, sig_max = rotor_sizing_tool(DL, N)
Preq = generate_Preq(A_eq, R, D_v, omega, T_level, sig_max, v)

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
