import numpy as np
import matplotlib.pyplot as plt
from const import rho0, MTOW
from aero.main_aero import b, S, cl, cd_comp


AR = b**2/S
step = 0.1
sigma = 30*0
sigma_rad = np.deg2rad(sigma)
e = 1.78*(1-0.045*AR**0.68)-0.64

CL, CD = cl[140:], cd_comp[140:]
cd0 = np.min(cd_comp)
print(CL)
V = np.sqrt(MTOW*2/(rho0*S*CL))
dd=1/2*rho0*V**2*S
D_par = dd*(cd0)
D_ind = dd*CL**2/(np.pi*AR*e)
D = D_par + D_ind

T_par = D_par/np.cos(sigma_rad)
T_ind = D_ind/np.cos(sigma_rad)
T = T_par + T_ind
# T = D/np.cos(sigma_rad)

L = MTOW-T*np.sin(sigma_rad)

P_req_par = T_par*V
P_req_ind= T_ind*V
P_req = P_req_par+P_req_ind


plt.plot(V,P_req_par, label = 'Parasite power')
plt.plot(V,P_req_ind, label = 'Induced power')
plt.plot(V,P_req, label = 'Power Required')
plt.xlabel("Velocity [ m/s]")
plt.ylabel("Power [W]")
plt.title("Power curve for transition mode")
plt.legend()
plt.show()



# sigma_l = np.arange(0,90.1,0.1)
# sigma_rad_l = np.deg2rad(sigma_l)
# L_max = []
# L_min = []
# P=[]
# V_p =[]
# for i in range(len(sigma_rad_l)):
#     T_l = D/np.cos(sigma_rad_l[i])
#     l_max = MTOW-np.min(T_l)*np.sin(sigma_rad_l[i])
#     l_min = MTOW-np.max(T_l)*np.sin(sigma_rad_l[i])
#     L_max.append(l_max)
#     L_min.append(l_min)
#     p = (T_l)#*V
#     P.append(np.max(p))
#     V_p.append(V[p==np.max(p)])

# print(P)
#
# sigma_rad_p = np.radians(30)
# Lift = MTOW-D/np.cos(sigma_rad_p)*np.sin(sigma_rad_p)
# # plt.plot(V, Lift, label='Lift MAX')
# plt.plot(sigma_l, P, label='Power at tilt angle')
# # plt.plot(V_p, P, label='Power at tilt angle')
# # plt.xlabel('Velocity')
# plt.xlabel('Tilt angle')
# plt.ylabel('Power')
# plt.legend()
# plt.show()