import numpy as np
import matplotlib.pyplot as plt
from main import DL, R, W, N_r
from plot_power import P_req_rotor
from const import rho0


V_tip = 140*(R*2)**0.171
RPM = V_tip/R*60/(2*np.pi)

R_400 = np.sqrt((W/N_r)/(400*np.pi))
V_tip400 = 140*(R_400*2)**0.171
RPM_400 = V_tip400/R_400*60/(2*np.pi)
print('RPM 400', RPM_400)

power = np.max(P_req_rotor) #PLACEHOLDER
# print(power)
torque = power / N_r / (RPM/60*2*np.pi)

print('Max. Torque = ', np.max(torque) )
plt.plot(DL,        torque)
plt.ylabel("Torque [Nm]")
plt.xlabel("DL [N/m^2]")
plt.show()

m_mot_ref = 97
torque_mot_ref = 60
m_mot = m_mot_ref*(np.max(torque)/torque_mot_ref)**(3/3.5)
print('Mass of motors = ', m_mot, (6-np.min(R)*2*N_r-0.8)/N_r)

r_mot = np.sqrt(np.max(torque)*2/(m_mot*np.min(RPM)))
density_mot = 2710
d_mot = 0.25

if d_mot<r_mot*2:
    d_mot = r_mot*2
v_mot = m_mot/density_mot
l_mot = v_mot/(np.pi*(d_mot/2)**2)
print("motor:", np.max(torque), np.min(RPM), power)
print("Dimensions:", d_mot, l_mot)


# torque = pd.DataFrame(torque)
# # Create a heatmap
# plt.figure(figsize=(10, 8))
# sns.heatmap(torque, annot=True, xticklabels=DL, yticklabels=R, cmap='viridis')
# plt.xlabel('Radius [m]')
# plt.ylabel('Disk loading [N/m^2]')
# plt.title('Torque [Nm]')
# plt.show()


# Propeller efficiency in the forward flight:
V=42
T=power/V
A_disk = np.pi*(d_mot/2)**2
u0=V

eta_p = (power/N_r)/(RPM*torque)
# eta_p = 2/(1+(T/(A_disk*u0**2*rho0/2)+1))
print('efficiency',eta_p,(power/N_r),(RPM*torque))
