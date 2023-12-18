import numpy as np
import matplotlib.pyplot as plt
# from main import DL, R, W, N_r
# from plot_power import P_req_rotor
#! CODE IS BROKEN, NOT SUPPOSED TO WORK RIGHT NOW
print("WARNING: CODE IS BROKEN, NOT SUPPOSED TO WORK RIGHT NOW")
assert False

V_tip = 140*(R*2)**0.171
RPM = V_tip/R*60/(2*np.pi)

R_400 = np.sqrt((W/N_r)/(400*np.pi))
V_tip400 = 140*(R_400*2)**0.171
RPM_400 = V_tip400/R_400*60/(2*np.pi)
print('RPM 400', RPM_400)
N_r = 4
power = np.max(P_req_rotor) #PLACEHOLDER
# print(power)
torque = power / N_r / (RPM/60*2*np.pi)

print('Max. Torque = ', np.max(torque) )
plt.plot(DL,        torque)
plt.ylabel("Torque [Nm]")
plt.xlabel("DL [N/m^2]")
plt.show()

m_mot_ref = 3.9
torque_mot_ref = 60
m_mot = m_mot_ref*(np.max(torque)/torque_mot_ref)**(3/3.5)
print('Mass of motors = ', m_mot)

# torque = pd.DataFrame(torque)
# # Create a heatmap
# plt.figure(figsize=(10, 8))
# sns.heatmap(torque, annot=True, xticklabels=DL, yticklabels=R, cmap='viridis')
# plt.xlabel('Radius [m]')
# plt.ylabel('Disk loading [N/m^2]')
# plt.title('Torque [Nm]')
# plt.show()