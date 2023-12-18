import matplotlib.pyplot as plt
import numpy as np
from plot_setting import tex_fonts, set_size

# Read Xfoil data from the file
# Dual Phase
data_4412 = np.loadtxt(r"aero/naca4412.txt")
data_23012 = np.loadtxt(r"aero/naca23012.txt")
data_632A015 = np.loadtxt(r"aero/naca63(2)A015.txt")
data_63215 = np.loadtxt(r"aero/naca63215.txt")
data_63512 = np.loadtxt(r"aero/naca63512.txt")
# Tilt rotor
data_2410 = np.loadtxt(r"aero/naca2410.txt")
data_2412 = np.loadtxt(r"aero/naca2412.txt")
data_16 = np.loadtxt(r"aero/ag16.txt")
data_26 = np.loadtxt(r"aero/ag26.txt")
data_1223 = np.loadtxt(r"aero/s1223.txt")

# Extract columns
alpha_4412 = data_4412[:, 0]      # Angle of attack
cl_4412 = data_4412[:, 1]         # Lift coefficient
cd_4412 = data_4412[:, 2]         # Drag coefficient
cdp_4412 = data_4412[:, 3]        # Pressure drag coefficient
cm_4412 = data_4412[:, 4]         # Moment coefficient

alpha_23012 = data_23012[:, 0]    # Angle of attack
cl_23012 = data_23012[:, 1]       # Lift coefficient
cd_23012 = data_23012[:, 2]       # Drag coefficient
cdp_23012 = data_23012[:, 3]      # Pressure drag coefficient
cm_23012 = data_23012[:, 4]       # Moment coefficient

alpha_632A015 = data_632A015[:, 0]    # Angle of attack
cl_632A015 = data_632A015[:, 1]       # Lift coefficient
cd_632A015 = data_632A015[:, 2]       # Drag coefficient
cdp_632A015 = data_632A015[:, 3]      # Pressure drag coefficient
cm_632A015 = data_632A015[:, 4]       # Moment coefficient

alpha_63215 = data_63215[:, 0]    # Angle of attack
cl_63215 = data_63215[:, 1]       # Lift coefficient
cd_63215 = data_63215[:, 2]       # Drag coefficient
cdp_63215 = data_63215[:, 3]      # Pressure drag coefficient
cm_63215 = data_63215[:, 4]       # Moment coefficient

alpha_63512 = data_63512[:, 0]    # Angle of attack
cl_63512 = data_63512[:, 1]       # Lift coefficient
cd_63512 = data_63512[:, 2]       # Drag coefficient
cdp_63512 = data_63512[:, 3]      # Pressure drag coefficient
cm_63512 = data_63512[:, 4]       # Moment coefficient

alpha_2410 = data_2410[:, 0]      # Angle of attack
cl_2410 = data_2410[:, 1]         # Lift coefficient
cd_2410 = data_2410[:, 2]         # Drag coefficient
cdp_2410 = data_2410[:, 3]        # Pressure drag coefficient
cm_2410 = data_2410[:, 4]         # Moment coefficient

alpha_2412 = data_2412[:, 0]      # Angle of attack
cl_2412 = data_2412[:, 1]         # Lift coefficient
cd_2412 = data_2412[:, 2]         # Drag coefficient
cdp_2412 = data_2412[:, 3]        # Pressure drag coefficient
cm_2412 = data_2412[:, 4]         # Moment coefficient

alpha_16 = data_16[:, 0]      # Angle of attack
cl_16 = data_16[:, 1]         # Lift coefficient
cd_16 = data_16[:, 2]         # Drag coefficient
cdp_16 = data_16[:, 3]        # Pressure drag coefficient
cm_16 = data_16[:, 4]         # Moment coefficient

alpha_26 = data_26[:, 0]      # Angle of attack
cl_26 = data_26[:, 1]         # Lift coefficient
cd_26 = data_26[:, 2]         # Drag coefficient
cdp_26 = data_26[:, 3]        # Pressure drag coefficient
cm_26 = data_26[:, 4]         # Moment coefficient

alpha_1223 = data_1223[:, 0]      # Angle of attack
cl_1223 = data_1223[:, 1]         # Lift coefficient
cd_1223 = data_1223[:, 2]         # Drag coefficient
cdp_1223 = data_1223[:, 3]        # Pressure drag coefficient
cm_1223 = data_1223[:, 4]         # Moment coefficient

# Plot the data Dual Phase
plt.rcParams.update(tex_fonts)
# size = set_size(fraction=1, subplots=(1,1))
# plt.figure(figsize=(size[0]*0.7,size[1]))
plt.figure(figsize=(10, 6))
plt.subplot(1,3,1)
plt.plot(alpha_4412, cl_4412, label='NACA 4412')
plt.plot(alpha_23012, cl_23012, label='NACA 23012')
plt.plot(alpha_632A015, cl_632A015, label='NACA 632A015')
plt.plot(alpha_63215, cl_63215, label='NACA 63215')
plt.plot(alpha_63512, cl_63512, label='NACA 63512')
plt.title('Xfoil Data - $C_l$ vs alpha')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift Coefficient ($C_l$)')
plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.tight_layout()
plt.legend()

plt.subplot(1,3,2)
plt.plot(cd_4412, cl_4412, label='NACA 4412')
plt.plot(cd_23012, cl_23012, label='NACA 23012')
plt.plot(cd_632A015, cl_632A015, label='NACA 632A015')
plt.plot(cd_63215, cl_63215, label='NACA 63215')
plt.plot(cd_63512, cl_63512, label='NACA 63512')
plt.title('Xfoil Data - $C_l$ vs $C_d$')
plt.xlabel('Drag Coefficient ($C_d$)')
plt.ylabel('Lift Coefficient ($C_l$)')
plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.tight_layout()
plt.legend()

plt.subplot(1,3,3)
plt.plot(alpha_4412, cm_4412, label='NACA 4412')
plt.plot(alpha_23012, cm_23012, label='NACA 23012')
plt.plot(alpha_632A015, cm_632A015, label='NACA 632A015')
plt.plot(alpha_63215, cm_63215, label='NACA 63215')
plt.plot(alpha_63512, cm_63512, label='NACA 63512')
plt.title('Xfoil Data - $C_m$ vs alpha')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Moment Coefficient ($C_m$)')
plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.tight_layout()
plt.legend()

plt.savefig('aero/dual_airfoil.pdf')
plt.show()

# Plot the data Tilt Wing
plt.figure(figsize=(10, 6))
plt.subplot(1,3,1)
plt.plot(alpha_4412, cl_4412, label='NACA 4412')
plt.plot(alpha_2410, cl_2410, label='NACA 2410')
plt.plot(alpha_2412, cl_2412, label='NACA 2412')
plt.plot(alpha_16, cl_16, label='AG 16')
plt.plot(alpha_26, cl_26, label='AG 26')
plt.title('Xfoil Data - $C_l$ vs alpha')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift Coefficient ($C_l$)')
plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.minorticks_on()
plt.tight_layout()
plt.legend()

plt.subplot(1,3,3)
plt.plot(alpha_4412, cm_4412, label='NACA 4412')
plt.plot(alpha_2410, cm_2410, label='NACA 2410')
plt.plot(alpha_2412, cm_2412, label='NACA 2412')
plt.plot(alpha_16, cm_16, label='AG 16')
plt.plot(alpha_26, cm_26, label='AG 26')
plt.title('Xfoil Data - $C_m$ vs alpha')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Moment Coefficient ($C_m$)')
plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.tight_layout()
plt.legend()

plt.subplot(1,3,2)
plt.plot(cd_4412, cl_4412, label='NACA 4412')
plt.plot(cd_2410, cl_2410, label='NACA 2410')
plt.plot(cd_2412, cl_2412, label='NACA 2412')
plt.plot(cd_26, cl_26, label='AG 26')
plt.plot(cd_16, cl_16, label='AG 16')
plt.title('Xfoil Data - $C_l$ vs $C_d$')
plt.xlabel('Drag Coefficient ($C_d$)')
plt.ylabel('Lift Coefficient ($C_l$)')
plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.minorticks_on()
plt.legend()


plt.tight_layout()
plt.savefig('aero/tilt_airfoil.pdf')
plt.show()
