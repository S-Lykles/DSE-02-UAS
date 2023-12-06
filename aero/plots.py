import matplotlib.pyplot as plt
import numpy as np

# Read Xfoil data from the file
data_4412 = np.loadtxt(r"C:\Users\Florian Tjepkema\Documents\Aerospace Engineering\AeroSpace 2023-2024\DSE\DSE-02-UAS\aero\naca4412.txt")
data_23012 = np.loadtxt(r"C:\Users\Florian Tjepkema\Documents\Aerospace Engineering\AeroSpace 2023-2024\DSE\DSE-02-UAS\aero\naca23012.txt")
data_632A015 = np.loadtxt(r"C:\Users\Florian Tjepkema\Documents\Aerospace Engineering\AeroSpace 2023-2024\DSE\DSE-02-UAS\aero\naca63(2)A015.txt")
data_63215 = np.loadtxt(r"C:\Users\Florian Tjepkema\Documents\Aerospace Engineering\AeroSpace 2023-2024\DSE\DSE-02-UAS\aero\naca63215.txt")
data_63512 = np.loadtxt(r"C:\Users\Florian Tjepkema\Documents\Aerospace Engineering\AeroSpace 2023-2024\DSE\DSE-02-UAS\aero\naca63512.txt")

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

# Plot the data
plt.figure(figsize=(10, 6))
plt.subplot(1,3,1)
plt.plot(alpha_4412, cl_4412, label='NACA 4412')
plt.plot(alpha_23012, cl_23012, label='NACA 23012')
plt.plot(alpha_632A015, cl_632A015, label='NACA 632A015')
plt.plot(alpha_63215, cl_63215, label='NACA 63215')
plt.plot(alpha_63512, cl_63512, label='NACA 63512')
plt.title('Xfoil Data - Lift Coefficient vs. Angle of Attack')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Lift Coefficient (Cl)')
plt.grid(True)
plt.legend()

plt.subplot(1,3,2)
plt.plot(cd_4412, cl_4412, label='NACA 4412')
plt.plot(cd_23012, cl_23012, label='NACA 23012')
plt.plot(cd_632A015, cl_632A015, label='NACA 632A015')
plt.plot(cd_63215, cl_63215, label='NACA 63215')
plt.plot(cd_63512, cl_63512, label='NACA 63512')
plt.title('Xfoil Data - Lift Coefficient vs. Drag Coefficient')
plt.xlabel('Drag Coefficient (Cd)')
plt.ylabel('Lift Coefficient (Cl)')
plt.grid(True)
plt.legend()

plt.subplot(1,3,3)
plt.plot(alpha_4412, cm_4412, label='NACA 4412')
plt.plot(alpha_23012, cm_23012, label='NACA 23012')
plt.plot(alpha_632A015, cm_632A015, label='NACA 632A015')
plt.plot(alpha_63215, cm_63215, label='NACA 63215')
plt.plot(alpha_63512, cm_63512, label='NACA 63512')
plt.title('Xfoil Data - Moment Coefficient vs. Angle of Attack')
plt.xlabel('Angle of Attack (degrees)')
plt.ylabel('Moment Coefficient (Cm)')
plt.grid(True)
plt.legend()

plt.show()
