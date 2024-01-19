import numpy as np
from DSE import const
from DSE.Locations import locations
from DSE.aero import aero_constants
from DSE.stability.tail_sizing import horizontal_tail_sizing, elevator_surface_sizing
from control.matlab import *
from matplotlib import pyplot as plt

print('File control_derivatives.py needs to be revisted and use the inputs from other files, instead of using the actual values.')
PRInt = False



V = 42 # placeholder
Vh = V
d_dt = 99999 # placeholder time step
T = 288.15 - 0.0065 * 500

##   !!! Imported Values !!!
##  !! still undifend and not placed !!
Theta_0 = 2*np.pi/180 # placeholder
CL_alpha_cruise = 9999 # placeholder, input from aerodynamics CL_alpaha at CL cruise.
CL_alpha_CL_0 = 9999 # placeholder, input from aerodynamics CL_alpha at CL=0
Cd0_h = 9999 # placeholder, input from aerodynamics

# Base
rho = const.rho0
m = const.total_mass
M =0.12

# Propulsion
C_t = 9999 # placeholder, input from propulsion
CT_alpha = 0 # placeholder, input from propulsion
T1 = -9999 # placeholder, input from propulsion
T2 = -9999 # placeholder, input from propulsion
T3 = -9999 # placeholder, input from propulsion
T4 = -9999 # placeholder, input from propulsion
Tp = -9999 # placeholder, input from propulsion
Lw = -9999 # placeholder, input from propulsion
CY_alpha_p = -9999
J = -9999
sigma = -9999
ks = 1.14
ka = 0.4
# T_C = C_t/J**2
# a = np.sqrt(1+8*T_C/np.pi-1)/2
# fa = (1+a)*((1+a)+(1+2*a)**2)/(1+(1+2*a)**2
mo = 0.95*2*np.pi
# I = 3/4*mo*integral_blade_c
# CY_alpha_p = ks*fa*sigma*I1/(1+ka*sigma*I1)
lp = -9999
Vp = -9999


# Wing properties
S = aero_constants.S
Sh = 1.26
de_da = 0.42313505699610365
b = aero_constants.b
c_bar = aero_constants.c_bar
Cd = aero_constants.CD_cruise[0] # placeholder, input from aerodyamics
CL_w = aero_constants.CL_cruise
CL0 = aero_constants.CL_0
CD0 = aero_constants.CD_0
sweep_ang_25_c = aero_constants.sweep_ang_25_c_rad
CL_alpha_w = aero_constants.CL_alpha_wing
Cd_alpha = aero_constants.CD_alpha_wing
Cd0_w = aero_constants.CD0_wing
Cr_w = aero_constants.c_root
taper_w = aero_constants.taper
e = aero_constants.e

# Tail properties
    # Horizontal
# Sh = horizontal_tail_sizing()[0] # placeholder horizontal tail surface
AR_h = 6.8 # import from horizontal
Sh = 0.538 # import from horizontal
bh = 2.3 # import from horizontal
CL_h = aero_constants.Cl_cruise_h
Cl_alpha_h = aero_constants.Cl_alpha_h
eta = 0.95
de_da = horizontal_tail_sizing()[4]
V_h = Vh
beta = np.sqrt(1-M**2)
sweep_ang_50_c_rad = aero_constants.sweep_ang_50_c_rad
CL_alpha_h = Cl_alpha_h*AR_h/(2+np.sqrt(4+(AR_h*beta/eta)**2*(1+np.tan(sweep_ang_50_c_rad)**2/beta**2)))
    # Vertical
bv =  0.60188057586457
Sv = 0.1811301138015332
lv = 2.452651830421447
AR_v = 1.9
sweep_v = 11
eta_v = 0.90    # assumption
Cl_alpha_v = aero_constants.Cl_alpha_v
Vv = V
lv = -9999
#CL_alpha_v1 = aero_constants.CL_alpha_v
#CL_alpha_v2 = CL_alpha_v1

# Initial Calucaltions
M0 = V/(np.sqrt(1.4*287.15*T))
CDM = Cd * M0 / (1-M0**2)
beta = np.sqrt(1-M**2)
eta = 0.95
CT = 0.01
CT_alpha = 0

CD_alpha_w = aero_constants.CL_alpha_wing * 2 * CL_w / (np.pi * b*b/S*aero_constants.e)
Ixx = -9999 # placeholder, input from structures
Iyy = -9999 # placeholder, input from structures
Ixz = -9999 # placeholder, input from structures
# Kxz = -9999
Jxy = Ixz/(m*b**2)
Ky_2 = Iyy/(m*c_bar**2)
mu_c = m/(rho*S*c_bar)
mu_b = m/(rho*S*b)
Kx_2 = Ixx/(m*b**2)
Dc = c_bar/V * d_dt
Db = b/V * d_dt
Lh = 0.8*Lw

T1 = -9999 # placeholder, input from propulsion
T2 = -9999 # placeholder, input from propulsion
T3 = -9999 # placeholder, input from propulsion
T4 = -9999 # placeholder, input from propulsion
Tp = -9999 # placeholder, input from propulsion
Lw = -9999 # placeholder, input from propulsion
Lh = -9999 # placeholder, input from propulsion
q_rad= -9999 # placeholder, input for control

l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h,z_h,X_lemac, Xcg, Zac, Zh = locations()
Z_m = -9999



zv = -9999

vtol=False
if vtol:
    CX0 = 0
    CZ0 = -1*(T1+T2+T3+T4)/(0.5*rho*S*V**2)
    CXu = 0
    CZu = 0
    CMu = 0
    CXalpha = 0
    CZalpha = 0
    Cmalpha = 0
    CXalphadott = 0
    CZalphadott = 0
    Cmalphadott = 0
    CZq = -9999
    Cmq = -9999

    CYr = -9999
    Cnr = -9999
    Cmq = -9999
    CYr = 0
    CYp = -1.87  # Ref(lit) : K.W. Booth. Effect of horizontal-tail chord on the calculated subsonic span loads and stability derivatives of isolated unswept tail assemblies in sideslip and steady roll. Technical report, NASA Memo 4-1-59 L, 1959.
    Clr = -9999
    CXu = 0
    CZu = 0
    Cmu = 0
    CXq = 0

else:
    CX0 = Tp / (0.5*rho*S*V**2) - Cd
    CZ0 = -CL_w - CL_h*(Sh/S) * (V_h/V**2)
    CXalpha = CL_alpha_w * aero_constants.alpha_0 + aero_constants.CL_initial_conditions - Cd_alpha + CT_alpha
    CZalpha = - aero_constants.CL_alpha_wing - aero_constants.Cl_alpha_h * Sh / S
    Cmalpha = aero_constants.CL_alpha_wing * l_acw / aero_constants.c_bar - aero_constants.CL_alpha_h * l_h * Sh / S / aero_constants.c_bar - CD_alpha_w * Zac / aero_constants.c_bar  + CT_alpha * Zh / aero_constants.c_bar
    CXalphadott = 0
    CZalphadott = - Cl_alpha_h * (V_h/V)**2 * de_da * aero_constants.S_h * l_h / S / aero_constants.c_bar
    Cmalphadott = - Cl_alpha_h * (V_h/V)**2 * de_da * aero_constants.S_h * l_h**2 / S / aero_constants.c_bar/ aero_constants.c_bar

    CZq = -CL_alpha_w - CL_alpha_h*l_h*Sh/(c_bar*S)*(Vh/V)**2
    Cmq = CL_alpha_w * l_acw**2/c_bar**2 +-CL_alpha_h*l_h*Sh/(S*c_bar)*(Vh/V)**2
    # Cmq = CL_alpha_w * l_acw**2/c_bar**2 - CL_alpha_h
    if PRInt == True:
        print('CZq = ', CZq,'Cmq = ',  Cmq)
        print('CXalpha = ', CXalpha, 'CZalpha = ', CZalpha, 'Cmalpha = ', Cmalpha)
        print('CXalphadott = ', CXalphadott, 'CZalphadott = ', CZalphadott, 'Cmalphadott = ', Cmalphadott)
    #CYr_v1 = 2*(Vv/V)**2*Sv*lv/(S*b)*CL_alpha_v1
    #CYr_v2 = 2*(Vv/V)**2*Sv*lv/(S*b)*CL_alpha_v2
    #CYr = CYr_v1+CYr_v2 + 2*CY_alpha_p*(Vp/V)**2*Sv*lp/(S*b)*0
    #Clr = CL_w+CL_h*Sh*bh/(S*b)*(Vh/V)**2 - zv/b*(CYr_v1+CYr_v2)
    #Cnr = lv/b*CYr_v1+lv/b*CYr_v2
    if PRInt == True:
        print('CYr:',CYr_v1, CYr_v2)
        print('Clr:',CL_w,CL_h,Sh,bh,S,b,Vh,V,zv,b,CYr_v1,CYr_v2)
        print('Cnr:',lv,b,CYr_v1,CYr_v2)
        print('CYr_v1:',Vv,V,Sv,lv,S,b,CL_alpha_v2)

    CYp = -2*  8/(np.pi*3) *eta_v**2 * (bv*Sv/(b*S))* (Cl_alpha_v * AR_v) / (2 + np.sqrt(4 + (((AR_v * beta) / eta) ** 2) * (((np.tan(sweep_v * const.deg2rad)) ** 2 / beta ** 2) + 1)))
    Clp = -1* (((CL_alpha_w + Cd0_w)*Cr_w*b)/(24*S) * (1+3*taper_w)) - (( (( (Cl_alpha_h*AR_h)/(2+np.sqrt(4+(AR_h*beta/eta)**2))) + Cd0_h))/6)  # Radians
    Cnp = -lv / b * CYp - 1 / 8 * (CL_w + CL_h * Sh / S * bh / b)
    CXu = -3 * CD0 - 3 * CL0 * np.tan(Theta_0) - M0 * CDM # Caughey, D. A., Introduction to Aircraft Stability and Control Course Notes for AE5070, 2011
    CZu = -M0**2 / (1 - M0**2)  * (CL_w + CL_h * (Sh/S))
    # CMu = (2/c_bar) * (CL_w * l_acw - CL_h * l_h - Cd0_w * Zac + C_t * Z_m) * ((2 * Z_m)/(V * c_bar))
    CMu = M0**2 / (1 - M0**2) * (CL_w * (l_acw/c_bar) - CL_h * ((l_h * S)/(c_bar * Sh)) - Cd0_w * (Zac/c_bar))
    CXq = 0
    CXdelt_e = 0
    CZdelt_e = -CL_alpha_h*0.95*Sh*l_h/S*c_bar* elevator_surface_sizing()[0] * (aero_constants.c_bar/l_h)
    CMdelt_e = -CL_alpha_h*0.95*Sh*l_h/S*c_bar* elevator_surface_sizing()[0]
    CXdelt_t = 0
    CYdelt_t = 0
    CZdelt_t = 0
    CMdelt_t = 0

if vtol:
    b = b

else:
    P_symm =np.matrix( [[-2 * mu_c * c_bar / V, 0, 0, 0],
         [0, (CZalphadott - 2 * mu_c) * c_bar / V, 0, 0],
         [0, 0, -c_bar / V, 0],
         [0, Cmalphadott * c_bar / V, 0, - 2 * mu_c * Ky_2 * c_bar / V]])

    Q_symm = np.matrix([[-CXu, -CXalpha, -CZ0, 0],
         [-CZu, -CZalpha, CX0, -1*(CZq + 2*mu_c)],
         [0, 0, 0, -1],
         [-CMu, -Cmalpha, 0, -Cmq]])

    R_symm = np.matrix([[-CXdelt_e,CXdelt_t],
         [-CZdelt_e,CYdelt_t],
         [0,0],
         [-CMdelt_e, CMdelt_t]])

    #P_inv = np.linalg.inv(P_symm)
    #A = np.matmul(P_inv,Q_symm)
    A_symm = np.linalg.inv(P_symm) @ Q_symm
    B_symm = np.linalg.inv(P_symm) @ R_symm
    eig_val_symm = np.linalg.eig(A_symm)[0]


    #  Y = u_dott, w_dott, theta_dott, thata, delta_ele, delta_trim,

    C_symm = [[1 , 0 , 0 , 0],
              [0 , 1 , 0 , 0],
              [0 , 0 , 1 , 0],
              [0 , 0 , 0 , V/aero_constants.c_bar]]

    D_symm = [[0 , 0],
              [0 , 0],
              [0 , 0],
              [1 , 0]]
    print(A_symm,B_symm,C_symm,D_symm)
damping = True
if damping:
    sys = ss(A_symm,B_symm,C_symm,D_symm)
    state_matrix = np.eye(4)
    input_matrix = np.eye(2)
    K, S, E = lqr(sys, state_matrix, input_matrix)
    A_cl = A_symm - B_symm @ K
    sys_cl = ss(A_cl,B_symm,C_symm,D_symm)
    damping_system = sys_cl.feedback(K)
    K[0,0] = 0
    print("K",K)
    # The time vector
    tend = 100
    dt = 0.1
    t = np.arange(0,tend+dt, dt)
    #==============================================================================
    # Step input or initial condition?
    step_input = False
    initial_condition = True
    if step_input:
        # The input vector.
        ttussen = 1  # ttussen is the time that 1 should be present.
        u01 = np.zeros(len(t))
        for i in range(int(ttussen/dt)):
            u01[i] = 10 *np.pi/180
        #print('check number of ones in x01',sum(x01))
        u02 = np.zeros(len(t))
        s = (len(t),2)
        u0 = np.zeros(s)
        for i in range(len(t)):
            u0[i] = u01[i],u02[i]
        y, time, x = lsim(damping_system, u0, t)
    if initial_condition:
        x0 = np.array([[5],          # initial codnitions for u, alpha, theta and q respectively
                       [np.pi/4],
                       [0],
                       [np.pi/2* aero_constants.c_bar / V]])

        y, time = initial(sys_cl, t, x0)
    plt.plot(t,y[:,0],label = "u")
    plt.plot(t,y[:,1],label = "alpha")
    plt.plot(t,y[:,2],label = "theta")
    plt.plot(t,y[:,3],label = "q")
    plt.title('Initial')
    plt.xlabel('t')
    plt.ylabel('y')
    plt.legend()
    plt.show()

    print(pole(sys_cl))

    poles = True
    if poles:
        poles = pole(sys_cl)

        # Compute transmission zeros
        zeros = zero(sys_cl)

        # Plot the pole-zero map
        plt.scatter(np.real(poles), np.imag(poles), marker='x', label='Poles')
        plt.scatter(np.real(zeros), np.imag(zeros), marker='o', label='Zeros')

        plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
        plt.axvline(0, color='black', linewidth=0.5, linestyle='--')

        plt.title('Pole-Zero Map')
        plt.xlabel('Real')
        plt.ylabel('Imaginary')
        plt.legend()
        plt.grid(True)
        plt.show()

# print("""t
#           CXu, CZu, CMu,
#           CXalpha, CZalpha, 0,
#           CXq, CZq, Cmq,
#           CXalphadott, CZalphadott, 0
# """, B_symm)



## Short period
# 1. V is constant which makes the udˆ equal to zero.
# 2. As it this motion applies to the steady horizontal and symmetric flight it can be assumed that γ , Cx0
# is zero and thus pitch rate is zero
# 3. All inputs are 0, thus the right hand side of the equations of motion equal 0

## Phugoid motion
# 1. Constant velocity, thus, α˙ = 0
# 2. Steady symmetric flight so Cx0 = 0
# 3. Since the period is very long it can be assumed all vertical accelerations are zero and thus the change
# in pitch angle q˙ = 0
# 4. Since 2µc >> CZq in the phugoid motion we can neglect CZq

## Aperiodic roll
# 1. There is only rolling motion present. β = 0 and r = 0

## Dutch roll
# 1. 4µb >> CYr thus CYr can be neglected
# 2. Since the roling motion is not very significant for simplicity it can be neglected. ϕ = 0 and p = 0
# 3. The change in sideslip angle does not have a big effect on the normal force, therefore, CYβ˙ = Cnβ˙ = 0
# 4. The aircraft’s centre of gravity moves in a straight line. This follows from the assumption that there is
# only rotation in yaw

## Aperiodic spiral
# 1. Again due to the slow nature of the oscilation it can be assumed that all linear and angular accelerations are zero. Thus Db = 0
# 2. The only relevant variables are yaw, roll and pitch. Therefore, CYr and CYp
# can be neglected.





