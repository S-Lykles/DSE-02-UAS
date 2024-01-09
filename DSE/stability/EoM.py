import numpy as np
from DSE import const
from DSE.structures.center_of_gravity import class_two_cg_estimation
from DSE.Locations import locations

def EoM():
    """Assumes tail above cg"""
    "Moment equations around the center of gravity"
    "T_p = the propellor thrust"
    "T_fr = the rotor thrust of the front two rotors"
    "T_aft = the rotor thrust of the aft two rotors"
    h_p = locations()[4]
    l_fr = locations()[0]
    l_aft = locations()[1]
    h_ac = locations()[5]
    l_ac = locations()[2]
    h_h = locations()[6]
    l_h = locations()[3]
    #Force equations
    F_x = T_p + (L_w + L_h + L_fus)*np.sin(alpha* deg2rad) - (D_w + D_h + D_fus)*np.cos(alpha* deg2rad) - W*np.sin(theta)
    F_z = -T_fr_left - T_fr_right - T_aft_left - T_aft_right - (L_w + L_h + L_fus)*np.cos(alpha* deg2rad) + (D_w + D_h + D_fus)*np.sin(alpha* deg2rad) + W*np.cos(theta)

    "Moment around axes of the C.G."
    #Cm_y / c_mean  = Cm_w + Cm_fus + Cm_h *(c_mean_h / c_mean) + (CL_w * h_ac - CD_w * l_ac + (CL_h * h_h + CD_h * l_h) * (S_h / S) *(c_mean_h / c_mean)* (v_h / v) ** 2) * np.sin(alpha * deg2rad) + (             CL_w * l_ac + CD_w * h_ac + (CD_h * h_h - CL_h * l_h) * (S_h / S) *(c_mean_h / c_mean)* (v_h / v) ** 2) * np.cos alpha * deg2rad) - T_p * h_p / (0.5*rho*S*v**2) + (2 * T_fr_left * l_fr - 2 * T_aft_left * l_aft) / (0.5* rho *S *v**2)
    M_x = T_fr_left*w_rot_left + T_aft_left* w_rot_left - T_fr_right * w_rot_right - T_aft_right * w_rot_right
    M_y = M_ac_w + M_ac_h - T_p*h_p + T_fr_right*l_fr + T_fr_left*l_fr - T_aft_right*l_aft -T_aft_left*l_aft +(L_w*h_ac +D_w*l_ac - L_h*h_h - D_h*l_h)*np.sin(alpha* deg2rad) +(L_w*l_ac +D_w*h_ac -L_h*l_h +D_h*h_h)*np.cos(alpha* deg2rad)

    return F_x, F_z, M_x, M_y



def rotation_matrix(theta, psy, phi, deg2rad):
    "A very nice matrix that converts the body axis angles into the earth axis -(0_0)- "

    r_theta = theta * deg2rad
    r_phi = phi * deg2rad
    r_psy = psy * deg2rad
 # R{Xd,Ye,Ze,Phi,Theta,Psy}^t
    R = np.array([[np.cos(r_theta)*np.cos(r_psy), np.cos(r_psy)*np.sin(r_theta)*np.sin(r_phi)-np.sin(r_psy)*np.cos(r_phi), np.cos(r_psy)*np.sin(r_theta)*np.cos(r_phi) + np.sin(r_psy)*np.sin(r_phi), 0,0,0],
                  [np.cos(r_theta)*np.sin(r_psy), np.sin(r_psy)*np.sin(r_theta)*np.sin(r_phi) + np.cos(r_psy)*np.cos(r_phi), np.sin(r_psy)*np.sin(r_theta)*np.cos(r_phi) - np.cos(r_psy)*np.sin(r_phi),0,0,0],
                  [-1*np.sin(r_theta), np.sin(r_phi)*np.cos(r_theta), np.cos(r_phi)*np.cos(r_theta),0,0,0],
                  [0,0,0,1,np.sin(r_phi)*np.sin(r_theta)/np.cos(r_theta) , np.sin(r_theta)*np.cos(r_phi)/np.cos(r_theta)],
                  [0,0,0,0,np.cos(r_phi) , -1*np.sin(r_phi)],
                  [0,0,0,0,np.sin(r_phi)/np.cos(r_theta)]])
    return R



def inertial_matrix(mass):
    """Yet another matrix! We blessed!
    Define inputs for class_two_cg_estimation in time"""
    I_xx, I_yy, I_zz = class_two_cg_estimation()[3]

    M = np.array([[mass, 0.0, 0.0, 0.0, 0.0, 0.0],
                  [0.0, mass, 0.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, mass, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, I_xx, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, I_yy, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0, I_zz]])
    return M

def gyroscopic_matrix(mass, velocity_vector):
    """Gyroscopic matrix, Define inputs for class_two_cg_estimation, velocity vector called p in lit, this might be a wee bit wrong due to symmetry assumptions"""

    I_xx, I_yy, I_zz = class_two_cg_estimation()[3]
    u, v, w, p, q, r = velocity_vector

    C = np.array([[0.0, -mass * r, mass * q, 0.0, 0.0, 0.0],
                  [mass * r, 0.0, -mass * p, 0.0, 0.0, 0.0],
                  [-mass * q, mass * p, 0.0, 0.0, 0.0, 0.0],
                  [0.0, -mass * w, mass * v, 0.0, I_yy * r, -I_zz * q],
                  [mass * w, 0.0, -mass * u, I_zz * r, 0.0, -I_xx * p],
                  [-mass * v, mass * u, 0.0, I_xx * q, -I_yy * p, 0.0]])
    return C
