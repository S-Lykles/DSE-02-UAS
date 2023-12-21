import numpy as np
from DSE import const

def distance_stability():
    "The current datum point is set at the nose of the craft in the x-direction. whilst in the Z-direction it is located"
    "at the botum of the main fuselage of the UAV"
    " l  is the distance between force and centre of gravity in the x-direction"
    " h  is the distance between force and centre of gravity in the z-direction"
    "All the X distance are measured from the datum point, the most forward point of the 6x6 ground surface"
    "All the Z distances are measured from the datum point, the center of the nose heigth"

    l_fr  = Xcg - Xfr   #'Xfr is the distance of the front rotor'
    l_aft = Xaft - Xcg  #'Xaft is the distance of the aft rotor'
    l_acw = Xcg - Xac  #'Xacw is the distance of the aerodynamic center of the wing'
    l_h   = Xh - Xcg    #'Xh is the distance of the aerodynamic center of the horizontal tail'
    h_p   = Zp - Zcg    #'Zp is the position of the propellor'
    h_acw = Zac - Zcg
    h_h   = Zh - Zcg

    return l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h


def EoM():
    "Moment equations around the center of gravity"
    "T_p = the propellor thrust"
    "T_fr = the rotor thrust of the front two rotors"
    "T_aft = the rotor thrust of the aft two rotors"
    #Force equations
    F_x = T_p + (L_w + L_h + L_fus)*np.sin(alpha* deg2rad) - (D_w + D_h + D_fus)*np.cos(alpha* deg2rad) - W*np.sin(theta)
    F_z = T_fr_left + T_fr_right + T_aft_left + T_aft_right + (L_w + L_h + L_fus)*np.cos(alpha* deg2rad) - (D_w + D_h + D_fus)*np.sin(alpha* deg2rad) - W*np.cos(theta)

    "Moment around y-axis of C.G."
    #Cm_y / c_mean  = Cm_w + Cm_fus + Cm_h *(c_mean_h / c_mean) + (CL_w * h_ac - CD_w * l_ac + (CL_h * h_h + CD_h * l_h) * (S_h / S) *(c_mean_h / c_mean)* (v_h / v) ** 2) * np.sin(alpha * deg2rad) + (             CL_w * l_ac + CD_w * h_ac + (CD_h * h_h - CL_h * l_h) * (S_h / S) *(c_mean_h / c_mean)* (v_h / v) ** 2) * np.cos alpha * deg2rad) - T_p * h_p / (0.5*rho*S*v**2) + (2 * T_fr_left * l_fr - 2 * T_aft_left * l_aft) / (0.5* rho *S *v**2)
    M_x = T_fr_left*w_rot_left + T_aft_left* w_rot_left - T_fr_right * w_rot_right - T_aft_right * w_rot_right
    M_y = M_ac_w + M_fus + M_ac_h - T_p*h_p + T_fr_right*l_fr + T_fr_left*l_fr - T_aft_right*l*aft -T_aft_left*l_aft +(L_w*h_ac +D_w*h_ac - L_h*h_h - D_h*l_h)*np.sin(alpha* deg2rad) +(L_w*l_ac +D_w*l_ac -L_h*l_h -D_h*h_h)*np.cos(alpha* deg2rad)

    return F_x, F_z, M_x, M_y



def rotation_matrix():
    r_theta = theta * deg2rad
    r_phi = phi * deg2rad
    r_psy = psy * deg2rad

    R = np.matrix{[np.cos(r_theta)*np.cos(r_psy), np.cos(r_psy)*np.sin(r_theta)*np.sin(r_phi)-np.sin(r_psy)*np.cos(r_phi), np.cos(r_psy)*np.sin(r_theta)*np.cos(r_phi) + np.sin(r_psy)*np.sin(r_phi), 0,0,0],
                  [np.cos(r_theta)*np.sin(r_psy), np.sin(r_psy)*np.sin(r_theta)*np.sin(r_phi) + np.cos(r_psy)*np.cos(r_phi), np.sin(r_psy)*np.sin(r_theta)*np.cos(r_phi) - np.cos(r_psy)*np.sin(r_phi),0,0,0],
                  [-1*np.sin(r_theta), np.sin(r_phi)*np.cos(r_theta), np.cos(r_phi)*np.cos(r_theta),0,0,0],
                  [0,0,0,1,np.sin(r_phi)*np.sin(r_theta)/np.cos(r_theta) , np.sin(r_theta)*np.cos(r_phi)/np.cos(r_theta)],
                  [0,0,0,0,np.cos(r_phi) , -1*np.sin(r_phi)],
                  [0,0,0,0,np.sin(r_phi)/np.cos(r_theta),np.cos(r_phi)/np.cos(r_theta)]}
