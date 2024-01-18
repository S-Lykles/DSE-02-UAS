from DSE import const
from DSE.structures.center_of_gravity import class_two_cg_estimation
from DSE.structures.center_of_gravity import cg_per_mission
import numpy as np
# from DSE.

def locations():
    "The current datum point is set at the nose of the craft in the x-direction. whilst in the Z-direction it is located"
    "at the botum of the main fuselage of the UAV"
    " l  is the distance between force and centre of gravity in the x-direction"
    " h  is the distance between force and centre of gravity in the z-direction"
    "All the X distance are measured from the datum point, the most forward point of the 6x6 ground surface"
    "All the Z distances are measured from the datum point, the center of the nose heigth"
    # Xcg = class_two_cg_estimation(True, False, False,  False)[1][0] #double check this function
    # Zcg = class_two_cg_estimation(True, False, False,  False)[1][2]#double check this function
    # print("Caution: The values about CG are approximate, need to be revisited. python file: DSE/structures/locations")

    cg1 = cg_per_mission(True, True, False, True, False)[1][0]
    cg2 = cg_per_mission(True, False, True, False, True)[1][0]

    if cg1>=cg2:
        Xcg, Ycg, Zcg = cg_per_mission(True, True, False, True, False)[1]
    else:
        Xcg, Ycg, Zcg = cg_per_mission(True, False, True, False, True)[1]
    # Xfr = 0.6 #PLACEHOLDER,  guessed value
    # Xaft = 3 #PLACEHOLDER,  guessed value
    # Xac =1.1 #PLACEHOLDER,  guessed value
    # Xac_v = 3.2  #PLACEHOLDER,  guessed value
    # Xh =4.3 #PLACEHOLDER,  guessed value
    # Zp = 0.7 #PLACEHOLDER,  guessed value
    # Zac = 1 #PLACEHOLDER,  guessed value
    # Zh = 1.6 #PLACEHOLDER,  guessed value
    Xfr = 0.5 #intial value
    Xaft = 3 #intial value
    Xac = 2.0828 #intial value
    Xac_v = 5  #intial value
    Xh = 5.1 #intial value
    Zp = -0.165 #intial value
    Zac = 0.17 #intial value
    Zh = 0.79 #intial value


    # print("Caution: The values about about distances are just guessed values, must be revisited. python file: DSE/structures/locations")

    X_lemac = 3 #PLACEHOLDER,  guessed value

    l_fr  = Xcg - Xfr   # Xfr is the distance of the front rotor
    l_aft = Xaft - Xcg  # Xaft is the distance of the aft rotor
    l_acw = Xcg - Xac   # Xacw is the distance of the aerodynamic center of the wing
    l_h   = Xh - Xcg    # Xh is the distance of the aerodynamic center of the horizontal tail
    h_p   = Zp - Zcg    # Zp is the position of the propellor
    h_acw = Zac - Zcg   # Zac is the position of the aerodynamic centre of the wing
    h_h   = Zh - Zcg    # Zh is the position of the horizontal tail
    z_h   = Zh - Zac
    # l_v   is found in vertical tail sizing function in tail_sizing file
    # Xac_v is the position of the aerodynamic center of the vertical tail.


    return l_fr, l_aft, l_acw, l_h, h_p, h_acw, h_h, z_h, X_lemac, Xcg, Zac, Zh
print('xcg', locations()[-3])