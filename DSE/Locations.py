from DSE import const
from DSE.structures.center_of_gravity import class_two_cg_estimation

def locations():
    "The current datum point is set at the nose of the craft in the x-direction. whilst in the Z-direction it is located"
    "at the botum of the main fuselage of the UAV"
    " l  is the distance between force and centre of gravity in the x-direction"
    " h  is the distance between force and centre of gravity in the z-direction"
    "All the X distance are measured from the datum point, the most forward point of the 6x6 ground surface"
    "All the Z distances are measured from the datum point, the center of the nose heigth"
    Xcg = class_two_cg_estimation()[1][0]
    Zcg = class_two_cg_estimation()[1][2]


    l_fr  = Xcg - Xfr   # Xfr is the distance of the front rotor
    l_aft = Xaft - Xcg  # Xaft is the distance of the aft rotor
    l_acw = Xcg - Xac   # Xacw is the distance of the aerodynamic center of the wing
    l_h   = Xh - Xcg    # Xh is the distance of the aerodynamic center of the horizontal tail
    h_p   = Zp - Zcg    # Zp is the position of the propellor
    h_acw = Zac - Zcg   # Zac is the position of the aerodynamic centre of the wing
    h_h   = Zh - Zcg    # Zh is the position of the horizontal tail
    z_h   = Zh - Zac
    l_v   = Xac_v - Xcg # Xac_v is the position of the aerodynamic center of the vertical tail.

    return l_fr, l_aft, l_acw,l_h,h_p,h_acw,h_h,z_h,l_v