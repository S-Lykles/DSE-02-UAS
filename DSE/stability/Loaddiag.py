import numpy as np
from DSE import const
# Cannot see Cg location function from structural_analysis import

def datum_location():
    "The current datum point is set at the nose of the craft in the x-direction. whilst in the Z-direction it is located"
    "at the botum of the main fuselage of the UAV"
    " l  is the distance between force and centre of gravity in the x-direction"
    " h  is the distance between force and centre of gravity in the z-direction"
    "All the X distance are measured from the datum point, the most forward point of the 6x6 ground surface"
    "All the Z distances are measured from the datum point, the center of the nose heigth"

    l_fr  = Xcg - Xfr   'Xfr is the distance of the front rotor'
    l_aft = Xaft - Xcg  'Xaft is the distance of the aft rotor'
    l_acw = Xcg - Xacw  'Xacw is the distance of the aerodynamic center of the wing'
    l_h   = Xh - Xcg    'Xh is the distance of the aerodynamic center of the horizontal tail'