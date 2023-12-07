from structural_analysis.component_weight_estimation import (class_2_tailsitter, class_two_dual_phase,
                                                             class_two_tilt_wing, class_two_compound_helicopter)
from parameters_weight_estimations import MTOW, l_sm, l_sn, W_rotor


OEW_dual_phase = class_two_dual_phase(MTOW, l_sm, l_sn, W_rotor)
print(OEW_dual_phase)