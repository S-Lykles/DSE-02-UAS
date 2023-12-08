from structural_analysis.component_weight_estimation import (class_2_tailsitter, class_two_dual_phase,
                                                             class_two_tilt_wing, class_two_compound_helicopter)

OEW_compound_helicopter = class_two_compound_helicopter(MTOW= 160 *2.20462262)
print(OEW_compound_helicopter)