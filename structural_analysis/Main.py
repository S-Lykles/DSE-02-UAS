from structural_analysis.component_weight_estimation import (class_2_tailsitter, class_two_dual_phase,
                                                             class_two_tilt_wing, class_two_compound_helicopter)


OEW_dual_phase = class_two_dual_phase(160, 1, 1, 12)
print(OEW_dual_phase)