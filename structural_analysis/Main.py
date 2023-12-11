from structural_analysis.component_weight_estimation import (class_two_dual_phase,
                                                             class_two_tilt_wing, class_two_compound_helicopter)
from structural_analysis.weight_estimation_ii import (class_2_Cessna, class_2_Torenbeek, class_2_USAF, class_1_WE)

OEW_tilt_wing = class_two_compound_helicopter(MTOW= 160 *2.20462262)


print(OEW_tilt_wing)

print(class_2_Cessna(160, 'Tiltwing', 3.673, 8, 2, 15, 2.5, 0.8))
print(class_2_Cessna(160, 'Compound', 3.673, 8, 2, 15, 2.5, 0.8))
print(class_2_Torenbeek(160,'Tiltwing', 3.673, 6, 0, 2, 15, 2.5, 0.8, 0.07))
print(class_2_Torenbeek(160,'Compound', 3.673, 6, 0, 2, 15, 2.5, 0.8, 0.07))
print(class_2_USAF(160, 'Tiltwing', 3.673,8,0, 0.25, 0.12, 2, 15, 2.5, 0.7, 0.7, 65, 42))
print(class_2_USAF(160, 'Compound', 3.673,8,0, 0.25, 0.12, 2, 15, 2.5, 0.7, 0.7, 65, 42))

print(' class 1 weight estimation OEW in kg + mass fractions:', class_1_WE(160))