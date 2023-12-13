import numpy as np

from structural_analysis.component_weight_estimation import (class_two_dual_phase,
                                                             class_two_tilt_wing, class_two_compound_helicopter)
from structural_analysis.weight_estimation_ii import (class_2_Cessna, class_2_Torenbeek, class_2_USAF, class_1_WE)
from parameters_weight_estimations import b
# OEW_tilt_wing = class_two_dual_phase(MTOW= 160 * 2.20462262)
#
#
# print(OEW_tilt_wing)



print(' class 1 weight estimation OEW in kg + mass fractions:', class_1_WE(160))

print ('Wstruc (lbm)', 'Wstruc (kg)', 'Wwing (kg)', 'Wfus (kg)', 'Wnac (kg)', 'Wlg (kg)', 'Wemp (kg)')
print(class_2_Cessna(150, 'Dual Phase', 3.2, 6**2/3.2, 2, 21.3, 2.15, 0.8))
print(class_2_Cessna(150, 'Compound', 3.2, 5**2/3.2, 2, 21.3, 2.15, 0.8))
print(class_2_Cessna(150, 'Tiltwing', 3.2, 5**2/3.2, 2, 21.3, 2.15, 0.8))

# print(class_2_Torenbeek(160,'Dual Phase', 3.763, 6, 0, 2, 36, 2.5, 0.8, 0.2234))
# print(class_2_Torenbeek(160,'Compound', 3.763, 5, 0, 2, 44, 2.5, 0.8, 0.1862))
# print(class_2_Torenbeek(160,'Tiltwing', 3.763, 5, 0, 2, 40, 2.5, 0.8, 0.1862))
#
# print(class_2_USAF(160, 'Dual Phase', 3.763, 9.567, 0, 0.45, 0.12, 2, 36, 2.5, 0.7, 0.7, 65, 42))
# print(class_2_USAF(160, 'Compound', 3.763, 6.644, 0, 0.45, 0.12, 2, 44, 2.5, 0.7, 0.7, 65, 42))
# print(class_2_USAF(160, 'Tiltwing', 3.763, 6.644, 0, 0.45, 0.12, 2, 40, 2.5, 0.7, 0.7, 65, 42))
