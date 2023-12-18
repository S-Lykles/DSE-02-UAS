from tilt_wing.powertrain import P_max
from dual_phase.powertrain import P_cruise, CL_max
from power_curves.mass_frac import fuel_weight
from aero.cl_cd import dragpolar_dual

for (DL, N) in [(100, 1), (200, 1), (300, 2), (300, 3), (400, 1), (400, 4)]:
    print(P_max(DL, N, k_dl=1.04))
    print()

for b, S in [(5, 3.76), (4.5, 3.5), (4, 3.2), (3, 2.7)]:
    print(P_cruise(b, S))
    print()
    
for b, S in [(5, 3.76), (4.5, 3.5), (4, 3.2), (3, 2.7)]:
    CL, CD = dragpolar_dual(b,S,CL_start=0.2,CL_end=CL_max,CL_step=1000)
    
    SFC = 320 / (1000 * 1000 * 3600)
    mf = fuel_weight(CL, CD, SFC, S)
    print(mf)
    print()