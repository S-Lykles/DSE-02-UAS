import numpy as np
from matplotlib import pyplot as plt
import const
from powertrain_weight_calculator import table_hybrid_propulsion_weights, E_rho_bat, E_rho_H, spec_tank_W, SFC_from_Pmax
from aero.cl_cd import dragpolar_dual
from power_curves.rotor_tool import rotor_sizing_tool, P_profile_drag, P_induced, delta_p_climb
from power_curves.mass_frac import fuel_weight

DL = 500
N = 4
b = 6
S = 3.763
k_dl = 1.01 # lower than normal because the rotors are in free air
vc = 1
Ploss_frac = 0.05

# All in kg/J
# PLACEHOLDER VALUES for all the fuel engines
SFC_2stroke = 320 / (1000 * 1000 * 3600)
SFC_4stroke = 320 / (1000 * 1000 * 3600)
SFC_rotarty = 320 / (1000 * 1000 * 3600)
SFC_turboshaft = 320 / (1000 * 1000 * 3600)
SFC_hydrogen = 1 / (E_rho_H*1e6/spec_tank_W)
SFC_battery = 1 / (E_rho_bat*1e6)


def P_max(DL, N):
    R, D_v, omega, T_level, sig_max = rotor_sizing_tool(const.MTOW, DL, N, const.v_cruise)
    v = 0
    P_p = P_profile_drag(v, const.MTOW, N, R, omega, sig_max)
    P_i = P_induced(0,DL,const.MTOW, k_dl=k_dl)
    dP_c = delta_p_climb(vc, const.MTOW)
    P_max = (P_p + P_i + dP_c) * (1 + Ploss_frac) + const.P_aux + const.P_pay_end
    return P_max

def P_cruise(b,S):
    CL, CD = dragpolar_dual(b,S,CL_start=0.2,CL_end=1.2,CL_step=1000)
    v_cruise = const.v_cruise * 1.1 # A margin on the minimum cruise speed
    v1 = np.sqrt(const.MTOW * 2 / (const.rho0 * S * CL))

    D = const.MTOW * CD / CL
    P = D * v1 + const.P_aux + const.P_pay_end
    E = const.R_cruise / v1 * P
    idx = np.argmin(E[np.where(v1 > v_cruise)])
    return P[idx]


if __name__ == '__main__':
    CL, CD = dragpolar_dual(b,S,CL_start=0.1,CL_end=1.2,CL_step=1000)
    P_m = P_max(DL, N)/1000
    df = table_hybrid_propulsion_weights(P_cruise(b,S)/1000, P_m, 10*60)
    df.loc['Battery'] = ['-'] * len(df.columns)

    SFC = SFC_from_Pmax(P_m,800) +  [SFC_hydrogen, SFC_hydrogen, SFC_battery]
    opts = ['2-Stroke', '4-Stroke', 'Rotary', 'Turboshaft', 'Liquid Hydrogen Current Tech', 'Liquid Hydrogen Future Tech', 'Battery']
    e_frac = [1,1,1,1,1/spec_tank_W,1/spec_tank_W,0]
    m_f_e = []
    m_f_p = []
    for o, e, s in zip(opts, e_frac, SFC):
        m_f_e.append(fuel_weight(CL,CD,s,S,energy_fraction=e,which='endurance'))
        m_f_p.append(fuel_weight(CL,CD,s,S,energy_fraction=e,which='payload'))

    df['Fuel Mass Endurance'] = m_f_e
    df['Fuel Mass Payload'] = m_f_p
    print(df.to_latex())
    print(df)
