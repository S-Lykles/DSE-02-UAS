import numpy as np

#Unit conversions
kw_to_hp = 1.341022
kg_to_lb = 2.20462262
m_to_ft = 3.28084
mps_to_kts = 1.94384

#Class I Weight Estimation methods, various found
#inputs in SI: MTOW (Maximum takeoff weight, kg)
def class_1_WE(MTOW):
    WF_UAV_Market = 0.699 * (MTOW ** -0.051)
    WF_UAV_RaymS = 1.53 * (MTOW ** -0.16)
    WF_UAV_RaymTaccRecc = 0.86 * (MTOW ** -0.06)
    WF_UAV_RaymHA = 2.48 * (MTOW ** -0.18)

    OEW_Market = MTOW * WF_UAV_Market
    OEW_RaymS = MTOW * WF_UAV_RaymS
    OEW_TaccRecc = MTOW * WF_UAV_RaymTaccRecc
    OEW_RaymHA = MTOW * WF_UAV_RaymHA

    return MTOW, OEW_Market, OEW_RaymS, OEW_TaccRecc, OEW_RaymHA, WF_UAV_Market, WF_UAV_RaymS, WF_UAV_RaymTaccRecc, WF_UAV_RaymHA

#Class II Weight estimation General Aviation Roskam Cessna method:
#inputs IN SI: MTOW (Maximum Take Off Weight, kg), config (to select between configuration), S (Wing Surface Area, m^2),
#              A (Aspect Ratio), n_ult (ultimate load factor), Pmax (maximum cruise power, kW),
#              lfus (fuselage length, m), dfus (fuselage diameter, m)
def class_2_Cessna(MTOW, config, S, A, n_ult, Pmax, lfus, dfus):
    if config == 'Tiltwing':
        facwing = 2
    else:
        facwing = 1

    #Unit Conversion
    MTOW *= kg_to_lb
    S *= (m_to_ft ** 2)
    lfus *= m_to_ft
    perfus = dfus * m_to_ft * np.pi
    Pmax *= kw_to_hp

    # Weight estimation of main wing in kg(full cantilever wing assumed)

    W_wing = facwing * (0.04674 * MTOW ** 0.397 * S ** 0.36 * n_ult ** 0.397 * A ** 1.712)

    # Weight estimation of fuselage (high wing configuration assumed)
    W_f = (14.86 * (MTOW ** 0.144) * ((lfus / perfus) ** 0.778) * (lfus ** 0.383))

    # Weight estimation of the nacelle (horizontally opposed engines assumed)
    W_nac = 0.24 * Pmax

    # Weight estimation of a non-retractable landing gear using a comparable design
    W_lg = 0.02 * MTOW

    #Weight estimation of empennage using a comparable design
    W_emp = 0.023 * MTOW

    # Summation to generate a final estimation for the structural weight
    W_strucIMP = W_wing + W_f + W_nac + W_lg + W_emp
    W_strucSI = W_strucIMP / kg_to_lb

    return 'Cessna', config, W_strucIMP, W_strucSI, W_wing / kg_to_lb, W_f/ kg_to_lb, W_nac / kg_to_lb, W_lg / kg_to_lb, W_emp / kg_to_lb

#Class II Weight estimation General Aviation Roskam USAF method:
#inputs IN SI: MTOW (Maximum Take Off Weight, kg), config (to select between configuration), S (Wing Surface Area, m^2),
#              A (Aspect Ratio), sweep (sweep angle at 1/4 chord, rad), taper (taper ratio), tcm (maximum thickness/chord ratio),
#              n_ult (ultimate load factor),  Pmax (maximum cruise power, kW), lf (fuselage length, m), wf (fuselage width, m),
#              hf (fuselage height, m), Vmax (maximum flight speed, m/s), Vcruise (cruise speed, m/s)
def class_2_USAF(MTOW, config, S, A, sweep, taper, tcm, n_ult, Pmax, lf, wf, hf, Vmax, Vcruise):
    if config == 'Tiltwing':
        facwing = 2
    else:
        facwing = 1

    # Unit Conversion
    MTOW *= kg_to_lb
    S *= (m_to_ft ** 2)
    Pmax *= kw_to_hp
    lf *= m_to_ft
    wf *= m_to_ft
    hf *= m_to_ft
    Vmax *= mps_to_kts
    Vcruise *= mps_to_kts

    # Weight estimation of main wing in lbm (full cantilever wing assumed)

    W_wing = facwing * 96.948 * ((((MTOW * n_ult / 10000) ** 0.65) * ((A/np.cos(sweep)) ** 0.57) * ((0.01*S) ** 0.61) * (((1 + taper)/(2 * tcm)) ** 0.36) * ((1 + (Vmax/500)) ** 0.5)) ** 0.993)

    # Weight estimation of fuselage in lbm (high wing configuration assumed)
    W_f = 200*(((MTOW * n_ult/10000) ** 0.286 * ((0.1 * lf) ** 0.857) * (0.1*(wf + hf) * ((0.01 * Vcruise) ** 0.338))) ** 1.1)

    # Weight estimation of the nacelle in lbm from CESSNA (horizontally opposed engines assumed)
    W_nac = 0.24 * Pmax

    # Weight estimation of a non-retractable landing gear using a comparable design in lbm
    W_lg = 0.02 * MTOW

    # Weight estimation of empennage using a comparable design in lbm
    W_emp = 0.023 * MTOW

    # Summation to generate a final estimation for the structural weight in lbm
    W_strucIMP = W_wing + W_f + W_nac + W_lg + W_emp
    W_strucSI = W_strucIMP / kg_to_lb

    return 'USAF', config, W_strucIMP, W_strucSI, W_wing / kg_to_lb, W_f/ kg_to_lb, W_nac / kg_to_lb, W_lg / kg_to_lb, W_emp / kg_to_lb

#Class II Weight estimation General Aviation Roskam Torenbeek method:
#inputs IN SI: MTOW (Maximum Take Off Weight, kg), config (to select between configuration), S (Wing Surface Area, m^2),
#              b (span, m), sweep (sweep angle at 1/2 chord, rad),  n_ult (ultimate load factor),
#              Pmax (maximum cruise power, kW), lfus (fuselage length, m), dfus (fuselage diameter, m), tr (max thickness at root, m)
def class_2_Torenbeek(MTOW, config, S, b, sweep, n_ult, Pmax, lfus, dfus, tr):
    if config == 'Tiltwing':
        facwing = 2
    else:
        facwing = 1

    # Unit Conversion
    MTOW *= kg_to_lb
    S *= (m_to_ft ** 2)
    b *= m_to_ft
    Pmax *= kw_to_hp
    lfus *= m_to_ft
    perfus = dfus * m_to_ft * np.pi
    tr *= m_to_ft


    # Weight estimation of main wing in lb
    W_wing = facwing * (0.00125 * MTOW*((b/np.cos(sweep)) ** 0.75) * ((1 + 6.3 *np.cos(sweep)/b) ** 0.5) * (n_ult ** 0.55) * ((b * S / (tr * MTOW * np.cos(sweep))) ** 0.3))

    # Weight estimation of fuselage from CESSNA (high wing configuration assumed)
    W_f = (14.86 * (MTOW ** 0.144) * ((lfus / perfus) ** 0.778) * (lfus ** 0.383))

    # Weight estimation of the nacelle (turboprop)
    W_nac = 0.14 * Pmax

    # Weight estimation of a non-retractable landing gear using a comparable design
    W_lg = 0.02 * MTOW

    # Weight estimation of empennage
    W_emp = 0.023 * MTOW

    # Summation to generate a final estimation for the structural weight
    W_strucIMP = W_wing + W_f + W_nac + W_lg + W_emp
    W_strucSI = W_strucIMP / kg_to_lb

    return 'Torenbeek', config, W_strucIMP, W_strucSI, W_wing / kg_to_lb, W_f/ kg_to_lb, W_nac / kg_to_lb, W_lg / kg_to_lb, W_emp / kg_to_lb




 # Rotor Weight Estimation
def Propsizing():
    rot = 1