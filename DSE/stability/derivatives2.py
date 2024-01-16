from DSE.aero.aero_constants import *
import numpy as np

# CYva
#S
#dsigdB?
#Vv = Sv * lv / (Sw * bw)
# V
# Sv
# CLav?
# tau_r
# br
# bv
# b
# lv
# CLva
# CLaw
# gammaw
# lambdaw
#ARw
#Zwf
# lfus
#wfus
#bw
#Sw
#Vw
#CLah
# gammah
# lambdah
#ARh
#Zhf
#bh
#Sh
# Vh
#CYBv
# zv
#CLw0
# CLh0
#CYpv
# Vr
#CDw0
# CDh0
#CL0
# Clda

#Lateral stability derivatives
#CYB
CYB = -CYva * (1-dsigdB) * (Vv/V)**2 * Sv/S
CYBdot = 0
CYp = (-8 / (3 * np.pi)) * (Vv / V)**2 * (bv * Sv) / (b * S) * CLav
CYr = 2 * CYva * ((Sv * lv) / (S * b)) * (Vv / V)**2
CYda = 0

CYdr = CLva * (Sv / S) * tau_r * br / bv

ClB = (- 0.25 * CLaw * gammaw * (2/3) * (1 + 2 * lambdaw) / (1 + lambdaw) - 1.2 * (ARw**0.5 * Zwf * (lfus + wfus)/bw**2)) * (Sw * bw / (S * b)) * (Vw / V)**2 \\
    (- 0.25 * CLah * gammah * (2/3) * (1 + 2 * lambdah) / (1 + lambdah) - 1.2 * (ARh**0.5 * Zhf * (lfus + wfus)/bh**2)) * (Sh * bh / (S * b)) * (Vh / V)**2+ CYBv * zv / b

#Clp

Clr = (zv / b) * CYr + 0.25 * (CLw0 * Sw * bw / (S * b) + ((CLh0 * Sh * bh) / (S*b)) * (Vh / V)**2)

#Clda

Cldr = (zv / b) * CLva * (Sv / S) * tau_r * br / bv

#CnB
CnBdot = 0
Cnp = (-lv / b)* CYpv - (1/8) * (CLw0 * Sw * bw / (S*b) + ((CLh0 * Sh * bh) / (S * b)) * (Vr / V)**2)
Cnr = - (lv / b) * CYr - 0.25 * (CDw0 * Sw * bw / (S * b) + CDh0 * Sh * bh * (Vh / V)**2 / (S * b))
Cnda = -0.2 * CL0 * Clda
Cndr = -CLva * (Sv * lv / (S * b)) * tau_r * br / bv

