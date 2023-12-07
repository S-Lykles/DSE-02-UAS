from compound_helicopter import *
from cl_cd import *
from titl_wing import *

# Plot all configurations in 1 plot
def plots():
    plt.plot(cd_dual,cl,label="Dual-Phase")
    plt.plot(cd_comp,cl,label="Compound Helicopter")
    plt.plot(cd_tilt_wing,cl,label="Tilt Wing")
    plt.xlabel("cd")
    plt.ylabel("cl")
    plt.title("cd vs cl")
    plt.legend()
    plt.show()

# Call separate configurations
b = 6                                          # Width
S = 3.763                                      # Surface Area
h = 500
v = 43
c = 0.5
S_f = 2
d_eng = 0.5
N_eng = 4
Lambda = 0.4
c_root = 0.5
                                  # Cl_max
cl, cd_comp = dragpolar_comp(5,S)
cl, cd_dual = dragpolar_dual(b,S,h,v,c,S_f)     # (b,S,h,v,c,S_fuselage)
cl, cd_tilt_wing = dragpolar_tilt_wing(b-1,S,h,v,c,S_f,d_eng,N_eng,Lambda,c_root)



def constants_aero():
    b = 6
    S = 3.763
    h = 500
    v = 43
    c = 0.5
    S_f = 2
    d_eng = 0.5
    N_eng = 4
    Lambda = 0.4
    c_root = 0.5
    cl_max_dual = 1.47
    return b, S, cl_max_dual
 
plots()
print(constants_aero())