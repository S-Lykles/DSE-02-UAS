from .compound_helicopter import *
from .cl_cd import *
from .titl_wing import *
from plot_setting import *
import const


# Plot all configurations in 1 plot
def plots(name=None):
    plt.rcParams.update(tex_fonts)
    size = set_size(fraction=1, subplots=(1,1))
    plt.figure(figsize=(size[0]*0.7,size[1]))
    plt.plot(cd_dual,cl,label="Dual-Phase")
    plt.plot(cd_comp,cl,label="Compound Helicopter")
    plt.plot(cd_tilt_wing,cl,label="Tilt Wing")
    plt.gca().grid(which='major', color='#DDDDDD', linewidth=0.8)
    plt.gca().grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
    plt.minorticks_on()
    plt.xlabel("$C_D$")
    plt.ylabel("$C_L$")
    plt.title("$C_D$ vs $C_L$")
    plt.tight_layout()
    plt.legend()
    
    if name is not None:
        plt.savefig('aero/'+name)
    plt.show()


b = 6                                          # Width
S = 3.763                                      # Surface Area
h = 500
v = const.v_cruise
c = 0.5
S_f = 2
d_eng = 0.5
N_eng = 4
Lambda = 0.45

# Call separate configurations
cl, cd_comp = dragpolar_comp(b-1,S,d_eng,2,Lambda,v,h,c=None,Sf=2)
cl, cd_dual = dragpolar_dual(b,S,h,v,c,S_f)                                 # (b,S,h,v,c,S_fuselage)
cl, cd_tilt_wing, cd_prop = dragpolar_tilt_wing(b-1,S,h,v,c,S_f,d_eng,N_eng,Lambda)
    
def find_cd0():    
    print("cd0 dual",min(cd_dual))
    print("cd0 comp",min(cd_comp))
    print("cd0 tilt",min(cd_tilt_wing))
    return

# Find max L_D
def maxL_D():
    for i, data in zip([cd_comp, cd_dual, cd_tilt_wing], ['Compound', 'Dual Phase', 'Tilt Wing']):
        L_D = cl / i
        max_L_D = max(L_D)
        print(f"L/D max for {data} is {max_L_D}")
    return

find_cd0()
maxL_D()
plots('cl_cd.pdf')

# def constants_aero():
#     b = 6
#     S = 3.763
#     h = 500
#     v = 43
#     c = 0.5
#     S_f = 2
#     d_eng = 0.5
#     N_eng = 4
#     Lambda = 0.4
#     cl_max_dual = 1.47
#     return b, S, cl_max_dual
 

# print(constants_aero())