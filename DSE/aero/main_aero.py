from DSE.aero.cl_cd import *
from DSE.plot_setting import *
from DSE import const
from DSE.plot_setting import report_tex, set_size


# Plot all configurations in 1 plot
def plots(name=None):
    plt.rcParams.update(report_tex)
    size = set_size(fraction=1, subplots=(1,1))
    plt.figure(figsize=(size[0]*0.7,size[1]))
    plt.plot(cd_dual,cl,label="Dual-Phase")
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
cl, cd_dual = dragpolar_dual(b,S,h,v,c,S_f)                                 # (b,S,h,v,c,S_fuselage)
    
def find_cd0():    
    print("cd0 dual",min(cd_dual))
    return

# Find max L_D
def maxL_D():
    for i, data in zip([cd_dual], ['Dual Phase']):
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