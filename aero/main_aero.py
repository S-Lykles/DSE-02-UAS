from compound_helicopter import *
from cl_cd import *

# Plot all configurations in 1 plot
def plots(cl,cd_dual,cd_comp):
    plt.plot(cd_dual,cl,label="Dual-Phase")
    plt.plot(cd_comp,cl,label="Compound Helicopter")
    plt.xlabel("cd")
    plt.ylabel("cl")
    plt.title("cd vs cl")
    plt.legend()
    plt.show()

# Call separate configurations
b = 6                                          # Width
S = 3.763                                      # Surface Area
cl_max = 1.47                                  # Cl_max
cl, cd_comp = dragpolar_comp(5,S)
cl, cd_dual = dragpolar_dual(b,S,500,43,0.6,2) # (b,S,h,v,c,S_fuselage)
print(cd0_fuselage)

def constants():
    
    return b, S, cl_max
 
plots(cl,cd_dual,cd_comp)