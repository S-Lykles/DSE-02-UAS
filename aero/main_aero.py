from compound_helicopter import *
from cl_cd import *

# Plot all configurations in 1 plot
def plots(cl,cd,cd_heli):
    plt.plot(cd,cl,label="Dual-Phase")
    plt.plot(cd_heli,cl,label="Compound Helicopter")
    plt.xlabel("cd")
    plt.ylabel("cl")
    plt.title("cd vs cl")
    plt.legend()
    plt.show()

# Call separate configurations
cl, cd_heli = dragpolar_heli(6,2)
cl,cd, cd0_wing = dragpolar_dual(6,3.763,500,43,0.6,2)
plots(cl,cd,cd_heli)