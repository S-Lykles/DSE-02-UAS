import numpy as np
from pathlib import Path
from DSE import const
from DSE.power_curves import rotor_tool
from matplotlib import pyplot as plt

file_dir = Path(__file__).parent
data = np.genfromtxt(file_dir / r"clarky.csv", delimiter=",", skip_header=1)
alpha_clarky = data[:, 0]
cl_clarky = data[:, 1]
cd_clarky = data[:, 2]
Cl_func_clarky = lambda a: np.interp(a, alpha_clarky, cl_clarky, left=0, right=0)
Cd_func_clarky = lambda a: np.interp(a, alpha_clarky, cd_clarky)

# r = np.linspace(0, 0.5, 1000)
# omega = (5000*2*np.pi) / 60  # 5000 rpm

def chord_dist(r, a=0.03, b=-0.02):
    return a + b*r

def solidity(R, r, chord, N=2):
    return N * np.trapz(chord, r) / (np.pi * R**2)

def twist_dist(r, a=14, b=-22):
    return a + b*r

def dT_vi(r, vi, rho=const.rho0):
    return 4 * np.pi * rho * vi**2 * r

def dT_dr_bem(r, omega, chord, twist, vi, N=2, Vc=0, rho=const.rho0, Cl_func=Cl_func_clarky):
    alpha = twist - np.rad2deg((Vc + vi) / (omega * r))
    return 0.5 * rho * (omega * r)**2 * chord * Cl_func(alpha) * N

def vi_bem(r, dT_dr, rho=const.rho0):
    return np.sqrt(dT_dr / (4 * np.pi * rho * r))

def solve_dT_dr(r, omega, chord, twist, N=2, Cl_func=Cl_func_clarky, max_iter=100, tol=1e-6):
    """twist in deg"""
    vi = np.zeros_like(r)

    for _ in range(max_iter):
        dT1 = dT_dr_bem(r, omega, chord, twist, vi, N=N, Cl_func=Cl_func)
        vi = vi_bem(r, dT1)
        dT2 = dT_dr_bem(r, omega, chord, twist, vi, N=N, Cl_func=Cl_func)
        if np.linalg.norm(dT1 - dT2)/len(r) < tol:
            break
        dT1 = dT2

    return dT2

def dD_bem(r, dT, omega, chord, twist, N=2, Vc=0, rho=const.rho0, Cd_func=Cd_func_clarky):
    vi = vi_bem(r, dT)
    alpha = twist - np.rad2deg((Vc + vi) / (omega * r))
    dDp = 0.5 * const.rho0 * (omega*r)**2 * chord * Cd_func(alpha) * N
    return (dT * (Vc + vi) / (omega * r) + dDp)

if __name__ == "__main__":
    rpm_range = np.arange(1000, 6100, 100)

    rpm_values = np.zeros_like(rpm_range, dtype=float)
    total_thrust_values = np.zeros_like(rpm_range, dtype=float)

    target_rpm = 4000

    for i, rpm in enumerate(rpm_range):
        R = 0.5
        r = np.linspace(0.1, R, 1000)
        N = 4
        omega = (rpm*2*np.pi) / 60
        chord = chord_dist(r)
        twist = twist_dist(r)

        dT = solve_dT_dr(r, omega, chord, twist, N=N)
        vi = vi_bem(r, dT)
        dD = dD_bem(r, dT, omega, chord, twist, N=N)
        dQ = dD * r

        # to account for tip losses we integrate only to 0.97 R
        T = np.trapz(np.where(r<=0.97*R, dT, 0), r)
        T_ideal = np.trapz(dT, r)

        Q = np.trapz(dQ, r)
        P = Q * omega

        T_total = np.trapz(np.where(r <= 0.97 * R, dT, 0), r)
        rpm_values[i] = rpm
        total_thrust_values[i] = T_total

        if rpm == target_rpm:
            plt.plot(r, dT)
            plt.fill_between(r, 0, dT, alpha=0.2)
            plt.title(f"Thrust distribution at rpm = {rpm}")
            plt.xlabel('Radius [m]')
            plt.ylabel('Thurst [N]')
            plt.grid()
            plt.show()

            plt.plot(r, vi)
            plt.title(f'Induced velocity distribution at rpm = {rpm}')
            plt.xlabel('Radius [m]')
            plt.ylabel('Indiced velocity [m/s]')
            plt.grid()
            plt.show()

            plt.plot(r, dQ)
            plt.title(f'Torque distribution at rpm = {rpm}')
            plt.xlabel('Radius [m]')
            plt.ylabel('Torque [N/m]')
            plt.grid()
            plt.show()

            plt.plot(r, twist - np.rad2deg((vi) / (omega * r)))
            plt.title(f'Effective angle of attack distribution at rpm = {rpm}')
            plt.xlabel('Radius [m]')
            plt.ylabel('Effective angle of attack [deg]')
            plt.grid()
            plt.show()


            print(f"\nResults for RPM = {rpm}:")
            print(f"Tip loss % {100 * (T_ideal - T) / T_ideal:.2f}")
            print(f"Thrust: {T:.2f} N")
            print(f"Disk loading: {T / (np.pi * R ** 2):.2f} N/m^2")
            print(f"Torque: {Q:.2f} Nm")
            print(f"Power: {P:.2f} W")

            plt.clf()

    plt.plot(rpm_values, total_thrust_values)
    plt.title('Total Thrust vs RPM')
    plt.xlabel('RPM')
    plt.ylabel('Total Thrust (N)')
    plt.grid()
    plt.show()

