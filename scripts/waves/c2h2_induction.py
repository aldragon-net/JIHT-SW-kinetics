
import cantera as ct
from sdtoolbox.thermo import soundspeed_fr, soundspeed_eq
from sdtoolbox.postshock import PostShock_fr, PostShock_eq,CJspeed
import math
import numpy as np
from scipy.interpolate import pchip # creates a PCHIP interpolation function that can then be given a query point
from scipy.optimize import fminbound
import matplotlib.pyplot as plt

## set initial state, composition, and gas object for driven section
P1 = 100000.
T1 = 300.
q1 = 'CH4:20 AR:80' 
driven_mech = 'mechs/gri30_highT.yaml'
driven_gas = ct.Solution(driven_mech)
driven_gas.TPX = T1,P1,q1

## Evaluate initial state 
rho1 = driven_gas.density
a1 = soundspeed_fr(driven_gas)

## Evaluate post-shock state (frozen) for a range of shock speeds 
print('Generating points on shock P-u curve')
Ustart = 1100
Ustop = 1500
Ustep = 30
nsteps = int((Ustop-Ustart)/Ustep)

a2 = []; P2 = []; T2 = []; rho2 = []; u2 = []; Us = []

for U in np.linspace(Ustart, Ustop, num=nsteps):
    shocked_gas = PostShock_fr(U, P1, T1, q1, driven_mech)
    a2.append(soundspeed_fr(shocked_gas))
    P2.append(shocked_gas.P)
    T2.append(shocked_gas.T)
    rho2.append(shocked_gas.density)
    w2 = rho1*U/shocked_gas.density
    u2.append(U - w2)
    Us.append(U)


def c2h2_induction_time(T: float, P: float, carbon_fraction: float) -> float:
    """from Eremin2012"""
    high = {'a': 2.194E-9, 'E': 27450}
    low = {'a': 1.0803E-4, 'E': 9375}

    def t_ind(a, E):
        return (a*T*math.exp(E/T)) / (P*carbon_fraction)

    t_low = t_ind(low['a'], low['E'])
    t_high = t_ind(high['a'], high['E'])
    return min(t_low, t_high)


carbon_fraction = 1

for i, u in enumerate(Us):
    t = c2h2_induction_time(T2[i], P2[i], carbon_fraction)
    print(f'u = {u:.0f}, T2 = {T2[i]:.0f}, P2= {P2[i]:.0f}, time={t:.3e}')



# Sort all State 2 arrays by P2; State 3 arrays by P3.
iP2 = P2.argsort()
iP3 = P3.argsort()

P2_u2_fun = pchip(P2[iP2],u2[iP2])
u_max = P2_u2_fun(P4)

pmin = 1.01
pmax = max(P2)/P1

# define intersection function
P2P1_u2_fun = pchip(P2[iP2]/P1,u2[iP2])
P3P1_u3_fun = pchip(P3[iP3]/P1,u3[iP3])
fun = lambda x: (P2P1_u2_fun(x)-P3P1_u3_fun(x))**2

# solve for intersection of P-u curves to find pressure at contact surface
P_contact = fminbound(fun, pmin, pmax)

# find shock speed
P2P1_Us_fun = pchip(P2[iP2]/P1,Us[iP2])
U_shock = P2P1_Us_fun(P_contact)
M_shock = U_shock/a1

# find velocity and states on driven side of contact surface
u2_contact = P2P1_u2_fun(P_contact)
P2P1_T2_fun = pchip(P2[iP2]/P1,T2[iP2])
T2_contact = P2P1_T2_fun(P_contact)
P2P1_a2_fun = pchip(P2[iP2]/P1,a2[iP2])
a2_contact = P2P1_a2_fun(P_contact)
P2P1_rho2_fun = pchip(P2[iP2]/P1,rho2[iP2])
rho2_contact = P2P1_rho2_fun(P_contact)
M2 = u2_contact/a2_contact

# find velocity and states on driver side of contact surface
u3_contact = P3P1_u3_fun(P_contact)
P3P1_T3_fun = pchip(P3[iP3]/P1,T3[iP3])
T3_contact = P3P1_T3_fun(P_contact)
P3P1_a3_fun = pchip(P3[iP3]/P1,a3[iP3])
a3_contact = P3P1_a3_fun(P_contact)
P3P1_rho3_fun = pchip(P3[iP3]/P1,rho3[iP3])
rho3_contact = P3P1_rho3_fun(P_contact)
M3 = u3_contact/a3_contact

##
# Display results on screen
print('Driven State (1)')
print('   Fill gas '+q1)
print('   Pressure '+str(P1/1e3)+' (kPa)')
print('   Temperature '+str(T1)+' (K)')
print('   Sound speed '+str(a1)+' (m/s)')
print('   Density '+str(rho1)+' (kg/m3)')

if CASE_DRIVER=='cj':
    print('Driver State (CJ)')
    print('   Fill gas '+q4)
    print('   Fill Pressure '+str(P_driver_fill/1E3)+' (kPa)')
    print('   Fill Temperature '+str(T_driver_fill)+' (K)')
    print('   CJ speed '+str(cj_speed)+' (m/s)')
elif CASE_DRIVER=='uv':
    print('Driver State (UV)')
    print('   Fill gas '+q4)
    print('   Fill Pressure '+str(P_driver_fill/1E3)+' (kPa)')
    print('   Fill Temperature '+str(T_driver_fill)+' (K)')
elif CASE_DRIVER=='gas':
    print('Driver State (4)')
    print('   Fill gas '+q4)

print('   Pressure '+str(P4/1e3)+' (kPa)')
print('   Temperature '+str(T4)+' (K)')
print('   Sound speed '+str(a3[0])+' (m/s)')
print('   Density '+str(rho3[0])+' (kg/m3)')
print('Solution matching P-u for states 2 & 3')
print('Shock speed '+str(U_shock)+' (m/s)')
print('Shock Mach number '+str(M_shock)+' ')
print('Shock Pressure ratio '+str(P_contact)+' ')
print('Postshock state (2) in driven section')
print('   Pressure '+str(P1*P_contact/1e3)+' (kPa)')
print('   Velocity '+str(u2_contact)+' (m/s)')
print('   Temperature '+str(T2_contact)+' (K)')
print('   Sound speed '+str(a2_contact)+' (m/s)')
print('   Density '+str(rho2_contact)+' (kg/m3)')
print('   Mach number '+str(M2)+' ')
print('Expansion state (3) in driver section')
print('   Pressure '+str(P1*P_contact/1e3)+' (kPa)')
print('   Velocity '+str(u3_contact)+' (m/s)')
print('   Temperature '+str(T3_contact)+' (K)')
print('   Sound speed '+str(a3_contact)+' (m/s)')
print('   Density '+str(rho3_contact)+' (kg/m3)')
print('   Mach number '+str(M3)+' ')

Z = rho2_contact*a2_contact/(rho3_contact*a3_contact)

print('Impedance ratio at contact surface (rho a)_2/(rho a)_3 '+str(Z))

##
# plot pressure-velocity diagram
plt.figure(1)

plt.plot(u2,P2/P1,linewidth=2) 
plt.plot(u3,P3/P1,linewidth=2)
plt.title('Shocktube P-u relations',fontsize=12)
plt.xlabel('Velocity (m/s)',fontsize=12)
plt.ylabel('Pressure P/P1',fontsize=12)
plt.axis([u3[0], max(u2), 0, P4/P1])
plt.legend(['Driven shock','Driver expansion'],loc='upper center')
plt.show()
