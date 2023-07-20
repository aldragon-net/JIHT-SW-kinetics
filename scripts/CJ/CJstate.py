import cantera as ct

from sdtoolbox.postshock import CJspeed, PostShock_eq
from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr

# Initial state specification:
# P1 = Initial Pressure
# T1 = Initial Temperature
# U = Shock Speed
# q = Initial Composition
# mech = Cantera mechanism File name

# stoichiometric H2-air detonation at nominal atmospheric conditions
P1 = 30000.
T1 = 293
q = 'C2H2:20 O2:10 AR:70'
mech = '../cti/CRECK_sink.yaml'
gas_initial = ct.Solution(mech)
gas_initial.TPX = T1, P1, q
rho_1 = gas_initial.density

# compute CJ speed
[cj_speed, R2, plot_data] = CJspeed(P1, T1, q, mech, fullOutput=True)

# compute equilibrium CJ state parameters
gas = PostShock_eq(cj_speed, P1, T1, q, mech)
ae = soundspeed_eq(gas)
af = soundspeed_fr(gas)
rho_2 = gas.density
gammae = ae**2*rho_2/gas.P
gammaf = af**2*rho_2/gas.P
w2 = cj_speed*rho_1/rho_2
u2 = cj_speed-w2
print('CJ computation for ' + mech + ' with composition ' + q )
print('Initial conditions: P1 = %.3e Pa & T1 = %.2f K'  % (P1,T1)  )
print('CJ Speed   %.1f m/s' % cj_speed)
print('CJ State')
print('   Pressure   %.3e Pa' % gas.P)
print('   Temperature  %.1f K' % gas.T)
print('   Density  %.3f kg/m3' % gas.density)
print('   Entropy  %.3e J/K' % gas.entropy_mass)
print('   w2 (wave frame) %.1f m/s' % w2)
print('   u2 (lab frame) %.1f m/s' % u2)
print('   c2 frozen %.1f m/s' % af)
print('   c2 equilbrium %.1f m/s' % ae)
print('   gamma2 frozen %.3f ' % gammaf)
print('   gamma2 equilbrium %.3f ' % gammae)


def pressure_CJ(pressures, temperature, mixture, mech):
    cj_speeds = []
    states = []
    for pressure in pressures:
        cj_speed, _, _ = CJspeed(pressure, temperature, mixture, mech)
        gas = PostShock_eq(cj_speed, P1, T1, q, mech)
        cj_speeds.append(cj_speed)
        states.append(gas)
    return cj_speeds, states


def write_output(filename: str, x: list, ys: list[list], x_label: str, y_labels: list[str]):
    path = 'output/' + filename
    with open(path, "w") as output:
        labels = [x_label].extend(y_labels).append('\n')
        output.write(labels)
        for i in range(len(x)):
            values = [x[i]].extend(ys[i])
            output.write(', '.join(values))
            output.write('\n')
