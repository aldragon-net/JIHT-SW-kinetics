from sdtoolbox.postshock import CJspeed, PostShock_eq

# Initial state specification:
# P1 = Initial Pressure
# T1 = Initial Temperature
# U = Shock Speed
# q = Initial Composition
# mech = Cantera mechanism File name

# # stoichiometric H2-air detonation at nominal atmospheric conditions
# P1 = 30000.
# T1 = 293
# q = 'C2H2:20 O2:10 AR:70'
# mech = '../cti/CRECK_sink.yaml'
# gas_initial = ct.Solution(mech)
# gas_initial.TPX = T1, P1, q
# rho_1 = gas_initial.density

# # compute CJ speed
# [cj_speed, R2, plot_data] = CJspeed(P1, T1, q, mech, fullOutput=True)

# # compute equilibrium CJ state parameters
# gas = PostShock_eq(cj_speed, P1, T1, q, mech)
# ae = soundspeed_eq(gas)
# af = soundspeed_fr(gas)
# rho_2 = gas.density
# gammae = ae**2*rho_2/gas.P
# gammaf = af**2*rho_2/gas.P
# w2 = cj_speed*rho_1/rho_2
# u2 = cj_speed-w2
# print('CJ computation for ' + mech + ' with composition ' + q )
# print('Initial conditions: P1 = %.3e Pa & T1 = %.2f K'  % (P1,T1)  )
# print('CJ Speed   %.1f m/s' % cj_speed)
# print('CJ State')
# print('   Pressure   %.3e Pa' % gas.P)
# print('   Temperature  %.1f K' % gas.T)
# print('   Density  %.3f kg/m3' % gas.density)
# print('   Entropy  %.3e J/K' % gas.entropy_mass)
# print('   w2 (wave frame) %.1f m/s' % w2)
# print('   u2 (lab frame) %.1f m/s' % u2)
# print('   c2 frozen %.1f m/s' % af)
# print('   c2 equilbrium %.1f m/s' % ae)
# print('   gamma2 frozen %.3f ' % gammaf)
# print('   gamma2 equilbrium %.3f ' % gammae)


def cj_state(pressure, temperature, mixture, mech):
    cj_speed = CJspeed(pressure, temperature, mixture, mech)
    gas = PostShock_eq(cj_speed, pressure, temperature, mixture, mech)
    return cj_speed, gas


def pressure_CJ(pressures, temperature, mixture, mech):
    cj_speeds = []
    states = []
    for pressure in pressures:
        cj_speed, state = cj_state(pressure, temperature, mixture, mech)
        cj_speeds.append(cj_speed)
        states.append(state)
    return cj_speeds, states


def write_output(filename: str, x: list, ys: list[list], x_label: str, y_labels: list[str]):
    path = 'output/' + filename
    with open(path, "w") as output:
        labels = [x_label]
        labels.extend(y_labels)
        labels.append('\n')
        output.write(', '.join(labels))
        for i in range(len(x)):
            values = [x[i]]
            values.extend([y[i] for y in ys])
            values = [str(x) for x in values]
            output.write(', '.join(values))
            output.write('\n')


def get_top_weight_fraction(gas, n=14):
    d = [(species_name, gas.Y[k])for k, species_name in enumerate(gas.species_names)]
    d.sort(key=lambda species: species[1], reverse=True)
    return d[:n]


def get_top_molar_fraction(gas, n=14):
    d = [(species_name, gas.X[k])for k, species_name in enumerate(gas.species_names)]
    d.sort(key=lambda species: species[1], reverse=True)
    return d[:n]


def write_top_species(filename: str, molar=None, weight=None):
    path = 'output/' + filename
    with open(path, "w") as output:
        if molar:
            output.write('Top molar fractions:\n')
            for name, fraction in molar:
                output.write(f'{name:>12}  {fraction}\n')
        if weight:
            output.write('Top weight fractions:\n')
            for name, fraction in weight:
                output.write(f'{name:>12}  {fraction}\n')


pressures = [20000, 40000, 80000]
T1 = 293


mechs = {
    'CRECK': 'mechs/CRECK/CRECK_2003_TOT_HT_LT.yaml',
    'CRECK-sink': 'mechs/CRECK/CRECK_sink.yaml',
}

mixtures = {
    '10C2H2': 'C2H2:10 AR:90',
    '20C2H2': 'C2H2:20 AR:80',
    '30C2H2': 'C2H2:30 AR:70',
}


for mech_label, mech in mechs.items():
    for mixture_label, mixture in mixtures.items():
        cj_velocities, states = pressure_CJ(pressures, T1, mixture, mech)
        cj_pressures = [state.P for state in states]
        filename = 'test_{mixture}_{mech}.out'.format(
            mixture=mixture_label,
            mech=mech_label
        )
        write_output(filename, pressures,
                     [cj_velocities, cj_pressures],
                     'Pressure[Pa]',
                     ['CJ_velocity[m/s]', 'CJ_pressure[Pa]'])
        gas_for_analysis = states[-1]
        top_molar = get_top_molar_fraction(gas_for_analysis)
        top_weight = get_top_weight_fraction(gas_for_analysis)
        analysis_filename = filename.split('.')[0] + '.spc'
        write_top_species(analysis_filename, top_molar, top_weight)
