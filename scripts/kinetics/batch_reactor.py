import numpy as np
import time
import cantera as ct
import matplotlib.pyplot as plt
from scripts.utils.mixture import get_trifuel_for_o2


def get_time_of_max_slope(states, species):
    times = states.t
    concentrations = states(species).Y
    max_slope = 0
    time_of_max = 0
    for i in range(1, len(concentrations)):
        if i == 0:
            pass
        slope = (concentrations[i] - concentrations[i-1]) / (times[i] - times[i-1])
        if slope > max_slope:
            max_slope = slope
            time_of_max = times[i]
    return time_of_max


reference_species = "OH"


def get_ignition_delay(gas, T, P, mixture, max_time=1):
    gas.X = mixture
    gas.TP = T, P
    r = ct.IdealGasReactor(contents=gas, name="Batch Reactor")
    reactor_network = ct.ReactorNet([r])
    time_history = ct.SolutionArray(gas, extra="t")
    t = 0
    counter = 1
    while t < max_time:
        t = reactor_network.step()
        time_history.append(r.thermo.state, t=t)
        counter += 1
    tau = get_time_of_max_slope(time_history, reference_species)
    return tau


def get_temperature_dependence(gas, Ts, P, mixture):
    taus = []
    for temperature in Ts:
        tau = get_ignition_delay(gas, temperature, P, mixture)
        taus.append(tau)
    return taus


def get_mixture_label(alpha, beta, tertiary=''):
    mixture_label = f'M_{alpha}_{beta}'
    if tertiary and beta > 0:
        mixture_label = mixture_label + f'({tertiary})'
    return mixture_label


def get_dependence_on_alfa(gas, oxygen_ratio, primary, secondary, tertiary, alphas, betas):
    dependencies = []
    for alpha in alphas:
        on_betas = []
        for beta in betas:
            mixture = get_trifuel_for_o2(oxygen_ratio, primary, secondary, tertiary, alpha, beta)
            mixture_label = get_mixture_label(alpha, beta, tertiary)
            print(f'Temperature dependency for mixture {mixture_label} ({mixture}) :')
            taus = get_temperature_dependence(gas, temperatures, pressure, mixture)
            on_betas.append(taus)
        dependencies.append(on_betas)
    return dependencies


def output_dependence_on_alfa(dependencies, alphas, betas, temperatures, beta_i=0, path='output/', prefix=''):
    beta = betas[beta_i]
    # direct export (to plot dependence on temperature)
    mixture_labels = []
    for alpha in alphas:
        mixture_labels.append(get_mixture_label(alpha, beta))
    file_labels = ['T[K]', '10000/T[K-1]']
    file_labels.extend(mixture_labels)
    file_labels.extend([label+'_log' for label in mixture_labels])
    filename = path + f'{prefix}_M_x_{beta}.out'
    with open(filename, 'w') as f:
        f.write(', '.join(file_labels))
        for t, temperature in enumerate(temperatures):
            f.write('\n')
            line = [f'{temperature:.0f}', f'{1e4/temperature:.4f}']
            for a in range(len(alphas)):
                delay = dependencies[a][beta_i][t]
                line.append(f'{delay*1e6:.0f}')
            for a in range(len(alphas)):
                log_delay = np.log(1e6*dependencies[a][beta_i][t])
                line.append(f'{log_delay:.4f}')
            f.write(', '.join(line))
    # transposed export (to plot dependence on alpha)
    file_labels = ['alpha[prcnt]']
    file_labels.extend([f'T[{T:.0f}K]' for T in temperatures])   
    filename = path + f'{prefix}_M_x_{beta}_tran.out'
    with open(filename, 'w') as f:
        f.write(', '.join(file_labels))
        for a, alpha in enumerate(alphas):
            f.write('\n')
            line = [f'{alpha:.0f}']
            for t in range(len(temperatures)):
                delay = dependencies[a][beta_i][t]
                line.append(f'{delay*1e6:.0f}')
            f.write(', '.join(line))


def output_dependence_on_beta(dependencies, alphas, betas, temperatures, alpha_i=0, path='output/', prefix=''):
    alpha = alphas[alpha_i]
    # direct export (to plot dependence on temperature)
    mixture_labels = []
    for beta in betas:
        mixture_labels.append(get_mixture_label(alpha, beta))
    file_labels = ['T[K]', '10000/T[K-1]']
    file_labels.extend(mixture_labels)
    file_labels.extend([label+'_log' for label in mixture_labels])
    filename = path + f'{prefix}M_{alpha}_x.out'
    with open(filename, 'w') as f:
        f.write(', '.join(file_labels))
        for t, temperature in enumerate(temperatures):
            f.write('\n')
            line = [f'{temperature:.0f}', f'{1e4/temperature:.4f}']
            for b in range(len(betas)):
                delay = dependencies[alpha_i][b][t]
                line.append(f'{delay*1e6:.0f}')
            for b in range(len(betas)):
                log_delay = np.log(1e6*dependencies[alpha_i][b][t])
                line.append(f'{log_delay:.4f}')
            f.write(', '.join(line))
    # transposed export (to plot dependence on beta)
    file_labels = ['beta[prcnt]']
    file_labels.extend([f'T[{T:.0f}K]' for T in temperatures])   
    filename = path + f'{prefix}M_{alpha}_x_tran.out'
    with open(filename, 'w') as f:
        f.write(', '.join(file_labels))
        for b, beta in enumerate(betas):
            f.write('\n')
            line = [f'{beta:.0f}']
            for t in range(len(temperatures)):
                delay = dependencies[alpha_i][b][t]
                line.append(f'{delay*1e6:.0f}')
            f.write(', '.join(line))


def idt_sensitivity(gas, T, P, mixture, dk=0.05):
    sensitivities = []
    gas.set_multiplier(1.0)
    t0 = get_ignition_delay(gas, T, P, mixture)
    print(f'Undisturbed tau = {t0*1e6:.2f} mks')
    for r in range(503, 504): # gas.n_reactions):
        gas.set_multiplier(1.0)  # reset all multipliers
        gas.set_multiplier(1 + dk, r)  # perturb reaction m
        t = get_ignition_delay(gas, T, P, mixture)
        print(f'Changing {gas.reaction_equations()[r]} t = {t*1e6:.2f} mks') 
        sensitivity = np.log(t/t0) / np.log(1+dk)
        sensitivities.append((gas.reaction_equations()[r], sensitivity))
        sensitivities.sort(key = lambda x: abs(x[1]), reverse=True)
    gas.set_multiplier(1.0)
    return sensitivities


MECHS = {
    'GRI': 'mechs/GRI/gri30_highT.yaml',
    'CRECK': 'mechs/CRECK/CRECK_2003_TPRF_HT_LT_ALC_ETHERS.yaml',
    'Aramco': 'mechs/Aramco/aramco2.yaml',
    'FFCM': 'mechs/FFCM/FFCM1.yaml',
    'BabuCRECK-NH3': 'mechs/modified/BabuCRECK-NH3.yaml',
    'Hong2011': 'mechs/Hong2011.yaml' 
}

temperature = 1700
mech = 'BabuCRECK-NH3'
gas = ct.Solution(MECHS[mech], 'gas')
pressure = 4.1*1e5
mixture = 'NH3:9.33 O2:7.000 CF3I:1 AR:82.67'


sensitivities = idt_sensitivity(gas, temperature, pressure, mixture)
for x in sensitivities[:20]:
    print(f'{x[0]}, {x[1]:.5f}')

# temperatures = np.linspace(1000, 1960, 49)
# mech = 'BabuCRECK-NH3'
# gas = ct.Solution(MECHS[mech], 'gas')
# pressure = 5.1*1e5
# mixture = 'H2:14 O2:7 AR:79'

# dependence = get_temperature_dependence(gas, temperatures, pressure, mixture)

# !! alpha-beta-dependence !
# mech = 'GRI'
# gas = ct.Solution(MECHS[mech], 'gri30')
# alphas = [20]
# betas = [0, 20, 50]
# temperatures = np.linspace(1000, 1960, 49)
# pressure = 4.9e5
# primary = 'CH4'
# secondary = 'H2'
# tertiary = 'CH3OH'

# gas = ct.Solution(MECHS[mech])
# print(gas)

# dependencies = get_dependence_on_alfa(gas, 7, primary, secondary, tertiary, alphas, betas)

# output_dependence_on_alfa(dependencies, alphas, betas, temperatures, prefix=f'{mech}_{tertiary}_')
# if len(betas) > 1:
#     for i in range(len(alphas)):
#         output_dependence_on_beta(dependencies, alphas, betas, temperatures, alpha_i=i, prefix=f'{mech}_{tertiary}_')