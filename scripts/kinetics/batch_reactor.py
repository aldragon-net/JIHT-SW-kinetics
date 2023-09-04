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
    gas.TP = T, P
    gas.X = mixture
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
        print(f'T = {temperature:.0f}, IDT = {tau*1e6:.0f} mks')
    return taus


def get_mixture_label(alpha, beta, tertiary=''):
    mixture_label = f'M_{alpha}_{beta}'
    if tertiary and beta>0:
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
            taus = get_temperature_dependence(gas, temperatures, 400000, mixture)
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


MECHS = {
    'GRI': 'mechs/GRI/gri30.yaml',
    'CRECK': 'mechs/CRECK/CRECK_2003_TPRF_HT_LT_ALC_ETHERS.yaml',
    'Aramco': 'mechs/Aramco/aramco2.yaml',
    'FFCM': 'mechs/FFCM/FFCM1.yaml'
}

temperatures = np.linspace(1000, 1900, 46)

betas = [0, 30, 100]
alphas = [0,]
primary = 'CH4'
secondary = 'H2'
tertiary = 'CH3OCH3'
mech = 'CRECK'

gas = ct.Solution(MECHS[mech])
print(gas)

dependencies = get_dependence_on_alfa(gas, 7, primary, secondary, tertiary, alphas, betas)

output_dependence_on_alfa(dependencies, alphas, betas, temperatures, prefix=f'{mech}_{tertiary}_')
if len(betas) > 1:
    output_dependence_on_beta(dependencies, alphas, betas, temperatures, alpha_i=0, prefix=f'{mech}_{tertiary}_')