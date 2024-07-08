import numpy as np
import pandas as pd
import time
import cantera as ct
import matplotlib.pyplot as plt
from pathlib import Path
from scripts.utils.mixture import get_trifuel_for_o2
from scripts.utils.csv_processors import sum_columns

from configs.constants import OUTPUT_DIR, BATCH_OUTPUT
from configs.species_slices import BIN_X_SLICES, Y_SOOT_SLICES


def get_induction_time(times, concentrations):
    max_slope = 0
    time_of_max = 0
    value_at_max = 0
    for i in range(1, len(concentrations)):
        slope = (concentrations[i] - concentrations[i-1]) / (times[i] - times[i-1])
        if slope > max_slope:
            max_slope = slope
            time_of_max = times[i]
            value_at_max = concentrations[i]
    time_of_tangent = time_of_max - value_at_max/max_slope
    return time_of_max, time_of_tangent


def get_induction_time_from_solution_array(states: ct.SolutionArray, species: str):
    times = states.t
    concentrations = [x[0] for x in states(species).X]
    return get_induction_time(times, concentrations)


def get_solution(gas, T, P, mixture, max_time=1):
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
    return time_history


def write_csv(time_history, output_path='default', label='default', cols=None):
    time_history.write_csv(output_path / f'{label}_X.csv', cols=cols, species="X")
    time_history.write_csv(output_path / f'{label}_Y.csv', cols=cols, species="Y")


def get_ignition_delay(gas, T, P, mixture, max_time=1, reference_species='OH'):
    time_history = get_solution(gas, T, P, mixture, max_time)
    t_max_slope, t_tangent = get_induction_time_from_solution_array(time_history, reference_species)
    return t_max_slope, t_tangent


def get_IDT_temperature_dependence(gas, Ts, P, mixture, max_time=0.1, reference_species='OH'):
    taus_max_slope = []
    taus_tangent = []
    for temperature in Ts:
        t_max_slope, t_tangent = get_ignition_delay(
            gas=gas,
            T=temperature,
            P=P,
            mixture=mixture,
            max_time=max_time,
            reference_species=reference_species
        )
        taus_max_slope.append(t_max_slope)
        taus_tangent.append(t_tangent)
        print(temperature, t_max_slope, t_tangent)
    return taus_max_slope, taus_tangent


def get_manymodel_idt_tempertaute_dependence(
        mechs: dict, Ts: list, P: float, mixture: str,
        max_time=0.1, reference_species='OH', mode='tangent') -> pd.DataFrame:
    output = pd.DataFrame(Ts, columns=['T[K]'])
    for label, mech in mechs.items():
        print(f'Obtainig IDT dependence with mech "{label}"')
        gas = ct.Solution(mech, 'gas')
        taus_max_slope, taus_tangent = get_IDT_temperature_dependence(
            gas, Ts, P, mixture, max_time, reference_species
        )
        if mode == 'max_slope':
            output[label] = taus_max_slope
        else:
            output[label] = taus_tangent
    return output


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
    for r in range(gas.n_reactions):
        if 'H2O2' not in gas.reaction_equations()[r]:
            print('skipped reaction', r)
            continue
        print(f'reaction {r} of {gas.n_reactions}')
        gas.set_multiplier(1.0)  # reset all multipliers
        gas.set_multiplier(1 + dk, r)  # perturb reaction m
        t = get_ignition_delay(gas, T, P, mixture)
        print(f'Changing {gas.reaction_equations()[r]} t = {t*1e6:.2f} mks') 
        sensitivity = np.log(t/t0) / np.log(1+dk)
        sensitivities.append((gas.reaction_equations()[r], sensitivity))
        sensitivities.sort(key=lambda x: abs(x[1]), reverse=True)
    gas.set_multiplier(1.0)
    return sensitivities


MECHS = {
    'GRI': 'mechs/GRI/gri30_highT.yaml',
    'CRECK': 'mechs/CRECK/CRECK_2003_TPRF_HT_LT_ALC_ETHERS.yaml',
    'Aramco': 'mechs/Aramco/aramco2.yaml',
    'FFCM': 'mechs/FFCM/FFCM1.yaml',
    'BabuCRECK-NH3': 'mechs/modified/BabuCRECK-NH3.yaml',
    'BabuCRECK-NH3-ALL': 'mechs/modified/BabuCLBRCRECK.yaml',
    'Hong2011': 'mechs/Hong2011.yaml',
    'AceHalo': 'mechs/AcetyleneHalo/AceHalo.yaml'
}

NH3MECHS = {
    'Glarborg': 'mechs/NH3/glarborg.yaml',
    'KAUST': 'mechs/NH3/KAUST_NH3.yaml',
    'Konnov': 'mechs/NH3/konnov.yaml',
    'Okafor': 'mechs/NH3/okafor.yaml',
    'CRECK': 'mechs/NH3/CRECK_2003_C1_C3_HT_NOX.yaml'
}

temperatures = [1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250, 1275, 
                1300, 1333, 1366, 1400, 1433, 1466, 
                1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850,
                1900, 1950, 2000, 2100]
mech = 'AceHalo'
gas = ct.Solution(MECHS[mech], 'gas')

pressure = 1.2e6
mixtures = {'pure': 'C2H2:10 AR:90',
            'C2F4Br2': 'C2H2:10 C2BR2F4:1 AR:89',
            'CCl4': 'C2H2:10 CCL4:1 AR:89',
            'CF3I': 'C2H2:10 CF3I:1 AR:89',
            }

for mixture_label, mixture in mixtures.items():
    report_path = Path() / OUTPUT_DIR / BATCH_OUTPUT / mixture_label
    with report_path.open("w", encoding ="utf-8") as f:
        f.write('T[K],t_max_C4H2, t_max_C6H6, t_max_C16H10, t_ind\n')
    for temperature in temperatures:
        label = f'{mixture_label}_{temperature}K'
        output_path = Path() / OUTPUT_DIR / BATCH_OUTPUT / label
        output_path.mkdir(parents=True, exist_ok=True)
        solution = get_solution(gas, temperature, pressure, mixture, max_time=0.01)
        write_csv(solution, output_path, label)
        sum_columns(
            output_path / f'{label}_Y.csv',
            old_columns=['t', 'T', 'density'],
            new_columns=Y_SOOT_SLICES)
        cols = ('C2H2', 'C4H2', 'C6H6', 'C16H10')
        write_csv(solution(*cols), output_path, label+'short')
        result = pd.read_csv(output_path / f'{label}_Y_sums.csv')
        _, t_ind =  get_induction_time(result['t'], result['Y5'])
        Path.unlink(output_path / f'{label}_X.csv')
        Path.unlink(output_path / f'{label}_Y.csv')
        result = pd.read_csv(output_path / f'{label}short_X.csv')
        t_max_c4h2 = result['t'][result['X_C4H2'].argmax()]
        t_max_c6h6 = result['t'][result['X_C6H6'].argmax()]
        t_max_c16h10 = result['t'][result['X_C16H10'].argmax()]
        with report_path.open("a", encoding ="utf-8") as f:
            f.write(f'{temperature},{t_max_c4h2},{t_max_c6h6},{t_max_c16h10},{t_ind}\n')

# sensitivities = idt_sensitivity(gas, temperature, pressure, mixture)
# for x in sensitivities[:20]:
#     print(f'{x[0]}, {x[1]:.5f}')

# temperatures = np.linspace(1400, 2000, 25)
# output = get_manymodel_idt_tempertaute_dependence(NH3MECHS, temperatures, 4.3e5, 'NH3:9.33 CH4:6.66 O2:20.33 AR:163.67')
# output.to_csv('output/BatchReactor/model-compare-NH3-CH4.csv')

# temperatures = np.linspace(1000, 1960, 49)
# mech = 'BabuCRECK-NH3-ALL'
# gas = ct.Solution(MECHS[mech], 'gas')
# print(gas.n_species)
# pressure = 5.0*1e5
# mixture = 'NH3:7.67 H2:2.8 O2:7 CCL4:1 AR:81.53'

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