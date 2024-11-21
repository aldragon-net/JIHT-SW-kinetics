import numpy as np
import pandas as pd
import enum
import time
import cantera as ct
import matplotlib.pyplot as plt
from pathlib import Path
from scripts.utils.mixture import get_trifuel_for_o2
from scripts.utils.csv_processors import sum_columns

from configs.constants import OUTPUT_DIR, BATCH_OUTPUT
from configs.species_slices import BIN_X_SLICES, Y_SOOT_SLICES

class InductionTimeType(enum.Enum):
    TANGENT = 'tangent'
    SLOPE = 'slope'


class InductionTime:
    def __init__(self, slope: float, tangent: float) -> None:
        self.slope = slope
        self.tangent = tangent


def get_induction_time(times, concentrations) -> InductionTime:
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
    induction_time = InductionTime(slope=time_of_max, tangent=time_of_tangent)
    return induction_time


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
    induction_time = get_induction_time_from_solution_array(time_history, reference_species)
    return induction_time


def get_IDT_temperature_dependence(gas, Ts, P, mixture, max_time=0.1, reference_species='OH'):
    induction_times = []
    for temperature in Ts:
        induction_time = get_ignition_delay(
            gas=gas,
            T=temperature,
            P=P,
            mixture=mixture,
            max_time=max_time,
            reference_species=reference_species
        )
        induction_times.append(induction_time)
    return induction_times


def get_IDT_pressure_dependence(gas, T, Ps, mixture, max_time=0.1, reference_species='OH'):
    taus_max_slope = []
    taus_tangent = []
    for pressure in Ps:
        induction_time = get_ignition_delay(
            gas=gas,
            T=T,
            P=pressure,
            mixture=mixture,
            max_time=max_time,
            reference_species=reference_species
        )
        taus_max_slope.append(induction_time.slope)
        taus_tangent.append(induction_time.tangent)
    return taus_max_slope, taus_tangent


def get_IDT_mixture_dependence(gas, T, P, mixtures, max_time=0.1, reference_species='OH'):
    taus_max_slope = []
    taus_tangent = []
    for mixture in mixtures:
        t_max_slope, t_tangent = get_ignition_delay(
            gas=gas,
            T=T,
            P=P,
            mixture=mixture,
            max_time=max_time,
            reference_species=reference_species
        )
        taus_max_slope.append(t_max_slope)
        taus_tangent.append(t_tangent)
    return taus_max_slope, taus_tangent


def get_manymodel_idt_temperature_dependence(
        mechs: dict, Ts: list, P: float, mixture: str,
        max_time=0.1, reference_species='OH', mode=InductionTimeType.TANGENT.value) -> pd.DataFrame:
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


def get_manymodel_idt_pressure_dependence(
        mechs: dict, T: float, Ps: list, mixture: str,
        max_time=0.1, reference_species='OH', mode=InductionTimeType.TANGENT.value) -> pd.DataFrame:
    output = pd.DataFrame(Ps, columns=['P[Pa]'])
    for label, mech in mechs.items():
        print(f'Obtainig IDT dependence with mech "{label}"')
        gas = ct.Solution(mech, 'gas')
        taus_max_slope, taus_tangent = get_IDT_pressure_dependence(
            gas, T, Ps, mixture, max_time, reference_species
        )
        if mode == 'max_slope':
            output[label] = taus_max_slope
        else:
            output[label] = taus_tangent
    return output

def get_manymodel_idt_mixture_dependence(
        mechs: dict, T: float, P: float, fractions: list, mixtures: list,
        max_time=0.1, reference_species='OH', mode=InductionTimeType.TANGENT.value) -> pd.DataFrame:
    output = pd.DataFrame(fractions, columns=['fraction'])
    output['mixtures'] = mixtures
    for label, mech in mechs.items():
        print(f'Obtainig IDT dependence with mech "{label}"')
        gas = ct.Solution(mech, 'gas')
        taus_max_slope, taus_tangent = get_IDT_mixture_dependence(
            gas, T, P, mixtures, max_time, reference_species
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


def idt_sensitivity(gas, T, P, mixture, dk=0.05, idt_mode: str = InductionTimeType.TANGENT.value):
    sensitivities = []
    gas.set_multiplier(1.0)
    t0 = getattr(get_ignition_delay(gas, T, P, mixture), idt_mode)
    print(f'Undisturbed tau = {t0*1e6:.2f} mks')
    for r in range(gas.n_reactions):
        print(f'reaction {r} of {gas.n_reactions}')
        gas.set_multiplier(1.0)  # reset all multipliers
        gas.set_multiplier(1 + dk, r)  # perturb reaction m
        try:
            t = getattr(get_ignition_delay(gas, T, P, mixture), idt_mode)
            print(f'Changing {gas.reaction_equations()[r]} t = {t*1e6:.2f} mks')
            sensitivity = np.log(t/t0) / np.log(1+dk)
        except:
            print(f'Error while changing {gas.reaction_equations()[r]}')
            sensitivity = 666
        sensitivities.append(((r, gas.reaction_equations()[r]), sensitivity))
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
    'CRECK': 'mechs/NH3/CRECK_2003_C1_C3_HT_NOX.yaml',
    'Faravelli': 'mechs/NH3/NH3-Faravelli.yaml',
    'LiHeZhu': 'mechs/NH3/NH3-LiHeZhu.yaml'
}

temperatures = [1300, 1333, 1366, 1400, 1433, 1466, 1500,
                1550, 1600, 1650, 1700, 1750, 1800, 1850,
                1900, 1950, 2000]

def investigate_nh3():
    """Анализ зависимостей задержки воспламенения в системах NH3+CH4/C2H2/C2H4/C2H6"""
    AMMONIA = 'NH3'
    O2_FRACTION = 7
    admixtures = ['CH4', 'C2H2', 'C2H4', 'C2H6']
    temperatures = np.linspace(1100, 2000, 37)
    temperatures_short = [1200, 1400, 1700]
    pressures = [p*1e5 for p in np.linspace(2, 15, 14)]
    pressures_short = [3e5, 12e5]
    alphas = [0, 5, 10, 20, 30, 50, 70, 80, 90, 95, 100]
    alphas_short = [0, 10, 30, 50]
    for admixture in admixtures:
        for pressure in pressures_short:
            for temperature in temperatures_short:
                mixtures = []
                for alpha in alphas:
                    mixture = get_trifuel_for_o2(o2_fraction=O2_FRACTION, primary=AMMONIA, secondary=admixture, tertiary='', alpha=alpha, beta=0)
                    mixtures.append(mixture)
                mixoutput = get_manymodel_idt_mixture_dependence(mechs=NH3MECHS, T=temperature, P=pressure, fractions=alphas, mixtures=mixtures)
                mixoutput.to_csv(f'output/BatchReactor/NH3/NH3-{admixture}-T{temperature:.0f}K-P{pressure/1e5:.0f}bar-alphas.csv')
        for alpha in alphas_short:
            mixture = get_trifuel_for_o2(o2_fraction=O2_FRACTION, primary=AMMONIA, secondary=admixture, tertiary='', alpha=alpha, beta=0)
            for pressure in pressures_short:
                tempoutput = get_manymodel_idt_temperature_dependence(
                    mechs=NH3MECHS, Ts=temperatures, P=pressure, mixture=mixture
                )
                tempoutput.to_csv(f'output/BatchReactor/NH3/NH3-{admixture}-alpha{alpha:.0f}-P{pressure/1e5:.0f}bar-temperatures.csv')
            for temperature in temperatures_short:
                presoutput = get_manymodel_idt_pressure_dependence(
                    mechs=NH3MECHS, T=temperature, Ps=pressures, mixture=mixture
                )
                presoutput.to_csv(f'output/BatchReactor/NH3.NH3-{admixture}-alpha{alpha:.0f}-T{temperature:.0f}K-pressures.csv')  

def analyze_idt_sensitivity(gas, temperature, pressure, mixture, label='sensitivity', limit=None, dk=0.5):
    sensitivities = idt_sensitivity(gas, temperature, pressure, mixture, dk=dk)
    with open(f'output/IDTsens_{label}.out', 'w') as output:
        output.write(f'IDT sensitivity analysis for {mixture} at {temperature}K {pressure/1e5:.3f} bar:\n')
        output.write(f'dk = {dk}\n\n')
        for x in sensitivities[:limit]:
            output.write(f'R{x[0][0]}: {x[0][1]} \t\t {x[1]:.5f}\n')


def manymodel_idt_sensitivity(mechs: dict, temperature: float, pressure: float, mixture: str, mixlabel: str = 'mix'):
    for mechlabel, mech in mechs.items():
        gas = ct.Solution(mech, 'gas')
        analyze_idt_sensitivity(gas, temperature, pressure, mixture, label=f'{mixlabel}_{mechlabel}')


mixtures_for_analysis = {
    # 'NH3': 'NH3:9.333 O2:7.000 AR:83.667',
    # 'CH4': 'CH4:3.5 O2:7.000 AR:89.50',
    # 'C2H2': 'C2H2:2.8 O2:7.000 AR:90.20',
    # 'C2H4': 'C2H4:2.333 O2:7.000 AR:90.667',
    # 'C2H6': 'C2H6:2.0 O2:7.000 AR:91',
    # 'CH4_a10': 'NH3:8.400 CH4:0.350 O2:7.000 AR:84.250',
    'CH4_a30': 'NH3:6.533 CH4:1.050 O2:7.000 AR:85.417',
    'C2H2_a10': 'NH3:8.400 C2H2:0.280 O2:7.000 AR:84.320',
    'C2H2_a30': 'NH3:6.533 C2H2:0.840 O2:7.000 AR:85.627',
    'C2H4_a10': 'NH3:8.400 C2H4:0.233 O2:7.000 AR:84.367',
    'C2H4_a30': 'NH3:6.533 C2H4:0.700 O2:7.000 AR:85.767',
    'C2H6_a10': 'NH3:8.400 C2H6:0.200 O2:7.000 AR:84.400',
    'C2H6_a30': 'NH3:6.533 C2H6:0.600 O2:7.000 AR:85.867'
}

test_mechs = {
    'HongI': 'mechs/Hong2011.yaml',
    'HongII': 'mechs/Hong2011.yaml'
}

# output = get_manymodel_idt_temperature_dependence(NH3MECHS, temperatures, 9e5, 'NH3:9.333 O2:7.000 AR:83.67')
# output.to_csv('output/BatchReactor/model-compare-pureNH3-7bar.csv')

# gas = ct.Solution('mechs/NH3/konnov.yaml', 'gas')
# times = get_IDT_temperature_dependence(gas, temperatures, 700000, 'CH4:3.5 O2:7.0 AR:89.5', 2e-3)
# for i in range(len(temperatures)):
#     print(f'{temperatures[i]}, {times[1][i]:.4e}')


def multimixtures_manymodel_idt_sensitivity(mixtures: dict, mechs: dict, temperature: float, pressure: float):
    for mixlabel, mixture in mixtures.items():
        manymodel_idt_sensitivity(mechs=mechs, temperature=temperature, 
                                  pressure=pressure, mixture=mixture, mixlabel=mixlabel
                                  )


multimixtures_manymodel_idt_sensitivity(mixtures_for_analysis, NH3MECHS, 1450, 8e5)




#pressure = 1.2e6
# mixtures = {'pure': 'C2H2:10 AR:90',
#             'C2F4Br2': 'C2H2:10 C2BR2F4:1 AR:89',
#             'CCl4': 'C2H2:10 CCL4:1 AR:89',
#             'CF3I': 'C2H2:10 CF3I:1 AR:89',
#             }

# mech = 'AceHalo'
# # gas = ct.Solution(MECHS[mech], 'gas')

# for mixture_label, mixture in mixtures.items():
#     report_path = Path() / OUTPUT_DIR / BATCH_OUTPUT / mixture_label
#     with report_path.open("w", encoding ="utf-8") as f:
#         f.write('T[K],t_max_C4H2, t_max_C6H6, t_max_C16H10, t_ind\n')
#     for temperature in temperatures:
#         label = f'{mixture_label}_{temperature}K'
#         output_path = Path() / OUTPUT_DIR / BATCH_OUTPUT / label
#         output_path.mkdir(parents=True, exist_ok=True)
#         solution = get_solution(gas, temperature, pressure, mixture, max_time=0.01)
#         write_csv(solution, output_path, label)
#         sum_columns(
#             output_path / f'{label}_Y.csv',
#             old_columns=['t', 'T', 'density'],
#             new_columns=Y_SOOT_SLICES)
#         cols = ('C2H2', 'C4H2', 'C6H6', 'C16H10')
#         write_csv(solution(*cols), output_path, label+'short')
#         result = pd.read_csv(output_path / f'{label}_Y_sums.csv')
#         _, t_ind =  get_induction_time(result['t'], result['Y5'])
#         Path.unlink(output_path / f'{label}_X.csv')
#         Path.unlink(output_path / f'{label}_Y.csv')
#         result = pd.read_csv(output_path / f'{label}short_X.csv')
#         t_max_c4h2 = result['t'][result['X_C4H2'].argmax()]
#         t_max_c6h6 = result['t'][result['X_C6H6'].argmax()]
#         t_max_c16h10 = result['t'][result['X_C16H10'].argmax()]
#         with report_path.open("a", encoding ="utf-8") as f:
#             f.write(f'{temperature},{t_max_c4h2},{t_max_c6h6},{t_max_c16h10},{t_ind}\n')


# for x in sensitivities[:20]:
#     print(f'{x[0]}, {x[1]:.5f}')

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