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


def write_csv(time_history, output_path='default', label='default'):
    time_history.save(output_path / f'{label}_X.csv', basis="mole")
    time_history.save(output_path / f'{label}_Y.csv', basis="mass")


def get_ignition_delay(gas, T, P, mixture, max_time=1, reference_species='OH'):
    time_history = get_solution(gas, T, P, mixture, max_time)
    induction_time = get_induction_time_from_solution_array(time_history, reference_species)
    return induction_time


def get_IDT_temperature_dependence(gas, Ts, P, mixture, max_time=0.1, reference_species='OH'):
    taus_max_slope = []
    taus_tangent = []
    for temperature in Ts:
        induction_time = get_ignition_delay(
            gas=gas,
            T=temperature,
            P=P,
            mixture=mixture,
            max_time=max_time,
            reference_species=reference_species
        )
        taus_max_slope.append(induction_time.slope)
        taus_tangent.append(induction_time.tangent)
    return taus_max_slope, taus_tangent


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
        induction_time = get_ignition_delay(
            gas=gas,
            T=T,
            P=P,
            mixture=mixture,
            max_time=max_time,
            reference_species=reference_species
        )
        taus_max_slope.append(induction_time.slope)
        taus_tangent.append(induction_time.tangent)
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


def idt_sensitivity(gas, T, P, mixture, dk=0.05, idt_mode: str = InductionTimeType.TANGENT.value, selected_reactions=None):
    sensitivities = []
    gas.set_multiplier(1.0)
    t0 = getattr(get_ignition_delay(gas, T, P, mixture), idt_mode)
    print(f'Undisturbed tau = {t0*1e6:.2f} mks')
    for r in range(gas.n_reactions):
        if selected_reactions:
            if r not in selected_reactions:
                continue
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
NH3_ALC_MECHS = {
    'KAUST': 'mechs/NH3/KAUST_NH3.yaml',
    'CRECK': 'mechs/NH3/CRECK_2003_C1_C3_HT_NOX.yaml',
    'Konnov': 'mechs/NH3/konnov.yaml',
    'Faravelli': 'mechs/NH3/NH3-Faravelli.yaml',
    'LiHeZhu': 'mechs/NH3/NH3-LiHeZhu.yaml',
    'Shrestha2022': 'mechs/NH3-ethers/shrestha2022.yaml',
    'Zhang2023': 'mechs/NH3-ethers/zhang2023.yaml',
    'Dai': 'mechs/NH3-ethers/Dai2024.yaml'
}

NH3_DME_MECHS = {
    # 'Konnov': 'mechs/NH3/konnov.yaml',
    'LiHeZhu': 'mechs/NH3/NH3-LiHeZhu.yaml',
    # 'Shrestha2022': 'mechs/NH3-ethers/shrestha2022.yaml',
    'Zhang2023': 'mechs/NH3-ethers/zhang2023.yaml',
    # 'Dai': 'mechs/NH3-ethers/Dai2024.yaml'
}

NH3_DEE_MECHS = {
    'Shrestha2022': 'mechs/NH3-ethers/shrestha2022.yaml',
    'Dai': 'mechs/NH3-ethers/Dai2024.yaml'
}

NSK_MODELS = {
    # 'MODEL_1': 'mechs/NSK/NSKmod_1.yaml',
    # 'MODEL_2': 'mechs/NSK/NSKmod_2.yaml',
    # 'MODEL_3': 'mechs/NSK/NSKmod_3.yaml',
    # 'MODEL_4': 'mechs/NSK/NSKmod_4.yaml',
    # 'MODEL_5': 'mechs/NSK/NSKmod_5.yaml',
    # 'MODEL_6': 'mechs/NSK/NSKmod_6.yaml',
    # 'MODEL_7': 'mechs/NSK/NSKmod_7.yaml',
    # 'MODEL_B': 'mechs/NSK/NSKmod_mechB.yaml',
    'Shrestha2025': 'mechs/NSK/shrestha2025.yaml',
    'Shrestha2025mod': 'mechs/NSK/shrestha2025mod.yaml',
    'Shrestha2025v2': 'mechs/NSK/shrestha2025v2.yaml',
    'Shrestha2025v3': 'mechs/NSK/shrestha2025v3.yaml',
    'Shrestha2025v4': 'mechs/NSK/shrestha2025v4.yaml',
    # 'Shrestha2021': 'mechs/NSK/shrestha2021.yaml',
    # 'POLIMI': 'mechs/NSK/polimi.yaml',
    # 'NUIG_Yin': 'mechs/NSK/NUIGMech1.1-M.yaml'
}

temperatures = [1250, 1275, 1300, 1333, 1366, 1400, 1433, 1466, 1500,
                1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 2000 ]

def investigate_nh3():
    """Анализ зависимостей задержки воспламенения в системах NH3+CH4/C2H2/C2H4/C2H6"""
    AMMONIA = 'NH3'
    O2_FRACTION = 7
    admixtures = ['CH4', 'C2H2', 'C2H4', 'C2H6']
    temperatures = np.linspace(1100, 2000, 37)
    temperatures_short = []
    pressures = [p*1e5 for p in np.linspace(2, 15, 14)]
    pressures_short = [7e5]
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
                presoutput.to_csv(f'output/BatchReactor/NH3/NH3-{admixture}-alpha{alpha:.0f}-T{temperature:.0f}K-pressures.csv')  

def investigate_C2H6_NSK():
    """Анализ зависимостей задержки воспламенения в системах NH3+CH4/C2H2/C2H4/C2H6"""
    AMMONIA = 'NH3'
    O2_FRACTION = 7
    admixtures = ['C2H6']
    temperatures = np.linspace(1100, 2000, 37)
    temperatures_short = []
    pressures = [p*1e5 for p in np.linspace(2, 15, 14)]
    pressures_short = [7e5]
    alphas = []
    alphas_short = [0, 10, 30, 100]
    for admixture in admixtures:
        for pressure in pressures_short:
            for temperature in temperatures_short:
                mixtures = []
                for alpha in alphas:
                    mixture = get_trifuel_for_o2(o2_fraction=O2_FRACTION, primary=AMMONIA, secondary=admixture, tertiary='', alpha=alpha, beta=0)
                    mixtures.append(mixture)
                mixoutput = get_manymodel_idt_mixture_dependence(mechs=NSK_MODELS, T=temperature, P=pressure, fractions=alphas, mixtures=mixtures)
                mixoutput.to_csv(f'output/BatchReactor/NH3-NSK/NH3-{admixture}-T{temperature:.0f}K-P{pressure/1e5:.0f}bar-alphas.csv')
        for alpha in alphas_short:
            mixture = get_trifuel_for_o2(o2_fraction=O2_FRACTION, primary=AMMONIA, secondary=admixture, tertiary='', alpha=alpha, beta=0)
            for pressure in pressures_short:
                tempoutput = get_manymodel_idt_temperature_dependence(
                    mechs=NSK_MODELS, Ts=temperatures, P=pressure, mixture=mixture
                )
                tempoutput.to_csv(f'output/BatchReactor/NH3-NSK/NH3-{admixture}-alpha{alpha:.0f}-P{pressure/1e5:.0f}bar-temperatures.csv')
            for temperature in temperatures_short:
                presoutput = get_manymodel_idt_pressure_dependence(
                    mechs=NSK_MODELS, T=temperature, Ps=pressures, mixture=mixture
                )
                presoutput.to_csv(f'output/BatchReactor/NH3-NSK/NH3-{admixture}-alpha{alpha:.0f}-T{temperature:.0f}K-pressures.csv')


def analyze_idt_sensitivity(gas, temperature, pressure, mixture, label='sensitivity', limit=None, dk=0.5, selected_reactions=None):
    sensitivities = idt_sensitivity(gas, temperature, pressure, mixture, dk=dk, selected_reactions=selected_reactions)
    with open(f'output/IDTsens_{label}.out', 'w') as output:
        output.write(f'IDT sensitivity analysis for {mixture} at {temperature}K {pressure/1e5:.3f} bar:\n')
        output.write(f'dk = {dk}\n\n')
        for x in sensitivities[:limit]:
            output.write(f'R{x[0][0]}: {x[0][1]} \t\t {x[1]:.5f}\n')


def manymodel_idt_sensitivity(mechs: dict, temperature: float, pressure: float, mixture: str, mixlabel: str = 'mix'):
    for mechlabel, mech in mechs.items():
        gas = ct.Solution(mech, 'gas')
        analyze_idt_sensitivity(gas, temperature, pressure, mixture, label=f'{mixlabel}_{mechlabel}')


def investigate_reaction_sensitivity():
    alphas = [0, 5, 10, 20, 30, 50, 70, 80, 90, 95, 100]
    mixtures = []
    for alpha in alphas:
        mixture = get_trifuel_for_o2(o2_fraction=7, primary='NH3', secondary='C2H2', tertiary='', alpha=alpha, beta=0)
        mixtures.append(mixture)
    print(mixtures)
    gas = ct.Solution(NH3MECHS['KAUST'], 'gas')
    result = []
    for i, mixture in enumerate(mixtures):
        result.append((alphas[i], idt_sensitivity(gas, T=1350, P=8e5, dk=0.5, mixture=mixture, selected_reactions=[415])))
    for s in result:
        print(f'{s[0]}, {s[1][0][1]}')


# investigate_reaction_sensitivity()

def multimixtures_manymodel_idt_sensitivity(mixtures: dict, mechs: dict, temperature: float, pressure: float):
    for mixlabel, mixture in mixtures.items():
        manymodel_idt_sensitivity(mechs=mechs, temperature=temperature, 
                                  pressure=pressure, mixture=mixture, mixlabel=mixlabel
  )

# gas =  ct.Solution(NH3MECHS['LiHeZhu'], 'gas')
# output_path = Path() / OUTPUT_DIR / BATCH_OUTPUT
# solution = get_solution(gas, 1500, 7e5, 'NH3:9.333 O2:7.000 AR:83.667', max_time=0.002)
# write_csv(solution, output_path, 'AMM_1500K')
# solution = get_solution(gas, 1600, 7e5, 'NH3:9.333 O2:7.000 AR:83.667', max_time=0.002)
# write_csv(solution, output_path, 'AMM_1600K')
# solution = get_solution(gas, 1750, 7e5, 'NH3:9.333 O2:7.000 AR:83.667', max_time=0.002)
# write_csv(solution, output_path, 'AMM_1750K')
# solution = get_solution(gas, 1500, 7e5, 'CH4:3.5 O2:7.000 AR:89.50', max_time=0.002)
# write_csv(solution, output_path, 'MET100_1500K')
# solution = get_solution(gas, 1600, 7e5, 'CH4:3.5 O2:7.000 AR:89.50', max_time=0.002)
# write_csv(solution, output_path, 'MET100_1600K')
# solution = get_solution(gas, 1750, 7e5, 'CH4:3.5 O2:7.000 AR:89.50', max_time=0.002)
# write_csv(solution, output_path, 'MET100_1750K')

mixtures_for_analysis = {}

mixture_ammonia = {'NH3': 'NH3:9.333 O2:7.000 AR:83.667'}
mixtures_with_CH4 = {
    'CH4': 'CH4:3.5 O2:7.000 AR:89.50',
    'CH4_a10': 'CH4:0.350 NH3:8.400 O2:7.000 AR:84.250',
    'CH4_a30': 'CH4:1.050 NH3:6.533 O2:7.000 AR:85.417'
}
mixtures_with_C2H2 = {
    'C2H2': 'C2H2:2.8 O2:7.000 AR:90.20',
    'C2H2_a10': 'C2H2:0.280 NH3:8.400 O2:7.000 AR:84.320',
    'C2H2_a30': 'C2H2:0.840 NH3:6.533 O2:7.000 AR:85.627'
}
mixtures_with_C2H4 = {
    'C2H4': 'C2H4:2.333 O2:7.000 AR:90.667',
    'C2H4_a10': 'C2H4:0.233 NH3:8.400 O2:7.000 AR:84.367',
    'C2H4_a30': 'C2H4:0.700 NH3:6.533 O2:7.000 AR:85.767'
}
mixtures_with_C2H6 = {
    'C2H6': 'C2H6:2.0 O2:7.000 AR:91',
    'C2H6_a10': 'C2H6:0.200 NH3:8.400 O2:7.000 AR:84.400',
    'C2H6_a30': 'C2H6:0.600 NH3:6.533 O2:7.000 AR:85.867'
}
mixtures_with_CH3OH = {
    'CH3OH': 'CH3OH:4.667 O2:7.000 AR:88.333',
    'CH3OH_a10': 'CH3OH:0.467 NH3:8.400 O2:7.000 AR:84.133',
    'CH3OH_a30': 'CH3OH:1.400 NH3:6.533 O2:7.000 AR:85.067'
}
mixtures_with_C2H5OH = {
    'C2H5OH': 'C2H5OH:2.333 O2:7.000 AR:90.667',
    'C2H5OH_a10': 'C2H5OH:0.233 NH3:8.400 O2:7.000 AR:84.367',
    'C2H5OH_a30': 'C2H5OH:0.7 NH3:6.533 O2:7.000 AR:85.767'
}
mixtures_with_DME = {
    'DME': 'CH3OCH3:2.333 O2:7.000 AR:90.667',
    'DME_a10': 'CH3OCH3:0.233 NH3:8.400 O2:7.000 AR:84.367',
    'DME_a30': 'CH3OCH3:0.7 NH3:6.533 O2:7.000 AR:85.767'
}
mixtures_with_DEE = {
    'DEE': 'DEE:1.167 O2:7.000 AR:91.833',
    'DEE_a10': 'DEE:0.117 NH3:8.400 O2:7.000 AR:84.483',
    'DEE_a30': 'NH3:6.533 DEE:0.35 O2:7.000 AR:86.117'
}
mixtures_with_FURAN = {
    'FURAN': 'FURAN:1.55 O2:7.00 AR:91.45',
    'FURAN_a10': 'FURAN:0.155 NH3:8.400 O2:7 AR:84.445',
    'FURAN_a30': 'FURAN:0.467 NH3:6.533 O2:7 AR:86.000'
}
mixtures_with_THFURAN = {
    'THFURAN': 'THFURAN:1.272 O2:7.00 AR:91.728',
    'THFURAN_a10': 'THFURAN:0.127 NH3:8.400 O2:7 AR:84.473',
    'THFURAN_a30': 'THFURAN:0.382 NH3:6.533 O2:7 AR:86.085'
}
mixtures_with_TOLUENE = {
    'TOLUENE': 'C6H5CH3:0.778 O2:7.000 AR:92.222',
    'TOLUENE_a10': 'C6H5CH3:0.0778 NH3:8.400 O2:7.000 AR:84.5222',
    'TOLUENE_a30': 'C6H5CH3:0.233 NH3:6.533 O2:7.000 AR:86.234'
}
mixtures_with_CYC6H12 = {
    'CYC6H12': 'CYC6H12:0.778 O2:7.000 AR:92.222',
    'CYC6H12_a10': 'CYC6H12:0.0778 NH3:8.400 O2:7.000 AR:84.5222',
    'CYC6H12_a30': 'CYC6H12:0.233 NH3:6.533 O2:7.000 AR:86.234'
}
mixtures_with_NHEPTHANE = {
    'NHEPTHANE': 'NC7H16:0.636 O2:7.000 AR:92.364',
    'NHEPTHANE_a10': 'NC7H16:0.0636 NH3:8.400 O2:7.000 AR:84.5364',
    'NHEPTHANE_a30': 'NC7H16:0.191 NH3:6.533 O2:7.000 AR:86.276'
}

FURAN_MECHS = {
    'Sirjean13': 'mechs/furans/Sirjean13.yaml',
    'Tran17': 'mechs/furans/Tran17.yaml',
    'Fenard19': 'mechs/furans/Fenard19.yaml',
    'Wu20': 'mechs/furans/Wu20.yaml',
    'Somers': 'mechs/furans/Somers.yaml',
    'Cheng23': 'mechs/furans/Cheng23.yaml',
}

FURAN_NH3_MECHS = {
    'Tran17_KAUST_NH3': 'mechs/furans/Tran17-KAUST_NH3.yaml',
    'Wu20_KAUST_NH3': 'mechs/furans/Wu20-KAUST_NH3.yaml'
}

C3Mech4_MECHS = {
    'C3Mech_3XQ9': 'mechs/C3Mech401/C3MechV4.0.1_3XQ9_C0-C7_C5CY_C6CY_N_PAH_HT.yaml'
}


# multimixtures_manymodel_idt_sensitivity(mixtures_for_analysis, {'Shrestha2025': 'mechs/NSK/shrestha2025.yaml'}, 1400, 8.0e5)
    

# output_path = Path() / OUTPUT_DIR / BATCH_OUTPUT
# solution = get_solution(gas, 1400, 3e5, 'THFURAN:1.272e-3 O2:7.00e-3 AR:91.728', max_time=0.002)
# write_csv(solution, output_path, 'modCheng')



# FURANS 


temperatures = [1333, 1366, 1400, 1433, 1466, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850]

output = get_manymodel_idt_temperature_dependence(
    FURAN_NH3_MECHS, temperatures, 6.4e5, 'NH3:9.333 O2:7.000 AR:83.67')
output.to_csv('output/BatchReactor/NH3-furans/furamm_pureNH3-6.4bara.csv')

temperatures = [1166, 1200, 1233, 1266, 1300, 1333, 1366,
                1400, 1433, 1466, 1500, 1550, 1600, 1650]

output = get_manymodel_idt_temperature_dependence(
    FURAN_NH3_MECHS, temperatures, 7.8e5, 'FURAN:1.55 O2:7.000 AR:91.45')
output.to_csv('output/BatchReactor/NH3-furans/furamm-pureFURAN-7.8bar.csv')
output = get_manymodel_idt_temperature_dependence(
    FURAN_NH3_MECHS, temperatures, 7.5e5, 'THFURAN:1.272 O2:7.000 AR:91.7287')
output.to_csv('output/BatchReactor/NH3-furans/furamm-pureTHFURAN-7.5bar.csv')

output = get_manymodel_idt_temperature_dependence(
    FURAN_NH3_MECHS, temperatures, 6.3e5, 'FURAN:0.155 NH3:8.400 O2:7 AR:84.445')
output.to_csv('output/BatchReactor/NH3-furans/furamm-FURAN_a10-6.3bar.csv')
output = get_manymodel_idt_temperature_dependence(
    FURAN_NH3_MECHS, temperatures, 6.7e5, 'FURAN:0.467 NH3:6.533 O2:7 AR:86.000')
output.to_csv('output/BatchReactor/NH3-furans/furamm-FURAN_a30-6.7bar.csv')

output = get_manymodel_idt_temperature_dependence(
    FURAN_NH3_MECHS, temperatures, 6.2e5, 'THFURAN:0.127 NH3:8.400 O2:7 AR:84.473')
output.to_csv('output/BatchReactor/NH3-furans/furamm-THFURAN_a10-6.2bar.csv')
output = get_manymodel_idt_temperature_dependence(
    FURAN_NH3_MECHS, temperatures, 6.8e5, 'THFURAN:0.382 NH3:6.533 O2:7 AR:86.085')
output.to_csv('output/BatchReactor/NH3-furans/furamm-THFURAN_a30-6.8bar.csv')



# # HEAVY
# # ammonia
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 6.4e5, 'NH3:9.333 O2:7.000 AR:83.67')
# output.to_csv('output/BatchReactor/NH3-heavy/pureNH3-6.4bara.csv')

# # toluene
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 7.2e5, 'C6H5CH3:0.778 O2:7.000 AR:92.222')
# output.to_csv('output/BatchReactor/NH3-heavy/pureTOLUENE-7.2bara.csv')
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 6.6e5, 'C6H5CH3:0.0778 NH3:8.400 O2:7.000 AR:84.5222')
# output.to_csv('output/BatchReactor/NH3-heavy/TOLUENE_a10-6.6bara.csv')
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 6.9e5, 'C6H5CH3:0.233 NH3:6.533 O2:7.000 AR:86.234')
# output.to_csv('output/BatchReactor/NH3-heavy/TOLUENE_a30-6.9bara.csv')

# # cyhexane
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 7.6e5, 'CYC6H12:0.778 O2:7.000 AR:92.222')
# output.to_csv('output/BatchReactor/NH3-heavy/pureCYHEXANE-7.6bara.csv')
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 6.7e5, 'CYC6H12:0.0778 NH3:8.400 O2:7.000 AR:84.5222')
# output.to_csv('output/BatchReactor/NH3-heavy/CYHEXANE_a10-6.7bara.csv')
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 7.0e5, 'CYC6H12:0.233 NH3:6.533 O2:7.000 AR:86.234')
# output.to_csv('output/BatchReactor/NH3-heavy/CYHEXANE_a30-7.0bara.csv')

# # n-hepthane
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 7.9e5, 'NC7H16:0.636 O2:7.000 AR:92.364')
# output.to_csv('output/BatchReactor/NH3-heavy/pureNHEPTHANE-7.9bara.csv')
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 7.0e5, 'NC7H16:0.0636 NH3:8.400 O2:7.000 AR:84.5364')
# output.to_csv('output/BatchReactor/NH3-heavy/NHEPTHANE_a10-7.0bara.csv')
# output = get_manymodel_idt_temperature_dependence(
#     C3Mech4_MECHS, temperatures, 7.2e5, 'NC7H16:0.191 NH3:6.533 O2:7.000 AR:86.276')
# output.to_csv('output/BatchReactor/NH3-heavy/NHEPTHANE_a30-7.2bara.csv')


# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 12.7e5, 'NH3:9.333 O2:7.000 AR:83.67')
# output.to_csv('output/BatchReactor/NH3-alc-eth/pureNH3-12.7bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 4.1e5, 'NH3:9.333 O2:7.000 AR:83.67')
# output.to_csv('output/BatchReactor/NH3-alc-eth/pureNH3-4.1bar.csv')

# pure fuels
temperatures = [1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250,
                1275, 1300, 1333, 1366, 1400, 1433, 1466, 1500 ]
# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 8.8e5, 'CH3OH:4.667 O2:7.000 AR:88.333')
# output.to_csv('output/BatchReactor/NH3-alc-eth/pure-CH3OH-8.8bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 8.1e5, 'C2H5OH:2.333 O2:7.000 AR:90.667')
# output.to_csv('output/BatchReactor/NH3-alc-eth/pure-C2H5OH-8.1bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_DME_MECHS, temperatures, 8.2e5, 'CH3OCH3:2.333 O2:7.000 AR:90.667')
# output.to_csv('output/BatchReactor/NH3-alc-eth/pure-DME-8.2bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_DEE_MECHS, temperatures, 8.3e5, 'DEE:1.167 O2:7.000 AR:91.833')
# output.to_csv('output/BatchReactor/NH3-alc-eth/pure-DEE-8.3bar.csv')

# 10% 
# temperatures = [1275, 1300, 1325, 1350, 1375, 1400, 1425, 1450,
#                 1475, 1500, 1525, 1550, 1575, 1600, 1633, 1666, 1700 ]
# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 6.9e5, 'NH3:8.400 CH3OH:0.467 O2:7.000 AR:84.133')
# output.to_csv('output/BatchReactor/NH3-alc-eth/CH3OH_a10_6.9bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 6.7e5, 'NH3:8.400 C2H5OH:0.233 O2:7.000 AR:84.367')
# output.to_csv('output/BatchReactor/NH3-alc-eth/C2H5OH_a10_6.7bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_DME_MECHS, temperatures, 7.0e5, 'NH3:8.400 CH3OCH3:0.233 O2:7.000 AR:84.367')
# output.to_csv('output/BatchReactor/NH3-alc-eth/DME_a10_7.0bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_DEE_MECHS, temperatures, 6.9e5, 'NH3:8.400 DEE:0.117 O2:7.000 AR:84.483')
# output.to_csv('output/BatchReactor/NH3-alc-eth/DEE_a10_6.9bar.csv')

# 30% 
# temperatures = [1175, 1200, 1225, 1250, 1275, 1300, 1325, 1350,
#                 1375, 1400, 1425, 1450, 1475, 1500, 1533, 1566, 1600, 1650 ]
# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 7.6e5, 'NH3:6.533 CH3OH:1.400 O2:7.000 AR:85.067')
# output.to_csv('output/BatchReactor/NH3-alc-eth/CH3OH_a30_7.6bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_ALC_MECHS, temperatures, 7.2e5, 'NH3:6.533 C2H5OH:0.7 O2:7.000 AR:85.767')
# output.to_csv('output/BatchReactor/NH3-alc-eth/C2H5OH_a30_7.2bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_DME_MECHS, temperatures, 7.5e5, 'NH3:6.533 CH3OCH3:0.7 O2:7.000 AR:85.767')
# output.to_csv('output/BatchReactor/NH3-alc-eth/DME_a30_7.5bar.csv')

# output = get_manymodel_idt_temperature_dependence(
#     NH3_DEE_MECHS, temperatures, 7.4e5, 'NH3:6.533 DEE:0.35 O2:7.000 AR:86.117')
# output.to_csv('output/BatchReactor/NH3-alc-eth/DEE_a30_7.4bar.csv')



# gas = ct.Solution('mechs/NH3/konnov.yaml', 'gas')
# times = get_IDT_temperature_dependence(gas, temperatures, 700000, 'CH4:3.5 O2:7.0 AR:89.5', 2e-3)
# for i in range(len(temperatures)):
#     print(f'{temperatures[i]}, {times[1][i]:.4e}')



# investigate_C2H6_NSK()

# temperatures = [1250, 1275, 1300, 1333, 1366, 1400, 1433, 1466, 1500,
#                 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 2000 ]

# output = get_manymodel_idt_temperature_dependence(NSK_MODELS, temperatures, 7.7e5, 'NH3:8.400 C2H6:0.200 O2:7.000 AR:84.400')
# output.to_csv('output/BatchReactor/Shrestha2025Vesrions-10C2H6-7.7bar.csv')
# output = get_manymodel_idt_temperature_dependence(NSK_MODELS, temperatures, 8.0e5, 'NH3:6.533 C2H6:0.600 O2:7.000 AR:85.867')
# output.to_csv('output/BatchReactor/Shrestha2025Versions-30C2H6-8.0bar.csv')
# output = get_manymodel_idt_temperature_dependence(NSK_MODELS, temperatures, 7.0e5, 'NH3:9.333 O2:7.000 AR:83.667')
# output.to_csv('output/BatchReactor/Shrestha2025Versions-100NH3-7.0bar.csv')

# temperatures = [1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300, 1333, 1366, 1400, 1450, 1500]
# output = get_manymodel_idt_temperature_dependence(NSK_MODELS, temperatures, 9.0e5, 'C2H6:2.0 O2:7.000 AR:91')
# output.to_csv('output/BatchReactor/Shrestha2025Versions-100C2H6-9.0bar.csv')

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