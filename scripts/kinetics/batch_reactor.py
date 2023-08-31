import numpy as np

import time

import cantera as ct

import matplotlib.pyplot as plt

from scripts.utils.mixture import get_trifuel_for_o2


gas = ct.Solution("mechs/gri30_highT.yaml")


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

temperatures = np.linspace(1100, 1900, 3)

alpha = 20
beta = 0

beta = 0
labels = []
alphas = [0, 10, 20, 100]
dependencies = []
for alpha in alphas:
    mixture = get_trifuel_for_o2(7, 'CH4', 'H2', 'CH3OH', alpha, beta)
    mixture_label = f'M_{alpha}_{beta}'
    labels.append(mixture_label)
    print(f'Temperature dependency for mixture {mixture_label} ({mixture}) :')
    taus = get_temperature_dependence(gas, temperatures, 400000, mixture)
    dependencies.append(taus)

output_labels = ['T']
output_labels.extend(labels)
print(', '.join(output_labels))
for j in range(len(temperatures)):
    print(f'{temperatures[j]:.0f}', end=', ')
    for i in range(len(alphas)):
        print(f'{dependencies[i][j]*1e6:.0f}', end=',')
    print()
