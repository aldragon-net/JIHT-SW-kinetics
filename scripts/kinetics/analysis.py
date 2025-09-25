import numpy as np
import pandas as pd
import enum
import time
import cantera as ct
import matplotlib.pyplot as plt
from pathlib import Path


def get_reaction_rate_constants(gas: ct.Kinetics, reaction_number, T, P):
    gas.TP = T, P
    if reaction_number < 1 or reaction_number > gas.n_reactions:
        raise ValueError('Reaction number out of range')
    reaction_number -= 1
    forward_rate_constant = float(gas.forward_rate_constants[reaction_number])
    reverse_rate_constant = float(gas.reverse_rate_constants[reaction_number])
    return forward_rate_constant, reverse_rate_constant

def get_reaction_equilibrium_constant(gas: ct.Kinetics, reaction_number, T, P):
    gas.TP = T, P
    if reaction_number < 1 or reaction_number > gas.n_reactions:
        raise ValueError('Reaction number out of range')
    reaction_number -= 1
    equilibrium_constant = float(gas.equilibrium_constants[reaction_number])
    return equilibrium_constant


def get_reaction_rate_constant_temperature_dependencies(gas, reaction_number, temperatures, P):
    result = []
    for T in temperatures:
        forward_rate_constant, reverse_rate_constant = get_reaction_rate_constants(
            gas=gas, reaction_number=reaction_number, T=T, P=P
        )
        result.append((T, forward_rate_constant, reverse_rate_constant))
    return result


def get_reaction_equilibrium_constant_temperature_dependence(gas, reaction_number, temperatures, P):
    result = []
    for T in temperatures:
        equilibrium_constant = get_reaction_equilibrium_constant(
            gas=gas, reaction_number=reaction_number, T=T, P=P
        )
        result.append((T, equilibrium_constant))
    return result
                      


mech = 'mechs/NH3-ethers/zhang2023.yaml'
gas = ct.Solution(mech, 'gas')

def output_reaction_rate_constant_temperature_dependence(temperature_dependence, label='', info=None):
    with open(f'output/rates_{label}.txt', 'w') as f:
        if info:
            f.write(f'{info}\n\n')
        f.write('temperature[K], forward[cm3|mol|s], reverse[cm3|mol|s]\n')
        for i in range(len(temperature_dependence)):
            item = temperature_dependence[i]
            f.write(f'{item[0]},{item[1]*1e3:.5e},{item[2]*1e3:.5e}\n')
            # rate constants converted to [cm3/mol/s] from Cantera's [m3/kmol/s]

def output_reaction_equilibrium_constant_temperature_dependence(temperature_dependence, label='', info=None):
    with open(f'output/K_eq_{label}.txt', 'w') as f:
        if info:
            f.write(f'{info}\n\n')
        f.write('temperature[K], K_eq, 1/K_eq\n')
        for i in range(len(temperature_dependence)):
            item = temperature_dependence[i]
            f.write(f'{item[0]},{item[1]:.5e},{item[1]**-1:.5e}\n')

NH3_DEE_MECHS = {
    'Shrestha2022': 'mechs/NH3-ethers/shrestha2022.yaml',
    'Dai2024': 'mechs/NH3-ethers/Dai2024.yaml'
}

NH3MECHS = {
    'Glarborg': 'mechs/NH3/glarborg.yaml',
    'KAUST': 'mechs/NH3/KAUST_NH3.yaml',
    'Konnov': 'mechs/NH3/konnov.yaml',
    'Okafor': 'mechs/NH3/okafor.yaml',
    'CRECK': 'mechs/NH3/CRECK_2003_C1_C3_HT_NOX.yaml',
    'Faravelli': 'mechs/NH3/NH3-Faravelli.yaml',
    'LiHeZhu': 'mechs/NH3/NH3-LiHeZhu.yaml',
    'Shrestha2022': 'mechs/NH3-ethers/shrestha2022.yaml',
    'Zhang2023': 'mechs/NH3-ethers/zhang2023.yaml',
    'Dai2024': 'mechs/NH3-ethers/Dai2024.yaml'
}

temperatures = np.arange(300, 3001, 100)

# reactions_for_analysis = {
#     'NH2_HO2_in_Glarborg': ('Glarborg', NH3MECHS['Glarborg'], 642),
#     'NH2_HO2_in_KAUST': ('KAUST', NH3MECHS['KAUST'], 416),
#     'NH2_HO2_in_CRECK': ('CRECK', NH3MECHS['CRECK'], 2014),
#     'NH2_HO2_in_Konnov': ('Konnov', NH3MECHS['Konnov'], 527),
#     'NH2_HO2_in_Faravelli': ('Faravelli', NH3MECHS['Faravelli'], 1996),
#     'NH2_HO2_in_LiHeZhu': ('LiHeZhu', NH3MECHS['LiHeZhu'], 663),
#     'NH2_HO2_in_Okafor': ('Okafor', NH3MECHS['Okafor'], 249),
#     'NH2_HO2_in_Shrestha2022': ('Shrestha2022', NH3MECHS['Shrestha2022'], 1461),
#     'NH2_HO2_in_Zhang2024': ('Zhang2023', NH3MECHS['Zhang2023'], 244),
#     'NH2_HO2_in_Dai2024': ('Dai2024', NH3MECHS['Dai2024'], 115)
# }

reactions_for_analysis = {
    'NH3_OH_in_Faravelli': ('Faravelli', NH3MECHS['Faravelli'], 1993),
    'NH3_OH_in_Shrestha2022': ('Shrestha2022', NH3MECHS['Shrestha2022'], 1459),
    'NH3_OH_in_Zhang2023': ('Zhang2023', NH3MECHS['Zhang2023'], 226),
    'NH3_OH_in_CRECK': ('Zhang2023', NH3MECHS['CRECK'], 2004)
}   

# reactions_for_analysis = {
#     'N2NO_OH_in_Dai2024': ('Dai2024', NH3MECHS['Dai2024'], 236),
#     'N2NO_OH_in_Shrestha2022': ('Shrestha2022', NH3MECHS['Shrestha2022'], 1141),
#     'N2NO_OH_in_Zhang2023': ('Zhang2023', NH3MECHS['Zhang2023'], 243),
#     'N2NO_OH_in_CRECK': ('Zhang2023', NH3MECHS['CRECK'], 2013),
#     'N2NO_OH_in_LiHeZhu': ('LiHeZhu', NH3MECHS['LiHeZhu'], 669),
# }


pressure = 5e5
for label, info in reactions_for_analysis.items():
    mech_label, mech, reaction_number = info
    gas = ct.Solution(mech, 'gas')
    text = (f'Reaction rate constants for R{reaction_number} '
            f'({gas.reactions()[reaction_number-1]}) in mech "{mech_label}" '
            f'at P = {pressure} Pa\n'
            f'((constant {gas.reactions()[reaction_number-1].rate.input_data['rate-constant']}), '
            f'{gas.reactions()[reaction_number-1].rate_coeff_units})')
    output_reaction_rate_constant_temperature_dependence(
        get_reaction_rate_constant_temperature_dependencies(gas, reaction_number, temperatures, pressure),
        label=label,
        info=text
    )
    text = (f'Equilibrium rate constant for R{reaction_number} '
            f'({gas.reactions()[reaction_number-1]}) in mech "{mech_label}" '
            f'at P = {pressure} Pa')
    output_reaction_equilibrium_constant_temperature_dependence(
        get_reaction_equilibrium_constant_temperature_dependence(gas, reaction_number, temperatures, pressure),
        label=label,
        info=text
    )