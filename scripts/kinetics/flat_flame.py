import cantera as ct
import numpy as np
from datetime import datetime

from pathlib import Path

from scripts.utils.mixture import get_trifule_for_air




MECHS = {
    'GRI': 'mechs/GRI/gri30_highT.yaml',
    'CRECK': 'mechs/CRECK/CRECK_2003_TPRF_HT_LT_ALC_ETHERS.yaml'
}


gas = ct.Solution(MECHS['GRI'])

width = 0.05


alphas = [0, 5, 10, 20, 30, 50, 70, 80, 90, 95, 100]
beta = 0


def solve_flame(gas, loglevel=1):
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=5, slope=0.05, curve=0.1, prune=0.005)
    f.set_refine_criteria()
    # Solve with mixture-averaged transport model
    f.transport_model = 'mixture-averaged'
    f.solve(loglevel=1, auto=True)

    output = Path() / 'output/temp/adiabatic_flame.yaml'
    output.unlink(missing_ok=True)

    # Solve with the energy equation enabled
    f.save(output, name="mix", description="solution with mixture-averaged transport")
    print(f"mixture-averaged flamespeed = {f.velocity[0]:7f} m/s")
    # Solve with multi-component transport properties
    f.transport_model = 'multicomponent'
    f.solve(loglevel)  # don't use 'auto' on subsequent solves
    print(f"multicomponent flamespeed = {f.velocity[0]:7f} m/s")
    f.save(output, name="multi", description="solution with multicomponent transport")
    # write the velocity, temperature, density, and mole fractions to a CSV file
    f.save('output/temp/adiabatic_flame.csv', basis="mole", overwrite=True)
    Su0 = f.velocity[0]
    print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")
    return f


temperature = 300
pressure = ct.one_atm

for alpha in alphas:
    mixture = get_trifule_for_air('CH4', 'H2', '', alpha, 0)
    gas.X = mixture
    gas.TP = temperature, pressure
    gas.SV = gas.s, gas.v/10
    print(mixture, gas.P/ct.one_atm, gas.T)

with open('output/flame_report.out', 'a') as f:
    f.write(f'Sequence started at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
for alpha in alphas:
    mixture = get_trifule_for_air('CH4', 'H2', 'CH3OH', alpha, 50)
    gas.X = mixture
    gas.TP = temperature, pressure
    flame = solve_flame(gas, loglevel=1)
    velocity = flame.velocity[0]
    with open('output/flame_report.out', 'a') as f:
        f.write(f'for alpha={alpha} beta={beta} ({mixture}) flame velocity is {velocity*100:.2f} cm/s\n')