import cantera as ct
import numpy as np
from datetime import datetime

from pathlib import Path

from scripts.utils.mixture import get_trifule_for_air


MECHS = {
    'GRI': 'mechs/GRI/gri30_highT.yaml',
    'CRECK': 'mechs/CRECK/CRECK_2003_TPRF_HT_LT_ALC_ETHERS.yaml'
}

def solve_flame(gas, loglevel=1):
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=5, slope=0.08, curve=0.16, prune=0.005)
    f.transport_model = 'mixture-averaged'
    f.solve(loglevel=1, auto=True)
    output = Path() / 'output/temp/adiabatic_flame.yaml'
    output.unlink(missing_ok=True)
    # Solve with the energy equation enabled
    f.save(output, name="mix", description="solution with mixture-averaged transport")
    # print(f"mixture-averaged flamespeed = {f.velocity[0]:7f} m/s")
    # # Solve with multi-component transport properties
    # f.transport_model = 'multicomponent'
    # f.solve(loglevel)  # don't use 'auto' on subsequent solves
    # print(f"multicomponent flamespeed = {f.velocity[0]:7f} m/s")
    # f.save(output, name="multi", description="solution with multicomponent transport")
    # # write the velocity, temperature, density, and mole fractions to a CSV file
    # f.save('output/temp/current_flame.csv', basis="mole", overwrite=True)
    Su0 = f.velocity[0]
    print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")
    return f


temperature = 300
pressure = ct.one_atm
compression = 10
mech_label = 'GRI'
width = 0.05
alphas = [0, 5, 10, 20, 30, 50, 70, 80, 90, 95, 100]
beta = 0

gas = ct.Solution(MECHS[mech_label], 'gri30')

with open('output/flame/report.out', 'a') as f:
    f.write(f'\nSequence started at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
    f.write(f'Mechanism is {mech_label}\n')
for alpha in alphas:
    mixture = get_trifule_for_air('CH4', 'H2', 'CH3OH', alpha, beta)
    mixture_label = f'M{alpha}_{beta}'
    gas.X = mixture
    gas.TP = temperature, pressure
    gas.SV = gas.s, gas.v/compression
    T_comp = gas.T
    P_comp = gas.P
    t0 = datetime.now()
    flame = solve_flame(gas, loglevel=1)
    t1 = datetime.now()
    velocity = flame.velocity[0]
    flame.save(f'output/flame/{mech_label}_{mixture_label}.csv', basis="mole", overwrite=True)
    with open('output/flame/report.out', 'a') as f:
        f.write(f'M{alpha}_{beta} ({mixture}) P={P_comp/1e5:.3f} bar T = {T_comp:.0f} K v = {velocity*100:.2f} cm/s (computed in {str(t1 - t0)})\n')
with open('output/flame/report.out', 'a') as f:
    f.write(f'Sequence completed at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')