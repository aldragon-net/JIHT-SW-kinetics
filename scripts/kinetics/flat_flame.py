import cantera as ct
import numpy as np

from pathlib import Path

# Inlet Temperature in Kelvin and Inlet Pressure in Pascals
# In this case we are setting the inlet T and P to room temperature conditions
temperature = 300
pressure = 1013250

# Define the gas-mixutre and kinetics
# In this case, we are choosing a GRI3.0 gas
gas = ct.Solution('mechs/CRECK/CRECK_2003_TPRF_HT_LT_ALC_ETHERS.yaml')
# gas = ct.Solution('mechs/GRI/gri30_highT.yaml')

width = 0.03

gas.X = 'CH4:3.5 O2:7 N2:0'
gas.TP = temperature, pressure

f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=5, slope=0.1, curve=0.2, prune=0.01)
f.set_refine_criteria()

loglevel = 1

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
f.save('adiabatic_flame.csv', basis="mole", overwrite=True)

Su0 = f.velocity[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")
