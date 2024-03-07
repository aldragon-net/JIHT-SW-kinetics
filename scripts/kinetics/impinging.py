from pathlib import Path
import cantera as ct
import numpy as np

print("Launchig")

# parameter values
p = 0.05 * ct.one_atm  # pressure
tburner = 373.0  # burner temperature
tsurf = 500.0

# each mdot value will be solved to convergence, with grid refinement, and
# then that solution will be used for the next mdot
mdot = [0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12]  # kg/m^2/s

rxnmech = 'mechs/Hong2011.yaml'  # reaction mechanism file
comp = 'H2:1.8, O2:1, AR:7'  # premixed gas composition

# The solution domain is chosen to be 20 cm
width = 0.1  # m

loglevel = 1  # amount of diagnostic output (0 to 5)

# Grid refinement parameters
ratio = 3
slope = 0.1
curve = 0.2
prune = 0.06

# Set up the problem
gas = ct.Solution(rxnmech)

# set state to that of the unburned gas at the burner
gas.TPX = tburner, p, comp

# Create the stagnation flow object with a non-reactive surface.  (To make the
# surface reactive, supply a surface reaction mechanism. See example
# catalytic_combustion.py for how to do this.)
sim = ct.ImpingingJet(gas=gas, width=width)

# set the mass flow rate at the inlet
sim.inlet.mdot = mdot[0]

# set the surface state
sim.surface.T = tsurf


zloc = np.array([
    0., 0.1, 0.9, 1
])

tvalues = np.array([
     373.7, 1400, 1400, 500
])


sim.set_grid_min(1e-4)
sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)

sim.set_initial_guess(products='equil')  # assume adiabatic equilibrium products
sim.show()

sim.solve(loglevel, auto=True)

output_path = Path() / "stagnation_flame_data"
output_path.mkdir(parents=True, exist_ok=True)

if "native" in ct.hdf_support():
    output = output_path / "stagnation_flame.h5"
else:
    output = output_path / "stagnation_flame.yaml"
output.unlink(missing_ok=True)

for m, md in enumerate(mdot):
    sim.inlet.mdot = md
    sim.energy_enabled = False
    sim.flame.set_fixed_temp_profile(zloc, tvalues)
    sim.solve(loglevel)
    sim.save(output, name=f"mdot-{m}", description=f"mdot = {md} kg/m2/s")

    # write the velocity, temperature, and mole fractions to a CSV file
    sim.save(output_path / f"stagnation_flame_{m}.csv", basis="mole", overwrite=True)

sim.show_stats()
