from pathlib import Path
import cantera as ct
from configs.burners.impinging_jet_data import ImpingingJetData

from configs.burners.flames import ethylene_flame

from configs.constants import OUTPUT_DIR, MCKENNA_OUTPUT 


class GridRefineCriteria:
    def __init__(self, ratio=3, slope=0.1, curve=0.2, prune=0.05):
        self.ratio = ratio
        self.slope = slope
        self.curve = curve
        self.prune = prune


print("Launching")

flame = ethylene_flame
# parameter values
p = ct.one_atm  # pressure
t_room = 293
t_burner = 400  # burner temperature
t_body = 600.0


rxnmech = 'mechs/GRI/gri30.yaml'  # CRECK/CRECK_2003_TPRF_HT_LT_ALC_ETHERS.yaml'
comp = 'C2H4:1 O2:1.43 N2:5.32 AR:0.068'  # premixed gas composition


loglevel = 1  # amount of diagnostic output (0 to 5)
# Grid refinement parameters
ratio = 3
slope = 0.1
curve = 0.2
prune = 0.06

# Set up the problem
gas = ct.Solution(rxnmech)

# set state to that of the unburned gas at the burner
gas.TPX = flame.T_room, flame.pressure, flame.mixture
md = flame.inlet_velocity * gas.density  # kg/m^2/s
gas.TP = flame.T_burner, flame.pressure

# Create the stagnation flow object with a non-reactive surface
sim = ct.ImpingingJet(gas=gas, width=flame.height)

# set the mass flow rate at the inlet
sim.inlet.mdot = md
sim.energy_enabled = False
sim.flame.set_fixed_temp_profile(
    flame.T_profile.positions, flame.T_profile.temperatures
)

# set the surface state
sim.surface.T = t_body


sim.set_grid_min(1e-4)
sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)

sim.set_initial_guess(products='inlet')  # assume adiabatic equilibrium products

sim.solve(loglevel, auto=False)

output_path = Path() / "output" / "stagnation_flame_data"
output_path.mkdir(parents=True, exist_ok=True)

output = output_path / f"{flame.label}.yaml"
output.unlink(missing_ok=True)

sim.save(output, name=f"{flame.label}_state", description=f"{flame.label}")
# write the velocity, temperature, and mole fractions to a CSV file
sim.save(output_path / f"{flame.label}_X.csv", basis="mole", overwrite=True)
sim.save(output_path / f"{flame.label}_Y.csv", basis="mass", overwrite=True)

sim.show_stats()


def solve_mckenna_stabilized(
        flame: ImpingingJetData,
        mech: str,
        grid_refine_criteria: GridRefineCriteria = GridRefineCriteria(),
        loglevel=1):
    gas = ct.Solution(mech)
    gas.TPX = flame.T_room, flame.pressure, flame.mixture
    md = flame.inlet_velocity * gas.density  # kg/m^2/s
    gas.TP = flame.T_burner, flame.pressure
    # Create the stagnation flow object with a non-reactive surface
    sim = ct.ImpingingJet(gas=gas, width=flame.height)
    sim.inlet.mdot = md
    sim.surface.T = flame.T_body
    sim.set_grid_min(1e-4)
    sim.set_refine_criteria(grid_refine_criteria)

    sim.set_initial_guess(products='equilib')  # assume adiabatic equilibrium products
    if flame.T_profile:
        sim.energy_enabled = False
        sim.flame.set_fixed_temp_profile(
            flame.T_profile.positions, flame.T_profile.temperatures
        )
        sim.solve(loglevel, auto=False)
    else:
        sim.solve(loglevel, auto=True)
    output_path = Path() / OUTPUT_DIR / MCKENNA_OUTPUT
    output_path.mkdir(parents=True, exist_ok=True)
    output = output_path / f"{flame.label}.yaml"
    output.unlink(missing_ok=True)
    sim.save(output, name=f"{flame.label}_state", description=f"{flame.label}")
    # write the velocity, temperature, and mole fractions to a CSV file
    sim.save(output_path / f"{flame.label}_X.csv", basis="mole", overwrite=True)
    sim.save(output_path / f"{flame.label}_Y.csv", basis="mass", overwrite=True)
    sim.show_stats()
    return
