from pathlib import Path
from datetime import datetime
import cantera as ct
from configs.burners.impinging_jet_data import ImpingingJetData

from configs.burners.flames import acetylene_flame, ethylene_flame, dme_flames

from configs.constants import OUTPUT_DIR, MCKENNA_OUTPUT

from scripts.utils.csv_processors import postprocess_creck_soot


class GridRefineCriteria:
    def __init__(self, ratio=3, slope=0.1, curve=0.2, prune=0.05):
        self.ratio = ratio
        self.slope = slope
        self.curve = curve
        self.prune = prune


flame = acetylene_flame


def get_report(flame, mech, grid_refine_criteria):
    return (
        f'SOLVING McKenna Stabilized at {datetime.now()}\n'
        f'Flame "{flame.label}" with mixture composition\n\t{flame.mixture}\n'
        f'using kinetic mechanism {mech}\n'
        f'with parameters: inlet velocity = {flame.inlet_velocity*100:.3f} cm/s\n'
        f'with grid refine criteria: ratio={grid_refine_criteria.ratio:.2f}, '
        f'slope={grid_refine_criteria.slope:.3f}, '
        f'curve={grid_refine_criteria.curve:.3f}, '
        f'prune={grid_refine_criteria.prune:.3f}\n'
        )


def solve_mckenna_stabilized(
        flame: ImpingingJetData,
        mech: str,
        grid_refine_criteria: GridRefineCriteria = GridRefineCriteria(),
        loglevel: int = 1):
    output_path = Path() / OUTPUT_DIR / MCKENNA_OUTPUT / flame.label
    output_path.mkdir(parents=True, exist_ok=True)
    report = get_report(flame, mech, grid_refine_criteria)
    print(report)
    with open(output_path / 'report.rxt', 'w') as f:
        f.write(report)
    gas = ct.Solution(mech)
    gas.TPX = flame.T_room, flame.pressure, flame.mixture
    md = flame.inlet_velocity * gas.density  # kg/m^2/s
    gas.TP = flame.T_burner, flame.pressure
    # Create the stagnation flow object with a non-reactive surface
    sim = ct.ImpingingJet(gas=gas, width=flame.height)
    sim.inlet.mdot = md
    sim.surface.T = flame.T_body
    sim.set_grid_min(1e-4)
    sim.set_refine_criteria(
        ratio=grid_refine_criteria.ratio,
        slope=grid_refine_criteria.slope,
        curve=grid_refine_criteria.curve,
        prune=grid_refine_criteria.prune
    )
    sim.set_initial_guess(products='equil')  # assume adiabatic equilibrium products
    if flame.T_profile:
        sim.energy_enabled = False
        sim.flame.set_fixed_temp_profile(
            flame.T_profile.positions, flame.T_profile.temperatures
        )
        sim.solve(loglevel, auto=False)
    else:
        sim.solve(loglevel, auto=True)
    output = output_path / f"{flame.label}.yaml"
    output.unlink(missing_ok=True)
    sim.save(output, name=f"{flame.label}_state", description=f"{flame.label}")
    # write the velocity, temperature, and mole fractions to a CSV file
    sim.save(output_path / f"{flame.label}_X.csv", basis="mole", overwrite=True)
    sim.save(output_path / f"{flame.label}_Y.csv", basis="mass", overwrite=True)
    if loglevel > 0: 
        sim.show_stats()
    with open(output_path / 'report.rxt', 'a') as f:
        f.write(f'Solved at {datetime.now()}')
    return


def multi_solve_mckenna_stabilized(
        flames: list[ImpingingJetData],
        mech: str,
        grid_refine_criteria: GridRefineCriteria = GridRefineCriteria(),
        loglevel: int = 1):
    print(f'Started processing set of {len(flames)} flames at {datetime.now()}')
    path = Path() / OUTPUT_DIR / MCKENNA_OUTPUT
    for i, flame in enumerate(flames):
        print(f'Task {i+1} of {len(flames)}')
        solve_mckenna_stabilized(flame, mech, grid_refine_criteria, loglevel)
        postprocess_creck_soot(path, flame)
    print(f'Finished at {datetime.now()}')


flames = [ethylene_flame, acetylene_flame]

rxnmech = 'mechs/GRI/gri30.yaml'#'mechs/CRECK/CRECK_2003_TOT_HT_LT_SOOT.yaml'

multi_solve_mckenna_stabilized(
    flames=dme_flames,
    mech=rxnmech,
    grid_refine_criteria=GridRefineCriteria(),
    loglevel=0
)
