import cantera as ct
from math import pi
from configs.burners.mckenna import mckenna_jiht_lab19


def get_inlet_mixture(fuel_mix, fuel_flow, oxydant_mix, oxydant_flow):
    gas = ct.Solution('gri30.yaml')
    fuel = ct.Quantity(gas, constant='HP')
    fuel.TPX = 300.0, ct.one_atm, fuel_mix
    fuel.moles = fuel_flow
    oxydant = ct.Quantity(gas, constant='HP')
    oxydant.TPX = 300.0, ct.one_atm, oxydant_mix
    oxydant.moles = oxydant_flow
    m = fuel + oxydant
    mixture = ' '.join([
        f'{species}:{float(fraction):.6f}' for species, fraction in m.mole_fraction_dict().items()
    ])
    return mixture


class ImpingingJetData:
    def __init__(self, label, fuel, fuel_flow_lph, oxydizer, oxydizer_flow_lph, height, T_burner, 
                 T_body, T_profile, T_room=293, pressure=ct.one_atm, burner=mckenna_jiht_lab19):
        self.label = label
        self.mixture = get_inlet_mixture(
            fuel_mix=fuel, fuel_flow=fuel_flow_lph,
            oxydant_mix=oxydizer, oxydant_flow=oxydizer_flow_lph)
        self.pressure = pressure
        self.inlet_velocity = (
            ((fuel_flow_lph + oxydizer_flow_lph) / (3600*1000)) / (pi*(burner.diameter**2)/4)
        )
        self.height = height
        self.T_room = T_room
        self.T_burner = T_burner
        self.T_body = T_body
        self.T_profile = T_profile
