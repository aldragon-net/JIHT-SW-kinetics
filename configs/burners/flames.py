import cantera as ct
from configs.burners.tprofile import read_temperature_profile


class FlameData:
    def __init__(self, label, mixture, inlet_velocity, height, T_burner, T_body,
                 T_profile, T_room=293, pressure=ct.one_atm):
        self.label = label
        self.mixture = mixture
        self.pressure = pressure
        self.inlet_velocity = inlet_velocity
        self.height = height
        self.T_room = T_room
        self.T_burner = T_burner
        self.T_body = T_body
        self.T_profile = T_profile


ethylene_flame = FlameData(
    label="ethylene_base",
    mixture='C2H4:1 O2:1.43 N2:5.32 AR:0.068',
    inlet_velocity=0.06,
    height=0.023,
    T_room=293,
    T_burner=400,
    T_body=600,
    T_profile=read_temperature_profile('mckenna-noflame.dat')
)
