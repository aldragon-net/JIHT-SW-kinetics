import cantera as ct


class ImpingingJetData:
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
