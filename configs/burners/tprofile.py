from os import path


class TemperatureProfile:
    def __init__(self, positions, temperatures):
        self.positions = positions
        self.temperatures = temperatures


def read_temperature_profile(filename):
    filepath = path.join(path.dirname(__file__), filename)
    z_locs = []
    t_values = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()
        if line:
            z_loc, t_value = line.strip().split()
            z_locs.append(float(z_loc))
            t_values.append(float(t_value))
    relative_z_locs = [z_loc / z_locs[-1] for z_loc in z_locs]
    profile = TemperatureProfile(relative_z_locs, t_values)
    return profile


temperature_profile = read_temperature_profile('mckenna-noflame.dat')
