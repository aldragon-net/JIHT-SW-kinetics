STOICHS = {
    '': 0,
    'O2': 1,
    'H2': 2,
    'CH4': 0.5,
    'CO': 2,
    'NH3': 4/3,
    'CH3OH': 2/3,
    'CH3OCH3': 1/3
}

admixtures = ['NH3', 'CH3OH', 'CH3OCH3']

PERCENT = 100
o2_fraction = 7
phi = 1.0

def fraction_string(species, fraction):
    return f' {species}:{fraction:.3f}' if fraction > 0 else ''


def get_trifuel_for_o2(o2_fraction, primary, secondary, tertiary, alpha, beta, phi=1.0, diluter='AR', TOTAL=100):
    alpha = alpha/TOTAL
    beta = beta/TOTAL
    if not all([species in STOICHS for species in [primary, secondary, tertiary]]):
        raise IndexError
    primary_fraction = ((1-alpha) * o2_fraction * STOICHS[primary]/phi) * (1 - beta)
    secondary_fraction = (alpha * o2_fraction * STOICHS[secondary]/phi) * (1 - beta)
    tertiary_fraction = o2_fraction * beta * STOICHS[tertiary]/phi
    diluter_fraction = TOTAL - o2_fraction - secondary_fraction - primary_fraction - tertiary_fraction
    result = ''.join([fraction_string(s, f) for (s, f) in (
        (primary, primary_fraction),
        (secondary, secondary_fraction),
        (tertiary, tertiary_fraction),
        ('O2', o2_fraction),
        (diluter, diluter_fraction)
    )]).strip()
    return result


def get_trifule_for_air(primary, secondary, tertiary, alpha, beta, phi=1.0):
    alpha = alpha/100
    beta = beta/100
    N2_TO_O2 = 3.77
    o2_fraction = 1
    n2_fraction = o2_fraction * N2_TO_O2
    if not all([species in STOICHS for species in [primary, secondary, tertiary]]):
        raise IndexError
    primary_fraction = ((1-alpha) * o2_fraction * STOICHS[primary]/phi) * (1 - beta)
    secondary_fraction = (alpha * o2_fraction * STOICHS[secondary]/phi) * (1 - beta)
    tertiary_fraction = o2_fraction * beta * STOICHS[tertiary]/phi
    total = sum((primary_fraction, secondary_fraction, tertiary_fraction, o2_fraction, n2_fraction))
    primary_fraction, secondary_fraction, tertiary_fraction, o2_fraction, n2_fraction = [
        100*fraction/total for fraction in [
            primary_fraction, secondary_fraction, tertiary_fraction, o2_fraction, n2_fraction
        ]
    ]
    result = ''.join([fraction_string(s, f) for (s, f) in (
        (primary, primary_fraction),
        (secondary, secondary_fraction),
        (tertiary, tertiary_fraction),
        ('O2', o2_fraction),
        ('N2', n2_fraction)
    )]).strip()
    return result

