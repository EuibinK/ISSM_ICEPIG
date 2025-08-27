import numpy as np

def gbsh2o(temperature, phi=0.01):
    """GBSH2O - figure out the rigidity of H2O ice for a given temperature.

    Usage:
        rigidity = gbsh2o(temperature, phi)
        rigidity = gbsh2o(temperature) # uses default phi = 0.01
    """

    # Declaring temperature and rigidity arrays
    if np.ndim(temperature) == 2:
        T = temperature.flatten()
    elif isinstance(temperature, float) or isinstance(temperature, int):
        T = np.array([temperature])
    else:
        T = np.array(temperature)

    # Check for negative temperatures
    if np.any(T < 0):
        raise ValueError('input temperature should be in Kelvin (positive)')

    # Beyond-melting-point case
    if np.any(T >= 273.15):
        melting_indices = np.where(T >= 273.15)[0]
        print(f'Warning: gbsh2o.py: H2O ICE - GUARANTEED MELTING. Some temperatures are beyond 273.15K.')
        print(f'Look at indices: {melting_indices}')


    # Coefficients
    Rg = 8.3144598    # J mol^-1 K^-1
    n = 1.8           # Glen's exponent

    rd = 1.5e-6
    f = 0.25
    A0 = 3.9e-3       # s^-1 MPa
    p = -1.4

    d = (13.4e-6 / np.sqrt(f)) * (rd / 1.5e-6) * ((0.05 / phi)**0.5)
    A1 = A0 * 0.5 * (3**((n + 1) / 2))
    A2 = A1 * np.exp(-2 * n * phi) * (d**p)

    Q = 49000.       # J mol^-1

    # Arrhenius Law
    A = A2 * np.exp(-Q / (T * Rg))    # s^-1 MPa
    rigidity = 1e6 * (A**(-1/n))     # s^(1/n) Pa

    # Return scalar if input was scalar
    if isinstance(temperature, (float, int)):
        return rigidity[0]

    return rigidity

