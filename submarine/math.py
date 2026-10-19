"""
Prelimenary calculations for an autonomous survey submarine
"""

IN_TO_M2 = 0.00064516

# csa in in^2
csa_1 = 12 * 18 
csa_2 = 16 * 16
csa_3 = 18 * 18

def pressure(rho, g, h):
    return rho * g * h

def pressure_force(pressure, csa, inches=True):
    return pressure * csa * IN_TO_M2 if inches else pressure * csa


rho_water = 1025 # kg/m³
g = 9.81 # m/s^2
h = 0 # m


while h <= 50:
    press = pressure(rho_water, g, h)
    force = pressure_force(press, csa_1, inches=True)
    print(f"Pressure at {h} m: {press} Pa")
    print(f"Force at {h} m on CSA {csa_1} in^2: {force}")
    h += 1
