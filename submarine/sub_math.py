"""
Prelimenary calculations for an autonomous survey submarine
"""

import matplotlib.pyplot as plt

FT2_TO_M2 = 0.09290304  # 1 ft^2 = 0.09290304 m^2
FT3_TO_M3 = 0.02831685  # 1 ft^3 = 0.02831685 m^3
P_ATM = 101325          # Standard atmospheric pressure in Pa (1 atm)

# sample cross sectional areas in ft^2
csa_1 = 12 * 18  # 216 ft^2
csa_2 = 16 * 16  # 256 ft^2
csa_3 = 18 * 18  # 324 ft^2

# sample submarine hull length in ft
sub_length_ft = 5.0  # 5 feet long

# gauge pressure formula: rho*g*h (Pa)
def pressure(rho, g, h):
    return rho * g * h

# absolute pressure formula: rho*g*h + P_atm (Pa)
def absolute_pressure(rho, g, h, p_atm=P_ATM):
    return pressure(rho, g, h) + p_atm

# force due to pressure: P * A (N)
def pressure_force(pressure_pa, csa_ft2):
    return pressure_pa * (csa_ft2 * FT2_TO_M2)

# buoyancy force: F_b = rho * g * V_displaced (N)
def buoyant_force(rho, g, volume_ft3):
    volume_m3 = volume_ft3 * FT3_TO_M3
    return rho * g * volume_m3

# neutral buoyancy mass: m = rho * V_displaced (kg)
def neutral_buoyancy_mass(rho, volume_ft3):
    volume_m3 = volume_ft3 * FT3_TO_M3
    return rho * volume_m3

# net buoyancy force: F_net = F_b - W_sub (N)
def net_buoyancy(f_buoyant, weight_sub_n):
    return f_buoyant - weight_sub_n

# relevant params
rho_water = 1025 # kg/m³ (seawater)
g = 9.81         # m/s^2
h = 0            # m

# Buoyancy calculations for CSA 1 (Volume = CSA * length)
v_displaced_ft3 = csa_1 * sub_length_ft
f_buoyancy = buoyant_force(rho_water, g, v_displaced_ft3)
m_neutral = neutral_buoyancy_mass(rho_water, v_displaced_ft3)

print("=== BUOYANCY ANALYSIS ===")
print(f"Submarine Length:              {sub_length_ft} ft")
print(f"Displaced Volume (CSA 1):      {v_displaced_ft3:.2f} ft³ ({v_displaced_ft3 * FT3_TO_M3:.2f} m³)")
print(f"Total Buoyant Force:           {f_buoyancy:.2f} N ({f_buoyancy/1000:.2f} kN)")
print(f"Required Neutral Mass:        {m_neutral:.2f} kg ({m_neutral * 2.20462:.2f} lbs)\n")

# lists to store results for graphing
depth_list = []
gauge_p_list = []
abs_p_list = []
hydro_force_list = []
crush_force_list = []

print("=== PRESSURE & FORCE VS DEPTH ===")
# calculations at different pressures from 0-50m
while h <= 50:
    p_gauge = pressure(rho_water, g, h)
    p_abs = absolute_pressure(rho_water, g, h)
    
    f_hydro = pressure_force(p_gauge, csa_1)
    f_crush = pressure_force(p_abs, csa_1)
    
    depth_list.append(h)
    gauge_p_list.append(p_gauge / 1000.0)    # kPa
    abs_p_list.append(p_abs / 1000.0)        # kPa
    hydro_force_list.append(f_hydro / 1000.0) # kN
    crush_force_list.append(f_crush / 1000.0) # kN
    
    if h % 10 == 0:  # print every 10 meters for conciseness
        print(f"Depth: {h} m")
        print(f"  Gauge Pressure:    {p_gauge:.2f} Pa ({p_gauge/1000:.2f} kPa)")
        print(f"  Absolute Pressure: {p_abs:.2f} Pa ({p_abs/1000:.2f} kPa)")
        print(f"  Hydrostatic Force: {f_hydro:.2f} N ({f_hydro/1000:.2f} kN)")
        print(f"  Crushing Force:    {f_crush:.2f} N ({f_crush/1000:.2f} kN)\n")
    
    h += 1

# graph pressure (kPa) and force (kN) vs. depth (m)
def graph_pressure_force(depths, hydro_forces, crush_forces):
    fig, ax1 = plt.subplots()

    color = 'tab:blue'
    ax1.set_xlabel('Depth (m)')
    ax1.set_ylabel('Hydrostatic Force (kN)', color=color)
    ax1.plot(depths, hydro_forces, color=color, label='Hydrostatic Force (Net)')
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  
    color = 'tab:red'
    ax2.set_ylabel('Crushing Force (kN)', color=color)
    ax2.plot(depths, crush_forces, color=color, linestyle='--', label='Crushing Force (Total)')
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()
    plt.title('Submarine Force vs. Depth (CSA = 216 ft^2)')
    plt.show()

graph_pressure_force(depth_list, hydro_force_list, crush_force_list)
