"""
CNC mill thermal / mechanical / vibration simulation suite.

This module is a fast, physics-based design explorer for a desktop sub-10-micron
precision CNC mill. Think of it as a "spreadsheet with physics": it gives
order-of-magnitude answers quickly, not final machine-shop numbers.

What it can do:
    1. Estimate how much a tool bends under cutting loads.
    2. Estimate cutting forces for different materials.
    3. Predict the tool's natural frequency and check for resonance.
    4. Compare machine-base materials for stiffness, damping, and cost.
    5. Sweep tool sizes and base materials to find combos that meet the 10 µm
       precision target.
    6. Optionally talk to a MATLAB/Simulink reference model.

Plain-English glossary used in the code:
    - Workpiece material: the block you are cutting.
    - Tool material: the spinning cutter (carbide, high-speed steel, etc.).
    - Base material: the machine frame or bed that holds everything.
    - Stick-out / overhang: how far the tool sticks out from the spindle.
    - Feed per tooth (f_z): how far the cutter advances each time a tooth bites.
    - Axial depth (a_p): how deep the cut is along the tool axis.
    - Tooth-passing frequency: how often teeth hit the workpiece per second.
    - Natural frequency: the speed at which the tool wants to vibrate.
    - Damping ratio (ζ): how quickly vibrations die out.
    - Dynamic amplification: how much worse vibration gets near resonance.

Terminology:
    - workpiece material is looked up in `material_properties` and
      `cutting_coefficients`.
    - tool material is looked up in `tool_material_properties`.
    - base material is looked up in `base_material_properties`.

Limitations:
    - Tool is modeled as a solid cylinder (no flutes, taper, collet compliance).
    - Base is modeled as a simple cantilever beam, not a real machine structure.
    - Cutting coefficients are rough averages, not calibrated to micro-milling.
    - Thermal model is a very simple lumped estimate.
    - Chatter prediction is a coarse flag, not a full stability-lobe analysis.
"""

# ---------------------------------------------------------------------------
# Imports and optional MATLAB Engine setup
# ---------------------------------------------------------------------------
import math
import os

# Pre-load Python's XML parser before matlab.engine starts. On macOS the MATLAB
# Engine runtime can shadow Conda Python's libexpat with MATLAB's own copy,
# causing symbol errors when Simulink later parses model XML. Loading pyexpat
# early binds it to Python's library before MATLAB changes the loader state.
import xml.etree.ElementTree  # noqa: F401

try:
    import pyexpat  # noqa: F401
except Exception:
    pass

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# matlab.engine is optional; the Python-native simulation works without it.
# On macOS, keep DYLD_LIBRARY_PATH unchanged. Prepending MATLAB's bin/maca64
# globally causes dylib conflicts (e.g. MATLAB's libexpat shadowing Python's).
# With MATLAB installed in the default /Applications/MATLAB_R####a.app path,
# matlabengine finds it automatically. If it cannot, run:
#     python setup_matlab_bridge.py
# and follow the instructions shown there.
try:
    import matlab.engine
except Exception as _matlab_import_error:
    _MATLAB_IMPORT_ERROR = str(_matlab_import_error)
else:
    _MATLAB_IMPORT_ERROR = None


# ---------------------------------------------------------------------------
# Global constants and defaults
# ---------------------------------------------------------------------------
MICRON_SCALE = 1e-6
MIN_DEFLECTION = 10 * MICRON_SCALE  # 10 µm, expressed in meters
MAX_SPINDLE_RPM = 95000  # rpm
MIN_SPINDLE_RPM = 50000  # rpm
MAX_CUTTING_LOAD = 300  # N
MIN_CUTTING_LOAD = 10  # N

DEFAULT_SPINDLE_RPM = 60000
DEFAULT_FEED_PER_TOOTH = 2e-6  # m/tooth
DEFAULT_AXIAL_DEPTH = 0.0001  # m
DEFAULT_STOP_TIME = 0.1  # s
DEFAULT_N_FLUTES = 2
DEFAULT_DAMPING_RATIO = 0.05
DEFAULT_TOOL_MATERIAL = "wc_co"
HEAT_PARTITION = 0.1  # fraction of cutting heat entering the tool

SIMULINK_MODEL_NAME = "cnc_mill_sim"


# ---------------------------------------------------------------------------
# Material databases
# ---------------------------------------------------------------------------
# Workpiece materials: the part being cut.
# E = Young's modulus (Pa), nu = Poisson's ratio (dimensionless).
# Values are typical room-temperature engineering averages.
material_properties = {
    "aluminum": {"rho": 2700.0,  "c_p": 900.0, "k": 237.0, "cost per kg": 3.26,    "E": 70e9,  "nu": 0.33},
    "brass":    {"rho": 8530.0,  "c_p": 380.0, "k": 120.0, "cost per kg": 7.00,    "E": 100e9, "nu": 0.35},
    "copper":   {"rho": 8960.0,  "c_p": 385.0, "k": 398.0, "cost per kg": 13.00,   "E": 119e9, "nu": 0.34},
    "gold":     {"rho": 19320.0, "c_p": 129.0, "k": 315.0, "cost per kg": 140225.00, "E": 77e9,  "nu": 0.42},
    "silver":   {"rho": 10500.0, "c_p": 235.0, "k": 429.0, "cost per kg": 2050.00,  "E": 76e9,  "nu": 0.37},
    "steel":    {"rho": 7850.0,  "c_p": 490.0, "k": 45.0,  "cost per kg": 2.50,    "E": 210e9, "nu": 0.29},
    "titanium": {"rho": 4500.0,  "c_p": 520.0, "k": 17.0,  "cost per kg": 35.00,   "E": 110e9, "nu": 0.34},
    "tungsten": {"rho": 19300.0, "c_p": 134.0, "k": 173.0, "cost per kg": 138.00,  "E": 400e9, "nu": 0.28},
}

# ---------------------------------------------------------------------------
# tool material properties
# ---------------------------------------------------------------------------
# These properties describe the cutting tool itself. They determine the tool's
# bending stiffness, natural frequency, and thermal response.
tool_material_properties = {
    "hss":   {"rho": 8200.0,  "c_p": 460.0, "k": 24.0,  "E": 210e9, "nu": 0.30},  # high-speed steel
    "wc_co": {"rho": 15160.0, "c_p": 250.0, "k": 80.0,  "E": 600e9, "nu": 0.24},  # tungsten carbide / cobalt
}

def build_materials_export_table():
    """Return a single flat table of all materials and their properties.

    This combines workpiece, tool, and base materials into one list of
    dictionaries, each with the same set of columns. Missing values are filled
    with None. Useful for exporting to CSV/Google Sheets.

    Columns:
        category, material, rho, c_p, k, E, nu, cost_per_kg, damping_ratio
    """
    columns = [
        "category", "material", "rho", "c_p", "k", "E", "nu",
        "cost_per_kg", "damping_ratio", "notes"
    ]

    rows = []
    for material, props in material_properties.items():
        rows.append({
            "category": "workpiece",
            "material": material,
            "rho": props.get("rho"),
            "c_p": props.get("c_p"),
            "k": props.get("k"),
            "E": props.get("E"),
            "nu": props.get("nu"),
            "cost_per_kg": props.get("cost per kg"),
            "damping_ratio": None,
            "notes": "",
        })

    for material, props in tool_material_properties.items():
        rows.append({
            "category": "tool",
            "material": material,
            "rho": props.get("rho"),
            "c_p": props.get("c_p"),
            "k": props.get("k"),
            "E": props.get("E"),
            "nu": props.get("nu"),
            "cost_per_kg": None,
            "damping_ratio": None,
            "notes": "",
        })

    for material, props in base_material_properties.items():
        rows.append({
            "category": "base",
            "material": material,
            "rho": props.get("rho"),
            "c_p": None,
            "k": None,
            "E": props.get("E"),
            "nu": props.get("nu"),
            "cost_per_kg": props.get("cost per kg"),
            "damping_ratio": props.get("damping_ratio"),
            "notes": props.get("notes", ""),
        })

    # Ensure every row has the same columns, in order.
    return [{col: row.get(col) for col in columns} for row in rows]


# ---------------------------------------------------------------------------
# machine base / frame material properties
# ---------------------------------------------------------------------------
# These describe the machine bed or structural loop that supports the spindle
# and/or workpiece. The base is modeled as a simple cantilever beam of
# rectangular cross-section; its stiffness and natural frequency are compared
# across common base materials.
base_material_properties = {
    # damping_ratio is a typical structural-damping estimate for each material.
    # Values are order-of-magnitude room-temperature averages from machinery
    # design handbooks and recent composite literature.
    "aluminum":      {
        "rho": 2700.0, "E":  70e9, "nu": 0.33, "cost per kg": 3.26,
        "damping_ratio": 0.0010,
        "notes": "light, stiff, low damping",
    },
    "steel":         {
        "rho": 7850.0, "E": 210e9, "nu": 0.29, "cost per kg": 2.50,
        "damping_ratio": 0.0015,
        "notes": "highest stiffness, low damping",
    },
    "cast_iron":     {
        "rho": 7300.0, "E": 110e9, "nu": 0.27, "cost per kg": 1.80,
        "damping_ratio": 0.0050,
        "notes": "good stiffness + damping balance",
    },
    "granite":       {
        "rho": 2700.0, "E":  50e9, "nu": 0.25, "cost per kg": 0.50,
        "damping_ratio": 0.0100,
        "notes": "cheaper, moderate damping",
    },
    "concrete":      {
        "rho": 2400.0, "E":  30e9, "nu": 0.18, "cost per kg": 0.10,
        "damping_ratio": 0.0150,
        "notes": "very cheap, good damping, low stiffness",
    },
    "epoxy_granite": {
        "rho": 2300.0, "E":  40e9, "nu": 0.25, "cost per kg": 2.00,
        "damping_ratio": 0.0202,
        "notes": "polymer-concrete composite; E and ζ from Abdellah et al. 2026, MDPI Polymers",
    },
}

# Representative base geometry: a vertical cantilever column / bridge beam.
DEFAULT_BASE_WIDTH = 0.10   # m
DEFAULT_BASE_HEIGHT = 0.10  # m
DEFAULT_BASE_LENGTH = 0.30  # m (overhang / unsupported length)
DEFAULT_BASE_DAMPING_RATIO = 0.05
# K_tc, K_rc, K_ac = specific cutting pressures in the tangential, radial and
# axial directions (Pa). K_te, K_re, K_ae = edge/ploughing coefficients (N/m).
# These are rough room-temperature averages for a sharp uncoated tool cutting
# the listed workpiece material.
cutting_coefficients = {
    "aluminum": {"K_tc": 8.0e8, "K_rc": 3.0e8, "K_ac": 1.5e8, "K_te": 1.2e5, "K_re": 5.0e4, "K_ae": 2.5e4},
    "brass":    {"K_tc": 1.0e9, "K_rc": 4.0e8, "K_ac": 2.0e8, "K_te": 1.5e5, "K_re": 7.0e4, "K_ae": 3.5e4},
    "copper":   {"K_tc": 1.2e9, "K_rc": 4.5e8, "K_ac": 2.0e8, "K_te": 1.8e5, "K_re": 8.0e4, "K_ae": 4.0e4},
    "gold":     {"K_tc": 7.0e8, "K_rc": 2.8e8, "K_ac": 1.4e8, "K_te": 1.0e5, "K_re": 4.5e4, "K_ae": 2.2e4},
    "silver":   {"K_tc": 7.5e8, "K_rc": 3.0e8, "K_ac": 1.5e8, "K_te": 1.1e5, "K_re": 4.8e4, "K_ae": 2.4e4},
    "steel":    {"K_tc": 1.8e9, "K_rc": 7.0e8, "K_ac": 3.5e8, "K_te": 2.5e5, "K_re": 1.2e5, "K_ae": 6.0e4},
    "titanium": {"K_tc": 2.5e9, "K_rc": 1.0e9, "K_ac": 5.0e8, "K_te": 3.5e5, "K_re": 1.7e5, "K_ae": 8.5e4},
    "tungsten": {"K_tc": 3.5e9, "K_rc": 1.4e9, "K_ac": 7.0e8, "K_te": 5.0e5, "K_re": 2.5e5, "K_ae": 1.2e5},
}


# graph material properties
def plot_material_properties(disp=True, save=False):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    
    materials = list(material_properties.keys())
    
    rho = [props["rho"] for props in material_properties.values()]
    c_p = [props["c_p"] for props in material_properties.values()]
    k = [props["k"] for props in material_properties.values()]
    
    ax.scatter(rho, c_p, k, s=100)  # type: ignore[arg-type]

    # label each point with its material name
    for i, material in enumerate(materials):
        ax.text(rho[i], c_p[i], k[i], material, fontsize=8)

    ax.set_xlabel("Density (kg/m^3)")
    ax.set_ylabel("Specific Heat (J/kg*K)")
    ax.set_zlabel("Thermal Conductivity (W/m*K)")
    ax.set_title("CNC Mill Material Thermal Properties")
    
    if save:
        plt.savefig("material_properties.png")
    if disp:
        plt.show()
    if not disp and not save:
        plt.close()


def plot_material_costs(disp=True, save=False):
    fig, ax = plt.subplots()

    materials = list(material_properties.keys())
    costs = [props["cost per kg"] for props in material_properties.values()]
    normalized_costs = [math.log10(cost) for cost in costs]

    ax.bar(materials, normalized_costs)
    ax.set_xlabel("Material")
    ax.set_ylabel("Cost per kg ($) (log_10 scale)")
    ax.set_title("CNC Mill Material Costs (normalized)")

    if save:
        plt.savefig("material_costs.png")
    if disp:
        plt.show()
    if not disp and not save:
        plt.close()


def plot_diameter_vs_deflection(force=26.65, tool_material="wc_co", disp=False, save=True):
    """Plot tool diameter vs. tip deflection for various stick-out lengths.

    Parameters
    ----------
    force : float
        Cutting force in N (default 26.65 N, representative peak cutting force).
    tool_material : str
        Tool material key (default 'wc_co' for tungsten carbide).
    disp : bool
        Whether to display the plot interactively.
    save : bool
        Whether to save the plot to 'tool_diameter_vs_deflection.png'.
    """
    diameters_mm = np.geomspace(0.1, 5.0, 300)
    lengths_mm = [5.0, 8.0, 10.0, 15.0]

    plt.figure(figsize=(9, 6))

    for L_mm in lengths_mm:
        deflections_um = []
        for D_mm in diameters_mm:
            D = D_mm * 1e-3
            L = L_mm * 1e-3
            delta_m = deflection(tool_material, D, L, force)
            deflections_um.append(delta_m * 1e6)
        plt.plot(diameters_mm, deflections_um, label=f"Carbide ({tool_material}), L = {L_mm:.0f} mm", linewidth=2)

    # Comparison curve for HSS at L = 5 mm
    deflections_hss_um = [
        deflection("hss", D_mm * 1e-3, 5e-3, force) * 1e6 for D_mm in diameters_mm
    ]
    plt.plot(diameters_mm, deflections_hss_um, "--", label="HSS, L = 5 mm", linewidth=2, color="gray")

    # 10 µm precision target threshold
    threshold_um = MIN_DEFLECTION * 1e6
    plt.axhline(threshold_um, color="red", linestyle=":", linewidth=2, label="10 µm Target Threshold")

    # Shade the passing region <= 10 µm
    plt.axhspan(0.01, threshold_um, color="green", alpha=0.1, label="Passing Zone (≤ 10 µm)")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(0.1, 5.0)
    plt.ylim(0.01, 1e6)

    plt.xlabel("Tool Diameter D (mm, log scale)", fontsize=11)
    plt.ylabel("Tip Deflection δ (µm, log scale)", fontsize=11)
    plt.title(f"Tool Diameter vs. Tip Deflection (Cutting Force F = {force:.2f} N)", fontsize=12)
    plt.grid(True, which="both", linestyle="--", alpha=0.5)

    plt.legend(loc="upper right", fontsize=10)
    plt.tight_layout()

    if save:
        filename = "tool_diameter_vs_deflection.png"
        plt.savefig(filename, dpi=300)
        print(f"Saved tool diameter vs. deflection plot to {filename}")
    if disp:
        plt.show()
    else:
        plt.close()


def tool_stiffness(tool_material, diameter, length):
    """Cantilever bending stiffness of a cylindrical tool.

    Parameters
    ----------
    tool_material : str
        Key into tool_material_properties.
    diameter : float
        Tool diameter in meters.
    length : float
        Tool overhang / stick-out length in meters.

    Returns
    -------
    float
        Bending stiffness k = 3 E I / L^3 in N/m.
    """
    props = tool_material_properties[tool_material]
    E = props["E"]
    second_moment = math.pi * diameter**4 / 64
    return 3 * E * second_moment / length**3

def deflection(tool_material, diameter, length, load):
    """Tip deflection of a cantilever cylindrical tool under a point load."""
    return load / tool_stiffness(tool_material, diameter, length)

def force_from_deflection(tool_material, diameter, length, deflection):
    """Force that produces a given tip deflection in a cantilever tool.

    This is the inverse of deflection(): F = k * delta.
    """
    return deflection * tool_stiffness(tool_material, diameter, length)

def chip_thickness(feed_per_tooth, rotation_angle):
    """Instantaneous uncut chip thickness for a milling tooth (m).

    Parameters
    ----------
    feed_per_tooth : float
        Feed per tooth in meters.
    rotation_angle : float
        Tooth rotation angle phi in radians; zero is at the top of the tool.

    Returns
    -------
    float
        Chip thickness h = f_z * sin(phi). Zero when the tooth is not engaged.
    """
    return max(0.0, feed_per_tooth * math.sin(rotation_angle))

def tooth_forces(material, feed_per_tooth, axial_depth, rotation_angle):
    """Tangential, radial, and axial forces on one milling tooth (N).

    Uses the mechanistic model:
        F_t = K_tc * h * b + K_te * b
        F_r = K_rc * h * b + K_re * b
        F_a = K_ac * h * b + K_ae * b
    """
    coeffs = cutting_coefficients[material]
    h = chip_thickness(feed_per_tooth, rotation_angle)
    b = axial_depth

    F_t = coeffs["K_tc"] * h * b + coeffs["K_te"] * b
    F_r = coeffs["K_rc"] * h * b + coeffs["K_re"] * b
    F_a = coeffs["K_ac"] * h * b + coeffs["K_ae"] * b
    return F_t, F_r, F_a

def resultant_cutting_force(material, feed_per_tooth, axial_depth, rotation_angle):
    """Magnitude of the total cutting force on one tooth (N)."""
    F_t, F_r, F_a = tooth_forces(material, feed_per_tooth, axial_depth, rotation_angle)
    return math.sqrt(F_t**2 + F_r**2 + F_a**2)


def sdof_tool_parameters(tool_material, diameter, length, damping_ratio=DEFAULT_DAMPING_RATIO):
    """Return consistent SDOF parameters for the tool as a cantilever.

    The tool is modeled as a single-degree-of-freedom oscillator:

        m_eq * x'' + c * x' + k * x = F(t)

    where the effective parameters are derived as follows:

    Bending stiffness (cantilever tip load):
        k = 3 E I / L^3
    with second moment of area:
        I = π D^4 / 64

    Total cylinder mass:
        m_total = ρ * (π D^2 / 4) * L

    Effective tip mass for the first bending mode of a uniform cantilever
    (Rayleigh method using the static deflection shape):
        m_eq ≈ 0.23 * m_total
    This gives the same first natural frequency as the exact Euler-Bernoulli
    solution to within a few percent.

    Natural circular frequency:
        ω_n = sqrt(k / m_eq)   [rad/s]

    Natural frequency in Hz:
        f_n = ω_n / (2 π)

    Viscous damping coefficient for a given damping ratio ζ:
        c = 2 ζ sqrt(k * m_eq) = 2 ζ m_eq ω_n

    Returns
    -------
    dict
        k, m_eq, c, omega_n, f_n, damping_ratio
    """
    k = tool_stiffness(tool_material, diameter, length)
    rho = tool_material_properties[tool_material]["rho"]
    area = math.pi * (diameter / 2) ** 2
    m_total = rho * area * length
    m_eq = 0.23 * m_total
    omega_n = math.sqrt(k / m_eq)
    f_n = omega_n / (2 * math.pi)
    c = 2 * damping_ratio * m_eq * omega_n
    return {
        "k": k,
        "m_total": m_total,
        "m_eq": m_eq,
        "c": c,
        "omega_n": omega_n,
        "f_n": f_n,
        "damping_ratio": damping_ratio,
    }


def tool_natural_frequency_exact(tool_material, diameter, length):
    """Exact first bending natural frequency of a uniform cantilever (Hz).

    Euler-Bernoulli beam theory gives the natural frequencies of a cantilever as

        f_i = (β_i^2 / 2π) * sqrt(E I / (ρ A L^4))

    For the first mode, β_1 = 1.875104. Substituting I = π D^4 / 64 and
    A = π D^2 / 4:

        f_1 = (β_1^2 / 2π) * sqrt((E π D^4 / 64) / (ρ (π D^2 / 4) L^4))
            = (β_1^2 / 2π) * sqrt(E D^2 / (16 ρ L^4))

    This is the benchmark against which the lumped-mass SDOF estimate is checked.
    """
    beta1 = 1.875104
    E = tool_material_properties[tool_material]["E"]
    rho = tool_material_properties[tool_material]["rho"]
    return (beta1**2 / (2 * math.pi)) * math.sqrt(
        E * diameter**2 / (16 * rho * length**4)
    )


def tool_natural_frequency(tool_material, diameter, length):
    """Estimate the first bending natural frequency of a cantilever tool (Hz).

    Uses the lumped-mass SDOF model:
        f_n = (1 / 2π) * sqrt(k / m_eq)
    where k = 3 E I / L^3 and m_eq = 0.23 * m_total.

    This is equivalent to the exact Euler-Bernoulli result to within a few
    percent for a uniform cylinder.
    """
    return sdof_tool_parameters(tool_material, diameter, length)["f_n"]


# ---------------------------------------------------------------------------
# machine base / frame structural model
# ---------------------------------------------------------------------------
# The base is treated as a simple rectangular cantilever beam representing a
# machine column or bridge. Bending is about the weaker axis (height h).


def base_second_moment(width, height):
    """Second moment of area for a rectangular beam, bending about its base axis.

        I = w * h^3 / 12
    """
    return width * height**3 / 12.0


def base_stiffness(base_material, width=DEFAULT_BASE_WIDTH,
                   height=DEFAULT_BASE_HEIGHT, length=DEFAULT_BASE_LENGTH):
    """Cantilever bending stiffness of a rectangular machine base (N/m).

        k_base = 3 * E * I / L^3

    with I = w * h^3 / 12.
    """
    E = base_material_properties[base_material]["E"]
    I_base = base_second_moment(width, height)
    return 3.0 * E * I_base / length**3


def base_deflection(base_material, force, width=DEFAULT_BASE_WIDTH,
                    height=DEFAULT_BASE_HEIGHT, length=DEFAULT_BASE_LENGTH):
    """Static tip deflection of the cantilever base under a point load (m)."""
    return force / base_stiffness(base_material, width, height, length)


def base_sdof_parameters(base_material, width=DEFAULT_BASE_WIDTH,
                         height=DEFAULT_BASE_HEIGHT, length=DEFAULT_BASE_LENGTH,
                         damping_ratio=None):
    """SDOF parameters for the machine base treated as a cantilever beam.

    Uses the full beam mass as the effective moving mass (conservative for a
    first-mode SDOF approximation):

        m_base = ρ * w * h * L
        k_base = 3 E I / L^3
        f_n = (1 / 2π) * sqrt(k_base / m_base)
        c_base = 2 ζ * sqrt(k_base * m_base)

    If damping_ratio is None, the material-specific value from
    base_material_properties is used.
    """
    if damping_ratio is None:
        damping_ratio = base_material_properties[base_material].get(
            "damping_ratio", DEFAULT_BASE_DAMPING_RATIO
        )
    rho = base_material_properties[base_material]["rho"]
    k = base_stiffness(base_material, width, height, length)
    m = rho * width * height * length
    omega_n = math.sqrt(k / m)
    f_n = omega_n / (2 * math.pi)
    c = 2 * damping_ratio * math.sqrt(k * m)
    return {
        "k": k,
        "m": m,
        "c": c,
        "omega_n": omega_n,
        "f_n": f_n,
        "damping_ratio": damping_ratio,
    }


def base_figure_of_merit(base_material, width=DEFAULT_BASE_WIDTH,
                         height=DEFAULT_BASE_HEIGHT, length=DEFAULT_BASE_LENGTH,
                         weights=None):
    """Return a normalized figure of merit (0-1 scale) for a base material.

    The score rewards high stiffness, high damping, and low cost:

        score = w_k * (k / k_max) + w_zeta * (ζ / ζ_max) + w_cost * (c_min / c)

    where k_max and ζ_max are taken over all base materials in the database,
    and c_min is the minimum cost per kg. Default weights:
        w_k = 0.5, w_zeta = 0.3, w_cost = 0.2

    A higher score means a better-balanced base material.
    """
    if weights is None:
        weights = {"k": 0.5, "damping": 0.3, "cost": 0.2}

    # gather database extrema
    k_vals = [base_stiffness(m, width, height, length) for m in base_material_properties]
    zeta_vals = [base_material_properties[m].get("damping_ratio", DEFAULT_BASE_DAMPING_RATIO)
                 for m in base_material_properties]
    cost_vals = [base_material_properties[m]["cost per kg"] for m in base_material_properties]

    k_max = max(k_vals)
    zeta_max = max(zeta_vals)
    cost_min = min(cost_vals)

    k = base_stiffness(base_material, width, height, length)
    zeta = base_material_properties[base_material].get(
        "damping_ratio", DEFAULT_BASE_DAMPING_RATIO
    )
    cost = base_material_properties[base_material]["cost per kg"]

    score = (
        weights["k"] * (k / k_max)
        + weights["damping"] * (zeta / zeta_max)
        + weights["cost"] * (cost_min / cost)
    )
    return score


def base_natural_frequency(base_material, width=DEFAULT_BASE_WIDTH,
                           height=DEFAULT_BASE_HEIGHT, length=DEFAULT_BASE_LENGTH):
    """First bending natural frequency of the machine base (Hz)."""
    return base_sdof_parameters(base_material, width, height, length)["f_n"]


def total_static_deflection(tool_material, base_material, force,
                            diameter, length,
                            base_width=DEFAULT_BASE_WIDTH,
                            base_height=DEFAULT_BASE_HEIGHT,
                            base_length=DEFAULT_BASE_LENGTH):
    """Sum of tool and base static deflections under a point load.

    The tool and base are treated as two springs in series:

        δ_total = δ_tool + δ_base = F/k_tool + F/k_base
    """
    delta_tool = deflection(tool_material, diameter, length, force)
    delta_base = base_deflection(base_material, force, base_width, base_height, base_length)
    return delta_tool + delta_base


# ---------------------------------------------------------------------------
# excitation / frequency helpers
# ---------------------------------------------------------------------------


def spindle_frequency(rpm):
    """Spindle rotation frequency in Hz.

        f_spindle = rpm / 60
    """
    return rpm / 60.0


def tooth_passing_frequency(rpm, n_flutes):
    """Tooth passing frequency in Hz.

        f_tooth = n_flutes * f_spindle = n_flutes * rpm / 60
    """
    return spindle_frequency(rpm) * n_flutes


def frequency_ratio(excitation_freq, natural_freq):
    """Frequency ratio r = f_excitation / f_natural."""
    return excitation_freq / natural_freq


def dynamic_amplification(r, damping_ratio):
    """Dynamic amplification factor of a single-DOF system.

    For harmonic forcing at frequency ratio r = f/f_n, the steady-state
    displacement amplification relative to the static deflection is

        H(r) = 1 / sqrt((1 - r^2)^2 + (2 ζ r)^2)

    At r = 1 (resonance), H ≈ 1 / (2 ζ).
    """
    return 1.0 / math.sqrt((1.0 - r**2) ** 2 + (2.0 * damping_ratio * r) ** 2)


def sdof_frequency_response(freqs, f_n, damping_ratio, k):
    """Compute the SDOF displacement FRF: X/F as a function of frequency.

    The complex frequency response function is

        H(ω) = 1 / (k - m ω^2 + i c ω)

    with magnitude

        |H(ω)| = 1 / (k * sqrt((1 - r^2)^2 + (2 ζ r)^2))

    where r = f/f_n and ω = 2π f.
    """
    import numpy as np
    r = np.array(freqs) / f_n
    mag = 1.0 / (k * np.sqrt((1.0 - r**2) ** 2 + (2.0 * damping_ratio * r) ** 2))
    return mag


def multi_flute_forcing(material, feed_per_tooth, axial_depth, rpm, n_flutes, t,
                        runout_amplitude=0.0, runout_phase=0.0):
    """Compute tangential, radial, and axial forces for a multi-flute end mill.

    Each tooth j enters the cut at angle

        phi_j(t) = omega * t + j * 2π / N_f

    where omega = 2π * f_spindle = 2π * rpm / 60. The uncut chip thickness for
    tooth j is

        h_j(t) = max(0, f_z * sin(phi_j(t)) + runout_amplitude * sin(phi_j(t) + runout_phase))

    The total force is the sum over all teeth currently engaged.

    For simplicity this function assumes full radial immersion (slotting), so a
    tooth is engaged whenever sin(phi_j) > 0, i.e. 0 < phi_j < π.

    Returns the total tangential, radial, and axial forces at time t.
    """
    coeffs = cutting_coefficients[material]
    omega = 2 * math.pi * rpm / 60.0
    b = axial_depth

    F_t_total = 0.0
    F_r_total = 0.0
    F_a_total = 0.0
    for j in range(n_flutes):
        phi_j = omega * t + j * 2 * math.pi / n_flutes
        # Chip thickness from feed plus optional runout displacement
        h_j = max(0.0, feed_per_tooth * math.sin(phi_j) + runout_amplitude * math.sin(phi_j + runout_phase))
        F_t_total += coeffs["K_tc"] * h_j * b + coeffs["K_te"] * b
        F_r_total += coeffs["K_rc"] * h_j * b + coeffs["K_re"] * b
        F_a_total += coeffs["K_ac"] * h_j * b + coeffs["K_ae"] * b

    return F_t_total, F_r_total, F_a_total


def vibration_analysis(results, n_flutes=2):
    """Frequency-domain analysis of simulated cutting force and deflection.

    Parameters
    ----------
    results : dict
        Output from run_python_simulation.
    n_flutes : int
        Number of tool flutes.

    Returns
    -------
    dict
        Frequencies, FFT spectra, SDOF frequency-response function, and key
        markers (spindle frequency, tooth passing frequency, natural frequency,
        frequency ratios, dynamic amplification factors).
    """
    import numpy as np
    
    t = results["time"]
    dt = float(np.mean(np.diff(t)))
    fs = 1.0 / dt
    n = len(t)

    freqs = np.fft.rfftfreq(n, dt)

    # Detrend (remove mean) and window to avoid spectral leakage.
    # A Hanning window tapers the signal to zero at the ends, reducing side
    # lobes in the FFT.
    force_detrended = results["F_res"] - np.mean(results["F_res"])
    deflection_detrended = results["deflection"] - np.mean(results["deflection"])
    window = np.hanning(n)

    # Windowed FFT; normalize by window area so amplitude is comparable.
    force_fft = np.abs(np.fft.rfft(force_detrended * window))
    deflection_fft = np.abs(np.fft.rfft(deflection_detrended * window))
    force_fft *= 2.0 / np.sum(window)
    deflection_fft *= 2.0 / np.sum(window)

    rpm = results["rpm"]
    f_spindle = spindle_frequency(rpm)
    f_tooth = tooth_passing_frequency(rpm, n_flutes)
    tool_material = results["tool_material"]
    diameter = results["diameter"]
    length = results["length"]

    sdof = sdof_tool_parameters(tool_material, diameter, length)
    f_natural = sdof["f_n"]
    zeta = sdof["damping_ratio"]
    k = sdof["k"]
    f_natural_exact = tool_natural_frequency_exact(tool_material, diameter, length)

    # SDOF analytical frequency response for the same frequency grid
    frf = sdof_frequency_response(freqs, f_natural, zeta, k)

    # Dynamic amplification at the forcing frequencies
    r_spindle = frequency_ratio(f_spindle, f_natural)
    r_tooth = frequency_ratio(f_tooth, f_natural)
    amp_spindle = dynamic_amplification(r_spindle, zeta)
    amp_tooth = dynamic_amplification(r_tooth, zeta)

    # Dominant peaks, excluding DC (index 0)
    force_peak_idx = np.argmax(force_fft[1:]) + 1
    deflection_peak_idx = np.argmax(deflection_fft[1:]) + 1

    return {
        "fs": fs,
        "freqs": freqs,
        "force_fft": force_fft,
        "deflection_fft": deflection_fft,
        "frf": frf,
        "f_spindle": f_spindle,
        "f_tooth": f_tooth,
        "f_natural": f_natural,
        "f_natural_exact": f_natural_exact,
        "r_spindle": r_spindle,
        "r_tooth": r_tooth,
        "amp_spindle": amp_spindle,
        "amp_tooth": amp_tooth,
        "damping_ratio": zeta,
        "force_peak_freq": float(freqs[force_peak_idx]),
        "force_peak_amp": float(force_fft[force_peak_idx]),
        "deflection_peak_freq": float(freqs[deflection_peak_idx]),
        "deflection_peak_amp": float(deflection_fft[deflection_peak_idx]),
        "n_flutes": n_flutes,
    }


def _near_any_harmonic(f_dom, base_freq, n_harmonics=5, tolerance=0.05):
    """Check whether f_dom is within tolerance of any harmonic of base_freq."""
    for n in range(1, n_harmonics + 1):
        f_harm = base_freq * n
        if abs(f_dom - f_harm) / f_harm < tolerance:
            return True
    return False


def chatter_indicator(results, n_flutes=2, tolerance=0.05, n_harmonics=5):
    """Improved chatter indicator based on frequency-domain peak matching.

    A regenerative chatter vibration typically appears at (or near) the tool's
    structural natural frequency and is sustained by the cutting process, not
    directly locked to spindle or tooth-passing harmonics.

    Logic:
        1. Find the dominant frequency in the deflection spectrum.
        2. If it is near f_natural and NOT near any spindle harmonic or tooth-
           passing harmonic up to n_harmonics, flag potential chatter.

    This remains a coarse indicator because true chatter requires a time-delay
    (regenerative) model and stability lobes.
    """
    freqs = vibration_analysis(results, n_flutes)
    f_dom = freqs["deflection_peak_freq"]
    f_natural = freqs["f_natural"]
    f_spindle = freqs["f_spindle"]
    f_tooth = freqs["f_tooth"]

    near_natural = abs(f_dom - f_natural) / f_natural < tolerance
    near_spindle_harmonics = _near_any_harmonic(f_dom, f_spindle, n_harmonics, tolerance)
    near_tooth_harmonics = _near_any_harmonic(f_dom, f_tooth, n_harmonics, tolerance)

    return near_natural and not (near_spindle_harmonics or near_tooth_harmonics)


def plot_vibration_analysis(vib, save=False):
    """Plot FFT spectra with analytical FRF and frequency markers."""
    import numpy as np
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

    freqs = vib["freqs"]
    f_max = max(3 * vib["f_tooth"], 2 * vib["f_natural"], 5000)
    f_max = min(f_max, vib["fs"] / 2)
    mask = freqs <= f_max

    axes[0].plot(freqs[mask], vib["force_fft"][mask], label="Force spectrum")
    axes[0].axvline(vib["f_spindle"], color="g", linestyle="--", label="spindle")
    axes[0].axvline(vib["f_tooth"], color="r", linestyle="--", label="tooth pass")
    axes[0].axvline(vib["f_natural"], color="m", linestyle="--", label="natural (SDOF)")
    axes[0].axvline(vib["f_natural_exact"], color="m", linestyle=":", label="natural (exact)")
    axes[0].set_ylabel("Force FFT amplitude")
    axes[0].set_title("Cutting force frequency spectrum")
    axes[0].legend()

    axes[1].plot(freqs[mask], vib["deflection_fft"][mask], label="Deflection spectrum")
    axes[1].axvline(vib["f_spindle"], color="g", linestyle="--", label="spindle")
    axes[1].axvline(vib["f_tooth"], color="r", linestyle="--", label="tooth pass")
    axes[1].axvline(vib["f_natural"], color="m", linestyle="--", label="natural (SDOF)")
    axes[1].axvline(vib["f_natural_exact"], color="m", linestyle=":", label="natural (exact)")
    axes[1].set_ylabel("Deflection FFT amplitude")
    axes[1].set_title("Tool deflection frequency spectrum")
    axes[1].legend()

    # Plot analytical FRF normalized to its peak for visibility
    frf_norm = vib["frf"] / np.max(vib["frf"])
    axes[2].plot(freqs[mask], frf_norm[mask], label="|X/F| normalized")
    axes[2].axvline(vib["f_spindle"], color="g", linestyle="--", label="spindle")
    axes[2].axvline(vib["f_tooth"], color="r", linestyle="--", label="tooth pass")
    axes[2].axvline(vib["f_natural"], color="m", linestyle="--", label="natural")
    axes[2].set_xlabel("Frequency (Hz)")
    axes[2].set_ylabel("Normalized |H(f)|")
    axes[2].set_title("SDOF frequency response function (analytical)")
    axes[2].legend()

    plt.tight_layout()
    if save:
        plt.savefig("vibration_analysis.png")
    else:
        plt.show()


def print_vibration_summary(vib):
    """Print a concise vibration analysis summary."""
    print("\n--- Vibration analysis (refined SDOF model) ---")
    print(f"  Sampling rate: {vib['fs']:.1f} Hz")
    print(f"  Spindle frequency: {vib['f_spindle']:.1f} Hz")
    print(f"  Tooth-passing frequency ({vib['n_flutes']} flutes): {vib['f_tooth']:.1f} Hz")
    print(f"  Estimated tool natural frequency (SDOF): {vib['f_natural']:.1f} Hz")
    print(f"  Exact cantilever natural frequency: {vib['f_natural_exact']:.1f} Hz")
    print(f"  Frequency ratio r_spindle = f_spindle/f_n: {vib['r_spindle']:.3f}")
    print(f"  Frequency ratio r_tooth = f_tooth/f_n: {vib['r_tooth']:.3f}")
    print(f"  Dynamic amplification at spindle freq: {vib['amp_spindle']:.3f}")
    print(f"  Dynamic amplification at tooth-pass freq: {vib['amp_tooth']:.3f}")
    print(f"  Dominant force frequency: {vib['force_peak_freq']:.1f} Hz")
    print(f"  Dominant deflection frequency: {vib['deflection_peak_freq']:.1f} Hz")


def run_vibration_comparison(diameter=0.001, length=0.010, rpm=None, feed_per_tooth=None,
                             axial_depth=None, stop_time=None, n_flutes=2,
                             tool_material=DEFAULT_TOOL_MATERIAL):
    """Run vibration analysis across all workpiece materials with cutting data.

    Tool material is held fixed.
    """
    rpm = rpm if rpm is not None else DEFAULT_SPINDLE_RPM
    feed_per_tooth = feed_per_tooth if feed_per_tooth is not None else DEFAULT_FEED_PER_TOOTH
    axial_depth = axial_depth if axial_depth is not None else DEFAULT_AXIAL_DEPTH
    stop_time = stop_time if stop_time is not None else DEFAULT_STOP_TIME

    materials = [m for m in material_properties if m in cutting_coefficients]
    rows = []
    for material in materials:
        results = run_python_simulation(
            material=material,
            diameter=diameter,
            length=length,
            rpm=rpm,
            feed_per_tooth=feed_per_tooth,
            axial_depth=axial_depth,
            stop_time=stop_time,
            n_flutes=n_flutes,
            tool_material=tool_material,
        )
        vib = vibration_analysis(results, n_flutes)
        rows.append({
            "material": material,
            "f_spindle": vib["f_spindle"],
            "f_tooth": vib["f_tooth"],
            "f_natural": vib["f_natural"],
            "f_natural_exact": vib["f_natural_exact"],
            "r_spindle": vib["r_spindle"],
            "r_tooth": vib["r_tooth"],
            "amp_spindle": vib["amp_spindle"],
            "amp_tooth": vib["amp_tooth"],
            "dominant_deflection_freq": vib["deflection_peak_freq"],
            "chatter_warning": chatter_indicator(results, n_flutes),
        })
    return rows


def print_vibration_comparison(rows):
    """Print vibration comparison results as a formatted markdown table."""
    print("\n--- Vibration comparison across materials ---")
    print("| Material | f_n SDOF (Hz) | f_n exact (Hz) | r_spindle | r_tooth | amp_spindle | amp_tooth | Chatter warning |")
    print("|---------|---------------|----------------|-----------|---------|-------------|-----------|-----------------|")
    for r in rows:
        print(
            f"| {r['material']:<8} "
            f"| {r['f_natural']:>13.1f} "
            f"| {r['f_natural_exact']:>14.1f} "
            f"| {r['r_spindle']:>9.3f} "
            f"| {r['r_tooth']:>7.3f} "
            f"| {r['amp_spindle']:>11.3f} "
            f"| {r['amp_tooth']:>9.3f} "
            f"| {str(r['chatter_warning']):<15} |"
        )


# ---------------------------------------------------------------------------
# Simulink bridge
# ---------------------------------------------------------------------------

def _ensure_simulink_model(engine, model_name=SIMULINK_MODEL_NAME):
    """Create the reference Simulink model if it does not already exist."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    simulink_dir = os.path.join(script_dir, "simulink")
    builder_path = os.path.join(simulink_dir, "build_cnc_simulink_model.m")
    runner_path = os.path.join(simulink_dir, "run_cnc_simulink_sim.m")

    if not os.path.exists(builder_path):
        raise FileNotFoundError(f"MATLAB model builder not found: {builder_path}")
    if not os.path.exists(runner_path):
        raise FileNotFoundError(f"MATLAB runner not found: {runner_path}")

    # Add project folder to MATLAB path so the .m files are discoverable
    engine.addpath(script_dir, nargout=0)
    engine.addpath(simulink_dir, nargout=0)

    # Check that Simulink is licensed before trying to create/load models.
    simulink_ok = engine.license("test", "Simulink", nargout=1)
    if not simulink_ok:
        raise RuntimeError(
            "Simulink license is not available in this MATLAB session. "
            "Make sure your MATLAB installation includes a valid Simulink license."
        )

    model_file = os.path.join(script_dir, f"{model_name}.slx")
    try:
        if not os.path.exists(model_file):
            print(f"Building Simulink model '{model_name}'...")
            engine.build_cnc_simulink_model(model_name, nargout=0)
        else:
            engine.load_system(model_name, nargout=0)
    except Exception as e:
        msg = str(e)
        if "requires simulink" in msg.lower() or "undefined function 'new_system'" in msg.lower():
            raise RuntimeError(
                "Simulink is licensed, but its toolbox files are not on the "
                "MATLAB Engine path. This usually happens on macOS when Simulink "
                "works in the interactive desktop but not in a fresh engine session. "
                "Workarounds: (1) open MATLAB desktop, run 'simulink', then use "
                "matlab.engine.connect_matlab() to connect to it; (2) run the model "
                "directly inside MATLAB; (3) use the Python-native simulation in sims.py."
            ) from e
        raise


def thermal_params(tool_material, diameter, length, heat_partition=HEAT_PARTITION):
    """Return a rough lumped thermal model for the tool.

    Model assumptions:
        - Heat is generated at the tool tip: Q = heat_partition * F_t * V_c.
        - Heat conducts away through the full tool length to the collet:
          R_th = L / (k * A) with A = pi * D^2 / 4.
        - Thermal capacitance is based on the full tool volume:
          C_th = rho * c_p * V.
    """
    props = tool_material_properties[tool_material]
    area = math.pi * diameter**2 / 4
    volume = area * length
    R_th = length / (props["k"] * area)
    C_th = props["rho"] * props["c_p"] * volume
    return {
        "thermal_resistance": float(R_th),
        "thermal_capacitance": float(C_th),
        "heat_partition": float(heat_partition),
    }


def _simulink_workspace_params(material, diameter, length, rpm=None, feed_per_tooth=None,
                                axial_depth=None, stop_time=None,
                                tool_material=DEFAULT_TOOL_MATERIAL):
    """Build a flat dict of parameters to push into the MATLAB base workspace.

    Parameters
    ----------
    material : str
        Workpiece material (key into cutting_coefficients).
    tool_material : str
        Tool material (key into tool_material_properties).
    """
    rpm = rpm if rpm is not None else DEFAULT_SPINDLE_RPM
    feed_per_tooth = feed_per_tooth if feed_per_tooth is not None else DEFAULT_FEED_PER_TOOTH
    axial_depth = axial_depth if axial_depth is not None else DEFAULT_AXIAL_DEPTH
    stop_time = stop_time if stop_time is not None else DEFAULT_STOP_TIME

    coeffs = cutting_coefficients[material]

    # Consistent SDOF tool parameters (same as the Python-native model)
    sdof = sdof_tool_parameters(tool_material, diameter, length)

    thermal = thermal_params(tool_material, diameter, length)

    params = {
        "rpm": float(rpm),
        "feed_per_tooth": float(feed_per_tooth),
        "axial_depth": float(axial_depth),
        "tool_diameter": float(diameter),
        "tool_length": float(length),
        "tool_stiffness": sdof["k"],
        "tool_mass": sdof["m_eq"],
        "tool_damping": sdof["c"],
        "sim_stop_time": float(stop_time),
        "K_tc": float(coeffs["K_tc"]),
        "K_rc": float(coeffs["K_rc"]),
        "K_ac": float(coeffs["K_ac"]),
        "K_te": float(coeffs["K_te"]),
        "K_re": float(coeffs["K_re"]),
        "K_ae": float(coeffs["K_ae"]),
    }
    params.update(thermal)
    return params


def run_simulink_simulation(material, diameter, length, rpm=None, feed_per_tooth=None,
                            axial_depth=None, stop_time=None,
                            tool_material=DEFAULT_TOOL_MATERIAL):
    """Run the reference Simulink CNC model from Python.

    Parameters
    ----------
    material : str
        Workpiece material (key into cutting_coefficients).
    tool_material : str
        Tool material (key into tool_material_properties).
    diameter : float
        Tool diameter in meters.
    length : float
        Tool overhang length in meters.
    rpm : float, optional
        Spindle speed in rev/min.
    feed_per_tooth : float, optional
        Feed per tooth in meters.
    axial_depth : float, optional
        Axial depth of cut in meters.
    stop_time : float, optional
        Simulation stop time in seconds.

    Returns
    -------
    dict
        Simulation results extracted from Simulink scopes.

    Notes
    -----
    This requires a live MATLAB + Simulink installation and the
    `matlabengine` Python package matching the MATLAB version.
    """
    if _MATLAB_IMPORT_ERROR is not None:
        raise RuntimeError(
            f"matlab.engine could not be imported ({_MATLAB_IMPORT_ERROR}). "
            "Run 'python setup_matlab_bridge.py' and install the matching "
            "matlabengine version."
        )

    print("Starting MATLAB engine...")
    eng = matlab.engine.start_matlab()
    try:
        params = _simulink_workspace_params(
            material, diameter, length, rpm, feed_per_tooth, axial_depth, stop_time,
            tool_material
        )

        # Push all parameters into the MATLAB base workspace
        for key, value in params.items():
            eng.workspace[key] = float(value)

        _ensure_simulink_model(eng, SIMULINK_MODEL_NAME)

        print(f"Running Simulink simulation for {material}...")
        results = eng.run_cnc_simulink_sim(SIMULINK_MODEL_NAME)
        return results
    finally:
        eng.quit()


def simulink_summary(material, diameter, length, rpm=None, feed_per_tooth=None,
                     axial_depth=None, stop_time=None,
                     tool_material=DEFAULT_TOOL_MATERIAL):
    """Convenience wrapper: run the Simulink model and print a short summary."""
    results = run_simulink_simulation(
        material, diameter, length, rpm, feed_per_tooth, axial_depth, stop_time,
        tool_material
    )

    summary = {}
    for name, signal in results.items():
        if signal["data"].size:
            summary[name] = float(signal["data"][-1])
        else:
            summary[name] = None

    return summary


def run_python_simulation(material, diameter, length, rpm=None, feed_per_tooth=None,
                          axial_depth=None, stop_time=None, num_points=1000,
                          n_flutes=1, runout_amplitude=0.0, runout_phase=0.0,
                          tool_material=DEFAULT_TOOL_MATERIAL):
    """Run a Python-native time-domain CNC simulation.

    This mirrors the reference Simulink model but uses scipy.integrate.solve_ivp
    instead of Simulink, so it works without MATLAB installed.

    Parameters
    ----------
    material : str
        Workpiece material (key into cutting_coefficients).
    tool_material : str
        Tool material (key into tool_material_properties).
    diameter : float
        Tool diameter in meters.
    length : float
        Tool overhang length in meters.

    States:
        x[0] = tool tip deflection (m)
        x[1] = tool tip velocity (m/s)
        x[2] = tool tip temperature rise above ambient (K)

    The tool is modeled as a single-degree-of-freedom oscillator:

        m_eq * x'' + c * x' + k * x = F_t(t)

    where:
        k      = 3 E I / L^3          (cantilever bending stiffness)
        m_eq   = 0.23 * m_total       (effective tip mass, first bending mode)
        c      = 2 ζ * sqrt(k * m_eq) (viscous damping, ζ default 0.05)
        F_t(t) = Σ_j [K_tc * h_j(t) + K_te] * a_p   (multi-flute forcing)

    For one tooth the chip thickness is:
        h_j(t) = max(0, f_z * sin(φ_j(t)) + runout_j(t))
    with tooth j angle:
        φ_j(t) = ω * t + j · 2π / N_f
    and spindle angular speed:
        ω = 2π * rpm / 60

    Thermal model:
        C_th * dT/dt = η * F_t(t) * V_c - T / R_th
    where η = HEAT_PARTITION, V_c = π D * rpm / 60.

    Returns
    -------
    dict
        time, deflection, velocity, temperature, force components,
        resultant force arrays, plus SDOF and thermal parameters.
    """
    rpm = rpm if rpm is not None else DEFAULT_SPINDLE_RPM
    feed_per_tooth = feed_per_tooth if feed_per_tooth is not None else DEFAULT_FEED_PER_TOOTH
    axial_depth = axial_depth if axial_depth is not None else DEFAULT_AXIAL_DEPTH
    stop_time = stop_time if stop_time is not None else DEFAULT_STOP_TIME

    # Consistent SDOF tool dynamics parameters
    sdof = sdof_tool_parameters(tool_material, diameter, length)
    k = sdof["k"]
    m_eq = sdof["m_eq"]
    c = sdof["c"]

    # Thermal parameters (lumped full-tool conduction model)
    thermal = thermal_params(tool_material, diameter, length)
    R_th = thermal["thermal_resistance"]
    C_th = thermal["thermal_capacitance"]
    heat_partition = thermal["heat_partition"]

    V_c = math.pi * diameter * rpm / 60.0  # m/s cutting speed

    def forcing(t):
        # Sum forces from all flutes; runout adds a sinusoidal perturbation to
        # the chip thickness of each tooth.
        return multi_flute_forcing(
            material, feed_per_tooth, axial_depth, rpm, n_flutes, t,
            runout_amplitude, runout_phase
        )

    def rhs(t, x):
        deflection, velocity, T = x
        F_t, F_r, F_a = forcing(t)

        # Mechanical: m_eq * x'' + c * x' + k * x = F_t
        accel = (F_t - c * velocity - k * deflection) / m_eq

        # Thermal: C_th * dT/dt = η * F_t * V_c - T / R_th
        heat_in = heat_partition * F_t * V_c
        dTdt = (heat_in - T / R_th) / C_th

        return [velocity, accel, dTdt]

    sol = solve_ivp(rhs, [0, stop_time], [0, 0, 0], t_eval=np.linspace(0, stop_time, num_points))

    F_t_arr = np.zeros_like(sol.t)
    F_r_arr = np.zeros_like(sol.t)
    F_a_arr = np.zeros_like(sol.t)
    for i, t in enumerate(sol.t):
        F_t_arr[i], F_r_arr[i], F_a_arr[i] = forcing(t)
    F_res_arr = np.sqrt(F_t_arr**2 + F_r_arr**2 + F_a_arr**2)

    return {
        "time": sol.t,
        "deflection": sol.y[0],
        "velocity": sol.y[1],
        "temperature": sol.y[2],
        "F_t": F_t_arr,
        "F_r": F_r_arr,
        "F_a": F_a_arr,
        "F_res": F_res_arr,
        "material": material,
        "tool_material": tool_material,
        "diameter": diameter,
        "length": length,
        "rpm": rpm,
        "feed_per_tooth": feed_per_tooth,
        "axial_depth": axial_depth,
        "n_flutes": n_flutes,
        "sdof": sdof,
        "thermal": thermal,
    }


def plot_simulation_results(results, save=False):
    """Plot the time-domain simulation results."""
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    t = results["time"]
    axes[0].plot(t, results["F_t"], label="F_t")
    axes[0].plot(t, results["F_r"], label="F_r")
    axes[0].plot(t, results["F_a"], label="F_a")
    axes[0].plot(t, results["F_res"], label="F_res", color="k", linestyle="--")
    axes[0].set_ylabel("Force (N)")
    axes[0].set_title(f"Cutting forces ({results['material']}, D={results['diameter']*1e3:.1f} mm)")
    axes[0].legend()

    axes[1].plot(t, results["deflection"] * 1e6)
    axes[1].axhline(MIN_DEFLECTION * 1e6, color="r", linestyle="--", label="min deflection")
    axes[1].set_ylabel("Deflection (µm)")
    axes[1].set_title("Tool tip deflection")
    axes[1].legend()

    axes[2].plot(t, results["velocity"] * 1e3)
    axes[2].set_ylabel("Velocity (mm/s)")
    axes[2].set_title("Tool tip velocity")

    axes[3].plot(t, results["temperature"])
    axes[3].set_ylabel("Temperature rise (K)")
    axes[3].set_title("Tool tip temperature rise")
    axes[3].set_xlabel("Time (s)")

    plt.tight_layout()
    if save:
        plt.savefig("simulation_results.png")
    plt.show()


def run_comparison(diameter=0.001, length=0.010, rpm=None, feed_per_tooth=None,
                     axial_depth=None, stop_time=None, n_flutes=2,
                     tool_material=DEFAULT_TOOL_MATERIAL):
    """Run the Python-native simulation across all workpiece materials with cutting data.

    Tool material is held fixed. Returns a list of dicts with key metrics per
    workpiece material.
    """
    rpm = rpm if rpm is not None else DEFAULT_SPINDLE_RPM
    feed_per_tooth = feed_per_tooth if feed_per_tooth is not None else DEFAULT_FEED_PER_TOOTH
    axial_depth = axial_depth if axial_depth is not None else DEFAULT_AXIAL_DEPTH
    stop_time = stop_time if stop_time is not None else DEFAULT_STOP_TIME

    materials = [m for m in material_properties if m in cutting_coefficients]
    rows = []
    for material in materials:
        stiffness = tool_stiffness(tool_material, diameter, length)
        deflection_10 = deflection(tool_material, diameter, length, MIN_CUTTING_LOAD)
        deflection_300 = deflection(tool_material, diameter, length, MAX_CUTTING_LOAD)
        results = run_python_simulation(
            material=material,
            diameter=diameter,
            length=length,
            rpm=rpm,
            feed_per_tooth=feed_per_tooth,
            axial_depth=axial_depth,
            stop_time=stop_time,
            n_flutes=n_flutes,
            tool_material=tool_material,
        )
        rows.append({
            "material": material,
            "stiffness_N_per_m": stiffness,
            "deflection_10_N_um": deflection_10 * 1e6,
            "deflection_300_N_um": deflection_300 * 1e6,
            "peak_force_N": results["F_res"].max(),
            "peak_deflection_um": results["deflection"].max() * 1e6,
            "peak_temp_rise_K": results["temperature"].max(),
            "natural_freq_Hz": results["sdof"]["f_n"],
        })
    return rows


def print_comparison(rows):
    """Print material comparison results as a formatted markdown table."""
    print("\n--- Material comparison (Python-native reference simulation) ---")
    header = (
        "| Material | Stiffness (N/m) | δ@10N (µm) | δ@300N (µm) | "
        "Peak F (N) | Peak δ (µm) | Peak ΔT (K) | f_n (Hz) |"
    )
    print(header)
    print("|" + "-" * 9 + "|" + "-" * 17 + "|" + "-" * 12 + "|" + "-" * 13 + "|" + "-" * 12 + "|" + "-" * 13 + "|" + "-" * 13 + "|" + "-" * 10 + "|")
    for r in rows:
        print(
            f"| {r['material']:<8} "
            f"| {r['stiffness_N_per_m']:.2e} "
            f"| {r['deflection_10_N_um']:>10.1f} "
            f"| {r['deflection_300_N_um']:>11.1f} "
            f"| {r['peak_force_N']:>10.2f} "
            f"| {r['peak_deflection_um']:>11.1f} "
            f"| {r['peak_temp_rise_K']:>11.2f} "
            f"| {r['natural_freq_Hz']:>8.1f} |"
        )


def run_base_comparison(peak_force, rpm=None, n_flutes=None,
                        tool_material=DEFAULT_TOOL_MATERIAL,
                        diameter=0.001, length=0.010,
                        base_width=DEFAULT_BASE_WIDTH,
                        base_height=DEFAULT_BASE_HEIGHT,
                        base_length=DEFAULT_BASE_LENGTH,
                        damping_ratio=None):
    """Compare machine base/frame materials under a representative cutting load.

    For each base material the model reports:
        * cantilever bending stiffness k_base = 3 E I / L^3
        * static tip deflection under the peak cutting force, δ_static = F / k_base
        * first bending natural frequency, f_n = (1/2π) sqrt(k_base / m_base)
        * dynamic amplification at the tooth-passing frequency,
          H(r) = 1 / sqrt((1 - r^2)^2 + (2 ζ r)^2)
        * approximate dynamic vibration amplitude,
          δ_dynamic = H(r) * δ_static (treating the peak force as a harmonic amplitude)
        * total static deflection budget (tool + base in series)
        * material damping ratio ζ and cost per kg
        * normalized figure of merit combining stiffness, damping, and cost

    The tool contribution is held fixed so the table isolates the base material
    effect. If damping_ratio is None, the material-specific damping ratio from
    base_material_properties is used.
    """
    rpm = rpm if rpm is not None else DEFAULT_SPINDLE_RPM
    n_flutes = n_flutes if n_flutes is not None else DEFAULT_N_FLUTES
    f_tooth = tooth_passing_frequency(rpm, n_flutes)

    rows = []
    for base_material in base_material_properties:
        k_base = base_stiffness(base_material, base_width, base_height, base_length)
        delta_static = base_deflection(base_material, peak_force,
                                       base_width, base_height, base_length)
        params = base_sdof_parameters(base_material, base_width, base_height,
                                        base_length, damping_ratio)
        f_n = params["f_n"]
        zeta = params["damping_ratio"]
        r = frequency_ratio(f_tooth, f_n)
        H = dynamic_amplification(r, zeta)
        delta_dynamic = delta_static * H
        delta_total = total_static_deflection(
            tool_material, base_material, peak_force,
            diameter, length, base_width, base_height, base_length
        )
        cost = base_material_properties[base_material]["cost per kg"]
        fom = base_figure_of_merit(base_material, base_width, base_height, base_length)
        rows.append({
            "base_material": base_material,
            "k_base_N_per_m": k_base,
            "static_delta_um": delta_static * 1e6,
            "dynamic_delta_um": delta_dynamic * 1e6,
            "f_n_Hz": f_n,
            "damping_ratio": zeta,
            "freq_ratio": r,
            "amplification": H,
            "cost_per_kg": cost,
            "figure_of_merit": fom,
            "total_delta_um": delta_total * 1e6,
        })
    return rows


def print_base_comparison(rows, peak_force, rpm=None, n_flutes=None):
    """Print base-material comparison as a formatted markdown table."""
    rpm = rpm if rpm is not None else DEFAULT_SPINDLE_RPM
    n_flutes = n_flutes if n_flutes is not None else DEFAULT_N_FLUTES
    f_tooth = tooth_passing_frequency(rpm, n_flutes)
    print("\n--- Base / frame material comparison ---")
    print(f"(peak transmitted force = {peak_force:.2f} N, tooth-passing freq = {f_tooth:.0f} Hz)")
    header = (
        "| Base material | k (N/m) | δ_static (µm) | δ_dynamic (µm) | "
        "f_n (Hz) | ζ | r | H | cost ($/kg) | FoM | δ_total (µm) |"
    )
    print(header)
    print("|" + "-" * 15 + "|" + "-" * 9 + "|" + "-" * 15 + "|" + "-" * 16 + "|" + "-" * 10 + "|" + "-" * 5 + "|" + "-" * 5 + "|" + "-" * 5 + "|" + "-" * 13 + "|" + "-" * 5 + "|" + "-" * 14 + "|")
    for r in rows:
        verdict = "OK" if r["total_delta_um"] <= MIN_DEFLECTION * 1e6 else "EXCEEDS"
        print(
            f"| {r['base_material']:<13} "
            f"| {r['k_base_N_per_m']:.2e} "
            f"| {r['static_delta_um']:>13.2f} "
            f"| {r['dynamic_delta_um']:>14.2f} "
            f"| {r['f_n_Hz']:>8.1f} "
            f"| {r['damping_ratio']:>3.4f} "
            f"| {r['freq_ratio']:>3.2f} "
            f"| {r['amplification']:>3.2f} "
            f"| {r['cost_per_kg']:>11.2f} "
            f"| {r['figure_of_merit']:>3.2f} "
            f"| {r['total_delta_um']:>12.2f} {verdict} |"
        )


def run_geometry_base_sweep(peak_force,
                          tool_material=DEFAULT_TOOL_MATERIAL,
                          diameters=None, lengths=None,
                          base_width=DEFAULT_BASE_WIDTH,
                          base_height=DEFAULT_BASE_HEIGHT,
                          base_length=DEFAULT_BASE_LENGTH):
    """Sweep tool diameter / stick-out and base material against deflection budget.

    For each combination of (diameter, length, base_material) the total static
    deflection is computed as two springs in series:

        δ_total(F) = F/k_tool(D,L) + F/k_base

    where k_tool = 3 E_tool I_tool / L^3 and k_base = 3 E_base I_base / L_base^3.

    Returns a list of result dictionaries sorted by total deflection (ascending).
    """
    if diameters is None:
        diameters = np.array([0.0005, 0.001, 0.0015, 0.002, 0.003])
    if lengths is None:
        lengths = np.array([0.005, 0.007, 0.010, 0.015])

    rows = []
    for diameter in diameters:
        for length in lengths:
            for base_material in base_material_properties:
                delta_total = total_static_deflection(
                    tool_material, base_material, peak_force,
                    diameter, length, base_width, base_height, base_length
                )
                rows.append({
                    "diameter_mm": diameter * 1e3,
                    "length_mm": length * 1e3,
                    "base_material": base_material,
                    "tool_stiffness_N_per_m": tool_stiffness(tool_material, diameter, length),
                    "base_stiffness_N_per_m": base_stiffness(base_material, base_width, base_height, base_length),
                    "total_delta_um": delta_total * 1e6,
                    "meets_threshold": delta_total <= MIN_DEFLECTION,
                })
    rows.sort(key=lambda r: r["total_delta_um"])
    return rows


def print_geometry_base_sweep(rows, peak_force, top_n=20):
    """Print the geometry/base sweep, showing the lowest-deflection options."""
    print("\n--- Tool geometry + base material sweep ---")
    print(f"(peak transmitted force = {peak_force:.2f} N, threshold = {MIN_DEFLECTION*1e6:.1f} µm)")
    print(f"Showing top {top_n} configurations by total deflection.")
    header = (
        "| D (mm) | L (mm) | base | k_tool (N/m) | k_base (N/m) | "
        "δ_total (µm) | status |"
    )
    print(header)
    print("|" + "-" * 8 + "|" + "-" * 8 + "|" + "-" * 14 + "|" + "-" * 14 + "|" + "-" * 14 + "|" + "-" * 14 + "|" + "-" * 8 + "|")
    for r in rows[:top_n]:
        status = "OK" if r["meets_threshold"] else "EXCEEDS"
        print(
            f"| {r['diameter_mm']:>6.2f} "
            f"| {r['length_mm']:>6.1f} "
            f"| {r['base_material']:<12} "
            f"| {r['tool_stiffness_N_per_m']:.2e} "
            f"| {r['base_stiffness_N_per_m']:.2e} "
            f"| {r['total_delta_um']:>12.2f} "
            f"| {status:<6} |"
        )
    n_ok = sum(1 for r in rows if r["meets_threshold"])
    print(f"\nCombinations meeting {MIN_DEFLECTION*1e6:.1f} µm threshold: {n_ok} / {len(rows)}")


def print_base_summary(rows, peak_force, rpm=None, n_flutes=None):
    """Print a concise narrative summary of the base-material comparison."""
    rpm = rpm if rpm is not None else DEFAULT_SPINDLE_RPM
    n_flutes = n_flutes if n_flutes is not None else DEFAULT_N_FLUTES
    f_tooth = tooth_passing_frequency(rpm, n_flutes)
    print("\n--- Base material narrative summary ---")
    print(
        f"At {rpm:.0f} rpm with {n_flutes} flutes, the cutter hammers the workpiece "
        f"{f_tooth:.0f} times per second ({f_tooth:.0f} Hz). The base's own natural "
        f"frequency is what it would 'ring' at if you tapped it."
    )
    print(
        "Because the base is much stiffer than the tool in this simplified model, "
        "the total deflection is mostly tool bending. The base still matters for "
        "damping and machine-level vibration."
    )
    for r in rows:
        print(
            f"  {r['base_material']:<13}: δ_static = {r['static_delta_um']:.3f} µm, "
            f"δ_dynamic = {r['dynamic_delta_um']:.3f} µm, f_n = {r['f_n_Hz']:.1f} Hz, "
            f"ζ = {r['damping_ratio']:.4f}, H = {r['amplification']:.2f}, "
            f"FoM = {r['figure_of_merit']:.2f}"
        )
    print("\n" + _interpret_base_material(rows))


def _format_analytical_results(tool_material, diameter, length, max_load, min_load,
                                material, feed_per_tooth, axial_depth):
    """Return human-readable strings for the static stiffness and cutting-force section.

    Plain-English idea:
        The tool acts like a diving board. A fatter or shorter diving board is
        stiffer. Stiffness is measured in N/m (Newtons of force per meter of
        bend). Deflection is how far the tip moves under a given push.
    """
    stiffness = tool_stiffness(tool_material, diameter, length)
    max_deflection = deflection(tool_material, diameter, length, max_load)
    min_deflection = deflection(tool_material, diameter, length, min_load)

    max_verdict = "exceeds" if max_deflection > MIN_DEFLECTION else "within"
    min_verdict = "exceeds" if min_deflection > MIN_DEFLECTION else "within"

    tool_stiffness_msg = (
        f"Tool stiffness ({tool_material}, D={diameter*1e3:.1f} mm, "
        f"L={length*1e3:.1f} mm): {stiffness:.2e} N/m"
    )
    max_deflection_msg = (
        f"Max deflection at {max_load} N: {max_deflection:.2e} m "
        f"({max_verdict} threshold)"
    )
    min_deflection_msg = (
        f"Min deflection at {min_load} N: {min_deflection:.2e} m "
        f"({min_verdict} threshold)"
    )

    phi = math.radians(90)
    F_t, F_r, F_a = tooth_forces(material, feed_per_tooth, axial_depth, phi)
    F_res = resultant_cutting_force(material, feed_per_tooth, axial_depth, phi)

    cutting_example_msg = (
        f"Cutting example: workpiece={material}, tool={tool_material}, "
        f"f_z={feed_per_tooth*1e6:.1f} µm/tooth, a_p={axial_depth*1e3:.2f} mm"
    )
    cutting_forces_msg = f"  F_t={F_t:.2f} N, F_r={F_r:.2f} N, F_a={F_a:.2f} N"
    resultant_cutting_msg = f"  Resultant force per tooth: {F_res:.2f} N"

    return {
        "tool_stiffness_msg": tool_stiffness_msg,
        "max_deflection_msg": max_deflection_msg,
        "min_deflection_msg": min_deflection_msg,
        "cutting_example_msg": cutting_example_msg,
        "cutting_forces_msg": cutting_forces_msg,
        "resultant_cutting_msg": resultant_cutting_msg,
    }


def _print_analytical_results(msgs):
    """Print the static analysis section with a plain-English intro."""
    print("\n--- Analytical results ---")
    print("How stiff is the tool, and how far does it bend under load?")
    for key in [
        "tool_stiffness_msg",
        "max_deflection_msg",
        "min_deflection_msg",
    ]:
        print(msgs[key])
    print("\nHow hard does the workpiece push back while cutting?")
    for key in [
        "cutting_example_msg",
        "cutting_forces_msg",
        "resultant_cutting_msg",
    ]:
        print(msgs[key])


def _run_reference_simulation(material, tool_material, diameter, length,
                              feed_per_tooth, axial_depth, rpm, stop_time, n_flutes,
                              run_simulink=False):
    """Run either the Simulink or Python-native reference simulation.

    Returns the time-domain results dict and prints a short summary.
    """
    if run_simulink:
        print("\n--- Simulink reference model ---")
        try:
            summary = simulink_summary(
                material=material,
                diameter=diameter,
                length=length,
                rpm=rpm,
                feed_per_tooth=feed_per_tooth,
                axial_depth=axial_depth,
                stop_time=stop_time,
                tool_material=tool_material,
            )
            print("Simulink final signal values:")
            for key, value in summary.items():
                print(f"  {key}: {value:.4e}")
            return None  # no Python time-domain results when Simulink succeeds
        except Exception as e:
            print(f"Simulink run failed: {e}")
            print("Falling back to Python-native simulation.")

    print(f"\n--- Python-native reference simulation (N_f={n_flutes}) ---")
    results = run_python_simulation(
        material=material,
        diameter=diameter,
        length=length,
        rpm=rpm,
        feed_per_tooth=feed_per_tooth,
        axial_depth=axial_depth,
        stop_time=stop_time,
        n_flutes=n_flutes,
        tool_material=tool_material,
    )
    print(f"  Peak resultant force: {results['F_res'].max():.2f} N")
    print(f"  Peak deflection: {results['deflection'].max()*1e6:.2f} µm")
    print(f"  Peak temperature rise: {results['temperature'].max():.2f} K")
    if not run_simulink:
        print("\n(Simulink demo disabled; set RUN_SIMULINK = True to run it.)")
    return results


def _summarize_with_ollama(results, base_rows, sweep_rows, msgs, vib_summary, chatter,
                           material, tool_material, diameter, length,
                           feed_per_tooth, axial_depth, n_flutes):
    """Ask a local LLM to write a plain-language summary of the results.

    This is purely for convenience. If Ollama is not running or the model is
    unavailable, the error is reported and the script continues.
    """
    try:
        import ollama
    except Exception as e:
        print(f"\n(Ollama summary skipped: {e})")
        return

    best_config = sweep_rows[0] if sweep_rows else None
    best_base_stiffness = min(base_rows, key=lambda r: r["total_delta_um"])
    best_base_fom = max(base_rows, key=lambda r: r["figure_of_merit"])

    payload_items = [
        f"Simulation goal: desktop CNC mill targeting sub-{MIN_DEFLECTION*1e6:.0f}-micron precision.",
        f"Default tool: {tool_material}, diameter={diameter*1e3:.1f} mm, length={length*1e3:.1f} mm, {n_flutes} flutes.",
        f"Default workpiece: {material}, feed_per_tooth={feed_per_tooth*1e6:.1f} µm/tooth, "
        f"axial_depth={axial_depth*1e3:.2f} mm, rpm={DEFAULT_SPINDLE_RPM}.",
        msgs["tool_stiffness_msg"],
        msgs["max_deflection_msg"],
        msgs["min_deflection_msg"],
        msgs["cutting_example_msg"],
        msgs["cutting_forces_msg"],
        msgs["resultant_cutting_msg"],
    ]

    if results is not None:
        payload_items += [
            f"Python-native peak resultant force: {results['F_res'].max():.2f} N",
            f"Python-native peak deflection: {results['deflection'].max()*1e6:.2f} µm",
            f"Python-native peak temperature rise: {results['temperature'].max():.2f} K",
        ]
    if vib_summary is not None:
        payload_items += [
            f"Tool natural frequency: {vib_summary['f_natural']:.1f} Hz; "
            f"exact cantilever: {vib_summary['f_natural_exact']:.1f} Hz.",
            f"Dominant vibration frequency: {vib_summary['deflection_peak_freq']:.1f} Hz.",
            f"Coarse chatter warning: {chatter}",
        ]

    payload_items += [
        "Base material comparison (best for stiffness / best balanced): "
        f"{best_base_stiffness['base_material']} (δ_total={best_base_stiffness['total_delta_um']:.2f} µm) / "
        f"{best_base_fom['base_material']} (FoM={best_base_fom['figure_of_merit']:.2f}).",
    ]

    if best_config is not None:
        payload_items += [
            "Top geometry/base configuration meeting target: "
            f"D={best_config['diameter_mm']:.2f} mm, L={best_config['length_mm']:.1f} mm, "
            f"base={best_config['base_material']}, "
            f"δ_total={best_config['total_delta_um']:.2f} µm, "
            f"status={('OK' if best_config['meets_threshold'] else 'EXCEEDS')}.",
            f"Out of {len(sweep_rows)} geometry/base combinations, "
            f"{sum(1 for r in sweep_rows if r['meets_threshold'])} meet the "
            f"{MIN_DEFLECTION*1e6:.0f} µm threshold.",
        ]

    payload_items += [
        "Simulink bridge status: matlabengine==24.2.2 matches MATLAB R2024b; engine starts; "
        "Simulink is licensed; but Simulink toolbox files are not on the MATLAB Engine path "
        "on this macOS+Conda setup, so the Python-Simulink bridge is not functional here. "
        "Recommended fallback: use the Python-native simulation or run the model directly "
        "inside the MATLAB desktop.",
    ]

    prompt = (
        "Provide a professional, concise, clear, and simple interpretation of these CNC mill "
        "simulation findings. Explain what matters most for a University of Pennsylvania "
        "Mechanical Engineering Senior Design project on a desktop sub-10-micron precision CNC mill. "
        "Include: (1) whether the default tool geometry meets the target, (2) what the dominant "
        "error source is, (3) which base material looks best and why, (4) what tool/base "
        "configuration meets the target, and (5) what the next engineering steps should be. "
        "Base all numerical claims strictly on the findings below; do not round or invent values."
        "\n\nFindings:\n" + "\n".join(f"- {item}" for item in payload_items)
    )

    client = ollama.Client()
    print("\n--- Ollama summary of findings ---")
    try:
        for part in client.chat("lfm2.5-thinking:latest", messages=[{"role": "user", "content": prompt}], stream=True):
            print(part.message.content, end="", flush=True)
        print("\n")
    except Exception as e:
        print(f"\n(Ollama summary failed during generation: {e})")


def run_all(run_simulink=False, summarize=True):
    """Run the full default analysis and print all comparison tables.

    Parameters
    ----------
    run_simulink : bool
        If True, attempt to use the Simulink bridge. Otherwise (default) use the
        Python-native simulation, which is reliable on this machine.
    summarize : bool
        If True, ask the local Ollama LLM to write a plain-language summary
        of the findings at the end. Set False to skip the LLM call and save
        time when you only need the raw tables.
    """
    plot_material_properties(disp=False, save=False)
    plot_material_costs(disp=False, save=False)

    # Default test case
    diameter = 0.001  # m
    length = 0.010    # m
    max_load = MAX_CUTTING_LOAD
    min_load = MIN_CUTTING_LOAD
    material = "aluminum"       # workpiece material
    tool_material = "wc_co"     # cutting tool material
    feed_per_tooth = 2e-6       # m/tooth (2 µm/tooth)
    axial_depth = 0.0001        # m (0.1 mm)
    n_flutes = DEFAULT_N_FLUTES

    # Static stiffness and instantaneous cutting forces
    msgs = _format_analytical_results(
        tool_material, diameter, length, max_load, min_load,
        material, feed_per_tooth, axial_depth
    )
    _print_analytical_results(msgs)

    # Time-domain reference simulation
    results = _run_reference_simulation(
        material, tool_material, diameter, length,
        feed_per_tooth, axial_depth, DEFAULT_SPINDLE_RPM, DEFAULT_STOP_TIME,
        n_flutes, run_simulink=run_simulink
    )

    # Vibration analysis for the primary material
    vib_summary = None
    chatter = None
    if results is not None:
        vib_summary = vibration_analysis(results, n_flutes=n_flutes)
        print_vibration_summary(vib_summary)
        chatter = chatter_indicator(results, n_flutes=n_flutes)
        print(f"  Coarse chatter warning (n_flutes={n_flutes}): {chatter}")
        plot_vibration_analysis(vib_summary, save=True)
        print("  Saved vibration spectrum plot to vibration_analysis.png")
    else:
        print("  (Vibration analysis skipped because no Python time-domain results are available.)")

    # Comparisons across workpiece materials
    vib_rows = run_vibration_comparison(
        diameter=diameter, length=length, n_flutes=n_flutes,
        tool_material=tool_material
    )
    print_vibration_comparison(vib_rows)

    rows = run_comparison(
        diameter=diameter, length=length, n_flutes=n_flutes,
        tool_material=tool_material
    )
    print_comparison(rows)

    # Base/frame material comparison
    peak_force = results["F_res"].max() if results is not None else MIN_CUTTING_LOAD
    base_rows = run_base_comparison(
        peak_force=peak_force,
        rpm=DEFAULT_SPINDLE_RPM,
        n_flutes=n_flutes,
        tool_material=tool_material,
        diameter=diameter,
        length=length,
    )
    print_base_comparison(base_rows, peak_force=peak_force,
                          rpm=DEFAULT_SPINDLE_RPM, n_flutes=n_flutes)
    print_base_summary(base_rows, peak_force=peak_force,
                       rpm=DEFAULT_SPINDLE_RPM, n_flutes=n_flutes)

    # Geometry/base sweep to find configurations that meet the target
    sweep_rows = run_geometry_base_sweep(
        peak_force=peak_force,
        tool_material=tool_material,
    )
    print_geometry_base_sweep(sweep_rows, peak_force=peak_force, top_n=20)
    print(
        "\nWhat this means: this sweep holds the cutting force fixed and tries many "
        "tool sizes and base materials. 'OK' means the static deflection budget is "
        "under the 10 µm target. Notice that short, thick tools pass; long, thin ones "
        "fail, no matter the base material."
    )

    print("\nThreshold: MIN_DEFLECTION =", MIN_DEFLECTION * 1e6, "µm")

    # Optional plain-language summary via local LLM.
    if summarize:
        _summarize_with_ollama(
            results, base_rows, sweep_rows, msgs, vib_summary, chatter,
            material, tool_material, diameter, length,
            feed_per_tooth, axial_depth, n_flutes
        )
    else:
        print("\n(LLM summary skipped; set summarize=True to enable.)")


def _interpret_deflection_status(deflection_m):
    """Return a plain-English sentence comparing deflection to the 10 µm target."""
    threshold_um = MIN_DEFLECTION * 1e6
    deflection_um = deflection_m * 1e6
    ratio = deflection_um / threshold_um
    if deflection_um <= threshold_um:
        return f"This is below the {threshold_um:.0f} µm target (ratio = {ratio:.2f})."
    else:
        return (
            f"This EXCEEDS the {threshold_um:.0f} µm target by a factor of "
            f"{ratio:.1f}×. Use a stiffer tool geometry (larger D, shorter L)."
        )


def _interpret_vibration(vib, n_flutes):
    """Return a plain-English summary of the vibration analysis."""
    f_tooth = vib["f_tooth"]
    f_n = vib["f_natural"]
    amp_tooth = vib["amp_tooth"]
    dominant_freq = vib["deflection_peak_freq"]

    if f_tooth / f_n < 0.5:
        proximity = "far below"
        risk = "low"
    elif f_tooth / f_n < 0.8:
        proximity = "below but approaching"
        risk = "moderate"
    else:
        proximity = "close to"
        risk = "potentially concerning"

    lines = [
        "What this means: the cutter hammers the tool at the tooth-passing "
        f"frequency ({f_tooth:.0f} Hz).",
        f"The tool's natural frequency is {f_n:.0f} Hz, so the hammering is "
        f"{proximity} resonance.",
        f"Dynamic amplification at the tooth-passing frequency is {amp_tooth:.3f} "
        "(1.0 = no extra vibration from resonance).",
        f"The strongest vibration in this simulation occurs at {dominant_freq:.0f} Hz.",
        f"Overall vibration risk: {risk}.",
    ]
    return "\n  ".join(lines)


def _interpret_base_material(rows):
    """Return a plain-English summary of the base-material comparison."""
    best_stiffness = min(rows, key=lambda r: r["total_delta_um"])
    best_damping = max(rows, key=lambda r: r["damping_ratio"])
    best_fom = max(rows, key=lambda r: r["figure_of_merit"])

    tool_dominates = all(r["total_delta_um"] - r["static_delta_um"] > 10 * r["static_delta_um"]
                         for r in rows)

    lines = [
        "What this means: the machine base acts like a second, much stiffer diving board.",
    ]
    if tool_dominates:
        lines.append(
            "For every base material, the total deflection is dominated by the tool, "
            "not the base. Base material mainly changes damping and natural frequency."
        )
    lines += [
        f"Stiffest choice: {best_stiffness['base_material']} "
        f"(total deflection {best_stiffness['total_delta_um']:.2f} µm).",
        f"Best damping: {best_damping['base_material']} "
        f"(damping ratio {best_damping['damping_ratio']:.4f}).",
        f"Best balanced score: {best_fom['base_material']} "
        f"(figure of merit {best_fom['figure_of_merit']:.2f}).",
    ]
    return "\n  ".join(lines)


def run_sub_ten_micron_demo():
    """Run a focused demo using the tool/base combo that meets the 10 µm budget.

    Think of this as: "show me the simplest design that actually passes."

    From the geometry/base sweep, the top configuration is roughly:
        - wc_co (tungsten carbide) tool, D = 3 mm, L = 5 mm
        - steel machine base
        - aluminum workpiece
    This gives a total static deflection well under 10 µm at the representative
    peak cutting force.
    """
    print("\n" + "=" * 70)
    print("Focused sub-10 µm precision demo")
    print("=" * 70)
    print(
        "Goal: find ONE tool/base combination that keeps total deflection under "
        f"{MIN_DEFLECTION*1e6:.0f} µm."
    )

    material = "aluminum"
    diameter = 0.003       # 3 mm
    length = 0.005         # 5 mm stick-out
    base_material = "steel"
    feed_per_tooth = 2e-6  # 2 µm/tooth
    axial_depth = 0.0001   # 0.1 mm
    rpm = DEFAULT_SPINDLE_RPM
    n_flutes = DEFAULT_N_FLUTES

    print(f"\nBase configuration: D={diameter*1e3:.1f} mm, L={length*1e3:.1f} mm, base={base_material}")
    print("Representative peak cutting force: 26.65 N (from the default 1 mm × 10 mm case)\n")

    print("| Tool material | E (GPa) | k_tool (N/m) | δ_tool (µm) | δ_total (µm) | f_n (Hz) |")
    print("|---------------|---------|--------------|-------------|--------------|----------|")

    for tool_material in tool_material_properties:
        k_tool = tool_stiffness(tool_material, diameter, length)
        delta_tool = deflection(tool_material, diameter, length, 26.65)
        delta_total = total_static_deflection(
            tool_material, base_material, 26.65, diameter, length
        )
        sdof = sdof_tool_parameters(tool_material, diameter, length)
        print(
            f"| {tool_material:<13} "
            f"| {tool_material_properties[tool_material]['E']/1e9:>7.0f} "
            f"| {k_tool:.2e} "
            f"| {delta_tool*1e6:>11.2f} "
            f"| {delta_total*1e6:>12.2f} "
            f"| {sdof['f_n']:>8.1f} |"
        )

    # Explain tool material trade-off in plain English.
    print(
        "\nInterpretation: for this short, thick geometry, BOTH tool materials pass. "
        "Carbide (wc_co) is stiffer, but even high-speed steel (hss) stays under 2 µm. "
        "Tool material matters, but geometry matters more."
    )

    # Generate and save Tool Diameter vs. Deflection plot
    plot_diameter_vs_deflection(force=26.65, tool_material="wc_co", disp=False, save=True)

    # Continue with the wc_co time-domain simulation
    tool_material = "wc_co"
    print(f"\n--- Running time-domain simulation for {tool_material} ---")
    results = run_python_simulation(
        material=material,
        tool_material=tool_material,
        diameter=diameter,
        length=length,
        rpm=rpm,
        feed_per_tooth=feed_per_tooth,
        axial_depth=axial_depth,
        stop_time=DEFAULT_STOP_TIME,
        n_flutes=n_flutes,
    )
    print(f"  Peak resultant force: {results['F_res'].max():.2f} N")
    print(f"  Peak deflection:      {results['deflection'].max()*1e6:.2f} µm")
    print(f"  Peak temperature rise: {results['temperature'].max():.2f} K")
    print("  " + _interpret_deflection_status(results["deflection"].max()))

    vib = vibration_analysis(results, n_flutes=n_flutes)
    print_vibration_summary(vib)
    print("  " + _interpret_vibration(vib, n_flutes))
    chatter = chatter_indicator(results, n_flutes=n_flutes)
    print(f"  Coarse chatter warning: {chatter}")
    print("  (A 'False' warning means the simulation did not see resonance-driven vibration.)")

    return results


if __name__ == "__main__":
    # Set True to attempt the Simulink bridge (requires a working MATLAB/Simulink setup).
    # On this macOS + Conda Python machine the Python-native simulation is the reliable path.
    RUN_SIMULINK = False
    # Set True to ask the local Ollama LLM for a plain-language summary at the end.
    SUMMARIZE = False
    # Set True to run the focused sub-10 µm demo instead of the full default report.
    SUB_TEN_MICRON_DEMO = True

    if SUB_TEN_MICRON_DEMO:
        run_sub_ten_micron_demo()
    else:
        run_all(run_simulink=RUN_SIMULINK, summarize=SUMMARIZE)
