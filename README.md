# MEAM Senior Design

## Primary Idea: Micron-Scale Mini CNC Mill

## Simulations

The project includes a Python simulation module (`sims.py`) for quick analytical estimates and a reference Simulink model for time-domain dynamics.

### Python-native reference simulation

Run directly from the project directory:

```bash
python sims.py
```

By default this runs the **sub-10 µm precision demo**, which shows the simplest
tool/base configuration that actually meets the 10 µm deflection target. To run
the full comparison report instead, set `SUB_TEN_MICRON_DEMO = False` near the
bottom of `sims.py`.

The output now includes plain-English "What this means" lines after each major
section, so the tables are easier to interpret without a mechanical-engineering
background.

Key tunable parameters in `sims.py`:

| Parameter | Default | Description |
|---|---|---|
| `material` | `"aluminum"` | Workpiece material being cut (key into `material_properties` and `cutting_coefficients`) |
| `tool_material` | `"wc_co"` | Cutting tool material (key into `tool_material_properties`) |
| `diameter` | `0.001` m | Tool diameter |
| `length` | `0.010` m | Tool stick-out / overhang |
| `feed_per_tooth` | `2e-6` m | Feed per tooth |
| `axial_depth` | `0.0001` m | Axial depth of cut |
| `DEFAULT_SPINDLE_RPM` | `60000` rpm | Spindle speed |
| `DEFAULT_N_FLUTES` | `2` | Number of flutes used in time-domain simulation |
| `DEFAULT_DAMPING_RATIO` | `0.05` | SDOF tool damping ratio |
| `HEAT_PARTITION` | `0.1` | Fraction of cutting heat entering the tool |
| `DEFAULT_BASE_WIDTH/HEIGHT/LENGTH` | `0.10 / 0.10 / 0.30` m | Default base cantilever dimensions |
| `DEFAULT_BASE_DAMPING_RATIO` | `0.05` | Fallback base damping ratio |
| `MIN_DEFLECTION` | `10e-6` m | Target deflection budget |

Material databases:

- `material_properties` — workpiece materials (`aluminum`, `brass`, `copper`, `gold`, `silver`, `steel`, `titanium`, `tungsten`).
- `tool_material_properties` — tool materials (`hss`, `wc_co`).
- `cutting_coefficients` — mechanistic cutting coefficients indexed by workpiece material.

### Python-native vibration analysis

`sims.py` includes a refined frequency-domain analysis of cutting force and tool deflection. The tool is treated as a single-degree-of-freedom (SDOF) oscillator:

```
m_eq * x'' + c * x' + k * x = F_t(t)
```

with parameters:

| Symbol | Formula | Meaning |
|---|---|---|
| `k` | `3 E I / L³` | Cantilever bending stiffness |
| `I` | `π D⁴ / 64` | Second moment of area |
| `m_total` | `ρ (π D² / 4) L` | Full cylinder mass |
| `m_eq` | `0.23 m_total` | Effective tip mass for first bending mode |
| `ω_n` | `sqrt(k / m_eq)` | Natural circular frequency |
| `f_n` | `ω_n / (2π)` | Natural frequency in Hz |
| `c` | `2 ζ sqrt(k m_eq)` | Viscous damping coefficient (ζ = 0.05 default) |

For comparison, the exact Euler-Bernoulli cantilever frequency is also computed:

```
f_1 = (β_1² / 2π) sqrt(E I / (ρ A L⁴)),   β_1 = 1.875
```

The FFT analysis reports:

- **Spindle frequency**: `f_spindle = rpm / 60`
- **Tooth-passing frequency**: `f_tooth = n_flutes * rpm / 60`
- **SDOF natural frequency** and **exact cantilever frequency**
- **Frequency ratios**: `r = f_excitation / f_n`
- **Dynamic amplification factor**: `H(r) = 1 / sqrt((1 - r²)² + (2 ζ r)²)`
- **Dominant frequencies** in force and deflection spectra
- **Coarse chatter warning**: flags when the dominant deflection frequency is near the natural frequency but not near any spindle or tooth-passing harmonic

Multi-flute forcing is included: each tooth j contributes with phase `φ_j(t) = ω t + j·2π/N_f`.

The script saves a spectrum plot to `vibration_analysis.png`.

You can run the analysis independently:

```python
from sims import run_python_simulation, vibration_analysis, print_vibration_summary

results = run_python_simulation(
    "aluminum",          # workpiece material
    diameter=0.001,
    length=0.010,
    rpm=60000,
    n_flutes=2,
    tool_material="wc_co",
)
vib = vibration_analysis(results, n_flutes=2)
print_vibration_summary(vib)
```

### Base / frame material selection

`sims.py` also compares candidate materials for the **machine base or frame** (`base_material_properties`):

| Material | E (GPa) | ζ (typical) | cost ($/kg) | Notes |
|---|---|---|---|---|
| aluminum | 70 | 0.0010 | 3.26 | light, stiff, low damping |
| steel | 210 | 0.0015 | 2.50 | highest stiffness, low damping |
| cast_iron | 110 | 0.0050 | 1.80 | good stiffness + damping balance |
| granite | 50 | 0.0100 | 0.50 | cheaper, moderate damping |
| concrete | 30 | 0.0150 | 0.10 | very cheap, good damping, low stiffness |
| epoxy_granite | 40 | 0.0202 | 2.00 | polymer-concrete composite; E and ζ from Abdellah et al. (2026, MDPI Polymers) |

The base is modeled as a rectangular cantilever beam with default geometry 100 mm × 100 mm × 300 mm. Key formulas:

```
I_base = w * h^3 / 12
k_base = 3 * E * I_base / L^3
m_base = ρ * w * h * L
f_n_base = (1 / 2π) * sqrt(k_base / m_base)
H(r) = 1 / sqrt((1 - r^2)^2 + (2 ζ r)^2)
δ_dynamic = H(r) * δ_static
δ_total = F/k_tool + F/k_base
```

`run_base_comparison()` returns a table of stiffness, static/dynamic deflection, natural frequency, damping ratio, dynamic amplification, cost, and a normalized **figure of merit (FoM)** that combines stiffness (50%), damping (30%), and cost (20%).

`run_geometry_base_sweep()` varies tool diameter, tool stick-out, and base material to find combinations that meet the `MIN_DEFLECTION` budget (default 10 µm). The sweep treats the tool and base as two springs in series and reports the total static deflection.

```python
from sims import run_base_comparison, print_base_comparison
rows = run_base_comparison(peak_force=26.65)
print_base_comparison(rows, peak_force=26.65)
```

### Simulink vibration analysis

When running Simulink, the model includes `To Workspace` blocks that log `F_res_log` and `deflection_log`. After a simulation, call:

```matlab
results = run_cnc_simulink_sim('cnc_mill_sim');
analyze_simulink_vibration(results, 2);
```

This produces the same kind of spectrum plot and summary from the Simulink data.

### Simulink reference model

MATLAB scripts live in `simulink/`:

- `simulink/build_cnc_simulink_model.m` — creates `cnc_mill_sim.slx` programmatically.
- `simulink/run_cnc_simulink_sim.m` — runs the model and returns logged signals.
- `simulink/analyze_simulink_vibration.m` — FFT and summary from Simulink logs.

You can run the model directly in MATLAB:

```matlab
addpath('simulink');
build_cnc_simulink_model('cnc_mill_sim');
results = run_cnc_simulink_sim('cnc_mill_sim');
analyze_simulink_vibration(results, 2);
```

#### MATLAB Engine setup (macOS)

1. Detect your MATLAB version and matching engine:

```bash
python diagnostics/setup_matlab_bridge.py
```

2. Install the matching `matlabengine` package. For MATLAB **R2024b**:

```bash
python -m pip install matlabengine==24.2.2
```

3. Test the bridge:

```bash
python diagnostics/test_matlab_bridge.py
```

This checks: MATLAB import, MATLAB startup, Simulink license, model build/load, and a short simulation.

⚠️ **macOS + Conda Python limitation:** The MATLAB Engine API can hit one of two problems:

1. A `libexpat` dylib conflict (MATLAB's `libexpat` shadows Conda Python's copy), causing `_XML_SetAllocTrackerActivationThreshold` errors.
2. Simulink is licensed, but its toolbox files are not on the MATLAB Engine path, causing `new_system requires Simulink` or `Undefined function 'new_system'`.

If `diagnostics/test_matlab_bridge.py` passes steps 1–3 but fails at step 4, try these workarounds in order:

- **Connect to an already-running MATLAB desktop** (often the only mode where Simulink is fully initialized on macOS):
  Open MATLAB, run `simulink` in its command window, then write a small Python script that calls `matlab.engine.connect_matlab()` instead of `start_matlab()`.
- **Run the Simulink model directly inside MATLAB** (not from Python):
  ```matlab
  addpath('simulink');
  build_cnc_simulink_model('cnc_mill_sim');
  results = run_cnc_simulink_sim('cnc_mill_sim');
  ```
- **Use Python from python.org or Homebrew** in a fresh venv instead of Conda.
- **Use the Python-native simulation in `sims.py` instead**; it implements the same physics.

You can call Simulink from Python in two ways:

#### 1. MATLAB Engine API (`sims.py`)

Set `RUN_SIMULINK = True` in the `__main__` block of `sims.py`. Requires `matlabengine`:

```bash
python -m pip install matlabengine==24.2.2
```

#### 2. File-based batch bridge (`run_simulink_batch.py`)

Does not keep a MATLAB engine process alive. Writes parameters to a `.mat` file, invokes MATLAB in batch mode, and reads results back:

```bash
python run_simulink_batch.py
```

If Simulink is unavailable, this script automatically falls back to the Python-native simulation.

### Known environment issues

- MATLAB/Simulink must be installed, licensed, and working for the Simulink path to run.
- The first time you run `python sims.py` with `RUN_SIMULINK = True`, MATLAB takes ~30–60 s to start.
- On macOS, long temporary-directory paths can cause MathWorks IPC socket errors. The batch bridge sets `TMPDIR=/tmp` to avoid this.
- If you see `Undefined function 'new_system'` or `License Error: Licensing shutdown`, Simulink is not loadable in your MATLAB session. Use the Python-native simulation instead until the license/environment is fixed.
- Some MATLAB licenses (especially student/trial or network licenses) do not allow Simulink in batch/headless mode, even though it works in interactive/engine mode. If `test_matlab_bridge.py` passes steps 1–3 but fails at step 4, this is likely the issue.

