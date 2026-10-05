"""
File-based Python -> MATLAB bridge for the CNC Simulink reference model.

This avoids keeping the MATLAB Engine process alive. Python writes a .mat file
with all simulation parameters, then invokes MATLAB in batch mode to build the
model (if needed), run it, and save the results to another .mat file. Python then
reads the results back.

Usage:
    from run_simulink_batch import simulink_batch_summary

    summary = simulink_batch_summary(
        material="aluminum",
        diameter=0.001,
        length=0.010,
        rpm=60000,
        feed_per_tooth=2e-6,
        axial_depth=0.0001,
        stop_time=0.1,
    )
    print(summary)
"""

import os
import subprocess
import tempfile

import scipy.io

from sims import (
    DEFAULT_TOOL_MATERIAL,
    SIMULINK_MODEL_NAME,
    _simulink_workspace_params,
    run_python_simulation,
)


def _matlab_executable():
    """Return a likely MATLAB executable path on macOS, Windows, or Linux."""
    candidates = [
        "/Applications/MATLAB_R2024b.app/bin/matlab",  # installed version
        "/Applications/MATLAB_R2026a.app/bin/matlab",
        "/Applications/MATLAB_R2025b.app/bin/matlab",
        "/Applications/MATLAB_R2025a.app/bin/matlab",
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate

    # Fall back to PATH lookup
    for executable in ["matlab", "matlab.exe"]:
        found = subprocess.run(["which", executable], capture_output=True, text=True)
        if found.returncode == 0 and found.stdout.strip():
            return found.stdout.strip()

    raise FileNotFoundError("MATLAB executable not found. Please add it to PATH.")


def simulink_batch_summary(material, diameter, length, rpm=None, feed_per_tooth=None,
                           axial_depth=None, stop_time=None, cleanup=True,
                           tool_material=DEFAULT_TOOL_MATERIAL):
    """Run the Simulink reference model via a file-based MATLAB batch call.

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
    rpm, feed_per_tooth, axial_depth, stop_time : float, optional
        Process parameters.
    cleanup : bool, optional
        If True, delete temporary .mat files after reading results.

    Returns
    -------
    dict
        Final values of each logged Simulink signal.
    """
    project_dir = os.path.dirname(os.path.abspath(__file__))
    simulink_dir = os.path.join(project_dir, "simulink")
    matlab = _matlab_executable()

    params = _simulink_workspace_params(
        material, diameter, length, rpm, feed_per_tooth, axial_depth, stop_time,
        tool_material
    )

    # Use /tmp for temp files because MATLAB's IPC socket can exceed the
    # max path length (104 chars) under macOS temp directories.
    tmpdir = tempfile.mkdtemp(prefix="cncsim_", dir="/tmp")
    try:
        inputs_path = os.path.join(tmpdir, "in.mat")
        outputs_path = os.path.join(tmpdir, "out.mat")
        matlab_log = os.path.join(tmpdir, "log.txt")
        runner_m = os.path.join(tmpdir, "run_cnc_batch.m")

        scipy.io.savemat(inputs_path, params)

        # Write a temporary MATLAB runner script instead of a long -batch
        # command; this avoids command-line length limits and quoting issues.
        runner_script = (
            f"try;\n"
            f"if ~license('test', 'Simulink');\n"
            f"    error('Simulink license is not available in this MATLAB session.');\n"
            f"end;\n"
            # explicitly load Simulink; harmless if already loaded
            f"simulink;\n"
            f"addpath('{project_dir}');\n"
            f"addpath('{simulink_dir}');\n"
            f"load('{inputs_path}');\n"
            f"if ~exist('{SIMULINK_MODEL_NAME}.slx', 'file');\n"
            f"    build_cnc_simulink_model('{SIMULINK_MODEL_NAME}');\n"
            f"end;\n"
            f"results = run_cnc_simulink_sim('{SIMULINK_MODEL_NAME}');\n"
            f"save('{outputs_path}', 'results');\n"
            f"catch ME;\n"
            f"disp(getReport(ME));\n"
            f"exit(1);\n"
            f"end;\n"
            f"exit(0);\n"
        )
        with open(runner_m, "w") as f:
            f.write(runner_script)

        cmd = [matlab, "-nodisplay", "-nodesktop", "-r", f"run('{runner_m}')"]

        print(f"Invoking MATLAB batch job for {material}...")

        # Short TMPDIR is critical on macOS to avoid MathWorks IPC socket
        # path-length errors (error 5201).
        env = os.environ.copy()
        env["TMPDIR"] = "/tmp"
        env["TEMP"] = "/tmp"
        env["TMP"] = "/tmp"

        with open(matlab_log, "w") as logfile:
            result = subprocess.run(
                cmd, stdout=logfile, stderr=subprocess.STDOUT, text=True, env=env
            )

        if result.returncode != 0:
            with open(matlab_log, "r") as logfile:
                log = logfile.read()
            extra = ""
            if "Unrecognized function or variable 'simulink'" in log:
                extra = (
                    "\n\nNote: Simulink could not be loaded in the MATLAB batch session. "
                    "Some MATLAB licenses do not allow Simulink in headless/batch mode. "
                    "Try running the model directly inside the MATLAB desktop, or use the "
                    "Python-native simulation in sims.py."
                )
            raise RuntimeError(
                f"MATLAB batch job failed (exit {result.returncode}). Log:\n{log}{extra}"
            )

        outputs = scipy.io.loadmat(outputs_path, simplify_cells=True)
        results = outputs["results"]

        summary = {}
        for name, signal in results.items():
            if name.startswith("__"):
                continue
            if isinstance(signal, dict) and "data" in signal:
                data = signal["data"]
                if hasattr(data, "size") and data.size:
                    summary[name] = float(data[-1])
                else:
                    summary[name] = None
            else:
                summary[name] = signal

        return summary
    finally:
        if cleanup:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    material = "aluminum"
    tool_material = "wc_co"
    diameter = 0.001
    length = 0.010
    rpm = 60000
    feed_per_tooth = 2e-6
    axial_depth = 0.0001
    stop_time = 0.1
    n_flutes = 2

    try:
        summary = simulink_batch_summary(
            material=material,
            tool_material=tool_material,
            diameter=diameter,
            length=length,
            rpm=rpm,
            feed_per_tooth=feed_per_tooth,
            axial_depth=axial_depth,
            stop_time=stop_time,
        )
        print("Simulink final signal values:")
        for key, value in summary.items():
            print(f"  {key}: {value}")
    except Exception as e:
        print(f"Simulink run failed: {e}")
        print("Falling back to Python-native reference simulation...")
        results = run_python_simulation(
            material=material,
            tool_material=tool_material,
            diameter=diameter,
            length=length,
            rpm=rpm,
            feed_per_tooth=feed_per_tooth,
            axial_depth=axial_depth,
            stop_time=stop_time,
            n_flutes=n_flutes,
        )
        print(f"  Peak resultant force: {results['F_res'].max():.2f} N")
        print(f"  Peak deflection: {results['deflection'].max()*1e6:.2f} µm")
        print(f"  Peak temperature rise: {results['temperature'].max():.2f} K")
