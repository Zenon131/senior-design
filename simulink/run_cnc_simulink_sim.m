function results = run_cnc_simulink_sim(modelName)
%RUN_CNC_SIMULINK_SIM Run the reference CNC mill Simulink model.
%
%   results = run_cnc_simulink_sim('cnc_mill_sim') runs the model, extracts
%   the last values of the logged signals, and returns them in a struct.
%   The function assumes the workspace already contains all parameters
%   (rpm, feed_per_tooth, axial_depth, tool_diameter, tool_length, material,
%   K coefficients, tool_stiffness, tool_mass, tool_damping,
%   thermal_resistance, thermal_capacitance).

    if nargin < 1 || isempty(modelName)
        modelName = 'cnc_mill_sim';
    end

    % Verify Simulink is licensed.
    if ~license('test', 'Simulink')
        error('Simulink license is not available on this MATLAB session.');
    end

    % Load the model. This call loads Simulink automatically.
    try
        load_system(modelName);
    catch ME
        error('Failed to load Simulink model %s: %s', modelName, ME.message);
    end

    % Set simulation stop time from workspace if available, otherwise 0.1 s
    if exist('sim_stop_time', 'var')
        set_param(modelName, 'StopTime', num2str(sim_stop_time));
    else
        set_param(modelName, 'StopTime', '0.1');
    end

    % Run the simulation
    simOut = sim(modelName, 'SaveOutput', 'on', 'OutputSaveName', 'yout');

    % Extract final values from logged output (struct format)
    results = struct();
    if isstruct(simOut.yout)
        names = fieldnames(simOut.yout);
        for i = 1:numel(names)
            signal = simOut.yout.(names{i});
            if isstruct(signal) && isfield(signal, 'Data') && isfield(signal, 'Time')
                results.(names{i}) = struct('time', signal.Time, 'data', signal.Data);
            end
        end
    end

    % Also pull To Workspace logs if they exist
    if exist('F_res_log', 'var')
        results.F_res_log = struct('time', F_res_log.Time, 'data', F_res_log.Data);
    end
    if exist('deflection_log', 'var')
        results.deflection_log = struct('time', deflection_log.Time, 'data', deflection_log.Data);
    end

    % Compute a few vibration markers. tool_mass is already the effective
    % lumped mass m_eq used in the SDOF model.
    if exist('tool_stiffness', 'var') && exist('tool_mass', 'var') && exist('rpm', 'var')
        f_natural = (1 / (2 * pi)) * sqrt(tool_stiffness / tool_mass);
        f_spindle = rpm / 60;
        results.vibration_markers = struct('f_natural_hz', f_natural, 'f_spindle_hz', f_spindle);
    end

    fprintf('Simulink simulation completed for model ''%s''.\n', modelName);
end
