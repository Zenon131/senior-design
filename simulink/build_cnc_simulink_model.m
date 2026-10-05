function build_cnc_simulink_model(modelName)
%BUILD_CNC_SIMULINK_MODEL Create a reference CNC mill Simulink model.
%
%   build_cnc_simulink_model('cnc_mill_sim') creates a Simulink model with:
%       - spindle RPM input
%       - mechanistic cutting-force generation (tangential, radial, axial)
%       - tool deflection dynamics (mass-spring-damper)
%       - simplified thermal dynamics (first-order RC)
%       - scopes for force, deflection, and temperature
%
%   The model is intentionally simple: it simulates one cutting tooth
%   passing through the workpiece once per revolution. Parameters are set
%   from the MATLAB base workspace by the Python caller before simulation.

    if nargin < 1 || isempty(modelName)
        modelName = 'cnc_mill_sim';
    end

    % Explicitly load Simulink so the model-creation API is available in
    % batch / nodisplay MATLAB sessions.
    simulink;

    % Verify Simulink is licensed.
    if ~license('test', 'Simulink')
        error('Simulink license is not available on this MATLAB session.');
    end

    % Delete any existing model file on disk.
    if exist([modelName '.slx'], 'file')
        delete([modelName '.slx']);
    end

    % If a model with this name is already in memory, close it first.
    % Use try/catch because close_system needs Simulink loaded.
    try
        close_system(modelName, 0);
    catch ME %#ok<NASGU>
        % Model was not open; continue.
    end

    % Create a new blank model. This call loads Simulink automatically.
    try
        new_system(modelName);
    catch ME
        error('Failed to create Simulink model %s: %s', modelName, ME.message);
    end
    open_system(modelName);
    set_param(modelName, 'StopTime', '0.1');
    set_param(modelName, 'Solver', 'ode45');
    set_param(modelName, 'MaxStep', '1e-5');

    % Add blocks
    % Clock / time source
    add_block('simulink/Sources/Clock', [modelName '/Clock'], ...
        'Position', [30 50 60 80]);

    % RPM constant
    add_block('simulink/Sources/Constant', [modelName '/RPM'], ...
        'Position', [30 120 60 150], ...
        'Value', 'rpm', ...
        'Commented', 'off');

    % Feed per tooth constant
    add_block('simulink/Sources/Constant', [modelName '/feed_per_tooth'], ...
        'Position', [30 190 60 220], ...
        'Value', 'feed_per_tooth');

    % Axial depth constant
    add_block('simulink/Sources/Constant', [modelName '/axial_depth'], ...
        'Position', [30 260 60 290], ...
        'Value', 'axial_depth');

    % Rotation angle = 2*pi*(rpm/60)*t
    add_block('simulink/Math Operations/Gain', [modelName '/rad_per_sec'], ...
        'Position', [120 120 160 150], ...
        'Gain', '2*pi/60');
    add_block('simulink/Math Operations/Product', [modelName '/phi'], ...
        'Position', [220 70 250 100], ...
        'Inputs', '**');

    % MATLAB Function block: cutting forces
    fcnLines = {
        'function [F_t, F_r, F_a] = fcn(phi, f_z, a_p, K_tc, K_rc, K_ac, K_te, K_re, K_ae)'
        'h = max(0.0, f_z * sin(phi));'
        'F_t = K_tc * h * a_p + K_te * a_p;'
        'F_r = K_rc * h * a_p + K_re * a_p;'
        'F_a = K_ac * h * a_p + K_ae * a_p;'
        'end'
    };
    fcnCode = strjoin(fcnLines, sprintf('\n'));

    add_block('simulink/User-Defined Functions/MATLAB Function', ...
        [modelName '/CuttingForces'], ...
        'Position', [320 70 470 200], ...
        'MATLABFunction', fcnCode);

    % Cutting coefficient constants (wired from workspace)
    add_block('simulink/Sources/Constant', [modelName '/K_tc'], ...
        'Position', [30 330 60 360], 'Value', 'K_tc');
    add_block('simulink/Sources/Constant', [modelName '/K_rc'], ...
        'Position', [30 370 60 400], 'Value', 'K_rc');
    add_block('simulink/Sources/Constant', [modelName '/K_ac'], ...
        'Position', [30 410 60 440], 'Value', 'K_ac');
    add_block('simulink/Sources/Constant', [modelName '/K_te'], ...
        'Position', [30 450 60 480], 'Value', 'K_te');
    add_block('simulink/Sources/Constant', [modelName '/K_re'], ...
        'Position', [30 490 60 520], 'Value', 'K_re');
    add_block('simulink/Sources/Constant', [modelName '/K_ae'], ...
        'Position', [30 530 60 560], 'Value', 'K_ae');

    % Tool deflection dynamics: 1/k mass-spring-damper
    % x'' = (F - c*x' - k*x)/m  =>  transfer function  1/(m*s^2 + c*s + k)
    add_block('simulink/Continuous/Transfer Fcn', [modelName '/ToolDynamics'], ...
        'Position', [550 90 650 140], ...
        'Numerator', '1', ...
        'Denominator', '[tool_mass tool_damping tool_stiffness]', ...
        'Commented', 'off');

    % Resultant force
    add_block('simulink/Math Operations/Vector Concatenate', ...
        [modelName '/F_concat'], ...
        'Position', [520 220 550 250], ...
        'NumInputs', '3');
    add_block('simulink/Math Operations/Fcn', [modelName '/F_res'], ...
        'Position', [620 220 660 250], ...
        'Expr', 'sqrt(u(1)^2 + u(2)^2 + u(3)^2)');

    % Thermal dynamics: C_th * dT/dt = Q - T/R_th
    % Use Q = F_t * V_c as heat input, first-order transfer function
    % T(s)/Q(s) = R_th / (R_th*C_th*s + 1)
    add_block('simulink/Math Operations/Product', [modelName '/HeatPower'], ...
        'Position', [750 90 780 120], ...
        'Inputs', '**');
    add_block('simulink/Continuous/Transfer Fcn', [modelName '/ThermalDynamics'], ...
        'Position', [820 90 920 140], ...
        'Numerator', '[thermal_resistance]', ...
        'Denominator', '[thermal_resistance*thermal_capacitance 1]');

    % Scopes
    add_block('simulink/Sinks/Scope', [modelName '/ForceScope'], ...
        'Position', [720 220 760 250], ...
        'NumInputPorts', '1');
    add_block('simulink/Sinks/Scope', [modelName '/DeflectionScope'], ...
        'Position', [720 90 760 140], ...
        'NumInputPorts', '1');
    add_block('simulink/Sinks/Scope', [modelName '/TempScope'], ...
        'Position', [960 90 1000 140], ...
        'NumInputPorts', '1');

    % To Workspace blocks for post-processing (FFT, stability analysis)
    add_block('simulink/Sinks/To Workspace', [modelName '/ForceLog'], ...
        'Position', [820 220 860 250], ...
        'VariableName', 'F_res_log', ...
        'SaveFormat', 'Array');
    add_block('simulink/Sinks/To Workspace', [modelName '/DeflectionLog'], ...
        'Position', [820 90 860 140], ...
        'VariableName', 'deflection_log', ...
        'SaveFormat', 'Array');

    % Connections
    add_line(modelName, 'Clock/1', 'phi/1');
    add_line(modelName, 'RPM/1', 'rad_per_sec/1');
    add_line(modelName, 'rad_per_sec/1', 'phi/2');
    add_line(modelName, 'phi/1', 'CuttingForces/1');
    add_line(modelName, 'feed_per_tooth/1', 'CuttingForces/2');
    add_line(modelName, 'axial_depth/1', 'CuttingForces/3');
    add_line(modelName, 'K_tc/1', 'CuttingForces/4');
    add_line(modelName, 'K_rc/1', 'CuttingForces/5');
    add_line(modelName, 'K_ac/1', 'CuttingForces/6');
    add_line(modelName, 'K_te/1', 'CuttingForces/7');
    add_line(modelName, 'K_re/1', 'CuttingForces/8');
    add_line(modelName, 'K_ae/1', 'CuttingForces/9');

    % Cutting force components to dynamics
    add_line(modelName, 'CuttingForces/1', 'ToolDynamics/1');   % F_t drives deflection
    add_line(modelName, 'CuttingForces/1', 'F_concat/1');
    add_line(modelName, 'CuttingForces/2', 'F_concat/2');
    add_line(modelName, 'CuttingForces/3', 'F_concat/3');
    add_line(modelName, 'F_concat/1', 'F_res/1');
    add_line(modelName, 'F_res/1', 'ForceScope/1');
    add_line(modelName, 'F_res/1', 'ForceLog/1');
    add_line(modelName, 'ToolDynamics/1', 'DeflectionScope/1');
    add_line(modelName, 'ToolDynamics/1', 'DeflectionLog/1');

    % Add a heat-partition gain so the thermal model only sees the fraction of
    % cutting energy that enters the tool.
    add_block('simulink/Math Operations/Gain', [modelName '/HeatPartition'], ...
        'Position', [700 90 730 120], ...
        'Gain', 'heat_partition');

    % Heat = F_t * V_c, where V_c = pi*D*rpm/60
    add_block('simulink/Sources/Constant', [modelName '/cutting_speed'], ...
        'Position', [650 170 680 200], ...
        'Value', 'pi * tool_diameter * rpm / 60');
    add_line(modelName, 'CuttingForces/1', 'HeatPower/1');
    add_line(modelName, 'cutting_speed/1', 'HeatPower/2');
    add_line(modelName, 'HeatPower/1', 'HeatPartition/1');
    add_line(modelName, 'HeatPartition/1', 'ThermalDynamics/1');
    add_line(modelName, 'ThermalDynamics/1', 'TempScope/1');

    % Save the model
    save_system(modelName);
    fprintf('Simulink model ''%s'' created and saved.\n', modelName);
end
