%ANALYZE_SIMULINK_VIBRATION Frequency-domain analysis of logged Simulink signals.
%
%   analyze_simulink_vibration(results, n_flutes) plots the FFT of the logged
%   cutting force and tool deflection signals, marks the spindle frequency,
%   tooth-passing frequency, and estimated tool natural frequency, and prints
%   a short summary.
%
%   Example:
%       results = run_cnc_simulink_sim('cnc_mill_sim');
%       analyze_simulink_vibration(results, 2);

function analyze_simulink_vibration(results, n_flutes)
    if nargin < 2
        n_flutes = 2;
    end

    % Pull logged signals
    if isfield(results, 'F_res_log')
        force = results.F_res_log;
    else
        error('No F_res_log found in results. Make sure the model has a To Workspace block.');
    end

    if isfield(results, 'deflection_log')
        deflection = results.deflection_log;
    else
        error('No deflection_log found in results. Make sure the model has a To Workspace block.');
    end

    t_force = force.time;
    y_force = force.data;
    t_defl = deflection.time;
    y_defl = deflection.data;

    % Sampling rates
    dt_force = mean(diff(t_force));
    dt_defl = mean(diff(t_defl));
    fs_force = 1 / dt_force;
    fs_defl = 1 / dt_defl;

    % Detrend and window
    y_force = y_force - mean(y_force);
    y_defl = y_defl - mean(y_defl);
    w_force = hanning(numel(y_force));
    w_defl = hanning(numel(y_defl));

    N_force = numel(y_force);
    N_defl = numel(y_defl);
    freqs_force = (0:N_force/2) * fs_force / N_force;
    freqs_defl = (0:N_defl/2) * fs_defl / N_defl;

    amp_force = abs(fft(y_force .* w_force));
    amp_force = amp_force(1:N_force/2 + 1) * 2 / sum(w_force);

    amp_defl = abs(fft(y_defl .* w_defl));
    amp_defl = amp_defl(1:N_defl/2 + 1) * 2 / sum(w_defl);

    % Frequency markers
    rpm = evalin('base', 'rpm');
    f_spindle = rpm / 60;
    f_tooth = n_flutes * f_spindle;

    % Natural frequency from the workspace. tool_mass is already the effective
    % lumped mass m_eq used in the SDOF model.
    tool_stiffness = evalin('base', 'tool_stiffness');
    tool_mass = evalin('base', 'tool_mass');
    f_natural = (1 / (2 * pi)) * sqrt(tool_stiffness / tool_mass);

    % Dominant peaks (skip DC)
    [~, idx_force] = max(amp_force(2:end));
    [~, idx_defl] = max(amp_defl(2:end));
    f_dom_force = freqs_force(idx_force + 1);
    f_dom_defl = freqs_defl(idx_defl + 1);

    % Print summary
    fprintf('\n--- Vibration analysis from Simulink ---\n');
    fprintf('  Sampling rate (force): %.1f Hz\n', fs_force);
    fprintf('  Spindle frequency: %.1f Hz\n', f_spindle);
    fprintf('  Tooth-passing frequency (%d flutes): %.1f Hz\n', n_flutes, f_tooth);
    fprintf('  Estimated tool natural frequency: %.1f Hz\n', f_natural);
    fprintf('  Dominant force frequency: %.1f Hz\n', f_dom_force);
    fprintf('  Dominant deflection frequency: %.1f Hz\n', f_dom_defl);

    % Coarse chatter warning
    tol = 0.05;
    near_natural = abs(f_dom_defl - f_natural) / f_natural < tol;
    near_spindle = abs(f_dom_defl - f_spindle) / f_spindle < tol;
    near_tooth = abs(f_dom_defl - f_tooth) / f_tooth < tol;
    chatter_warning = near_natural && ~(near_spindle || near_tooth);
    fprintf('  Coarse chatter warning: %s\n', string(chatter_warning));

    % Plot
    figure('Name', 'Simulink vibration analysis');

    subplot(2, 1, 1);
    plot(freqs_force, amp_force);
    hold on;
    xline(f_spindle, 'g--', 'spindle');
    xline(f_tooth, 'r--', 'tooth pass');
    xline(f_natural, 'm--', 'natural');
    xlabel('Frequency (Hz)');
    ylabel('Force FFT amplitude');
    title('Cutting force spectrum');
    xlim([0, min(max(3*f_tooth, 2*f_natural, 5000), fs_force/2)]);

    subplot(2, 1, 2);
    plot(freqs_defl, amp_defl);
    hold on;
    xline(f_spindle, 'g--', 'spindle');
    xline(f_tooth, 'r--', 'tooth pass');
    xline(f_natural, 'm--', 'natural');
    xlabel('Frequency (Hz)');
    ylabel('Deflection FFT amplitude');
    title('Tool deflection spectrum');
    xlim([0, min(max(3*f_tooth, 2*f_natural, 5000), fs_defl/2)]);
end
