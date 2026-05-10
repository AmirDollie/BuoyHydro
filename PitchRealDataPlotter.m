%% plot_pitch_data.m
% Reads pitch angle (θ), shifts release to t = 0, and truncates at specified end point.
clear; clc; close all;

%% ---- 1. File path --------------------------------------------------
filename = 'ActualPitchData.txt';   

%% ---- 2. Read raw data -----------------------------------------------
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Skip the two header lines
fgetl(fid);
fgetl(fid);

t_raw = [];
theta_raw = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line), continue; end
    
    % Handle European decimal commas and scientific notation
    line = strrep(line, ',', '.');
    tokens = strsplit(line, '\t');
    if numel(tokens) < 2, continue; end
    
    t_raw(end+1,1)     = str2double(tokens{1});      
    theta_raw(end+1,1) = str2double(tokens{2});  
end
fclose(fid);

%% ---- 3. Shift to t = 0 and Truncate at specified End Point ---------
% A. Find release point (where it stops being flat)
flatThresh = 0.5; 
startIdx = find(abs(theta_raw - theta_raw(1)) > flatThresh, 1, 'first');
if isempty(startIdx), startIdx = 1; end

% B. Shift time relative to release
t_shifted = t_raw(startIdx:end) - t_raw(startIdx);
theta_shifted = theta_raw(startIdx:end);

% C. Truncate at the specified X coordinate from Tracker
targetEnd = 3.43031; 
[~, endIdx] = min(abs(t_shifted - targetEnd));

t_final = t_shifted(1:endIdx);
theta_final = theta_shifted(1:endIdx);

%% ---- 4. Period calculation (Command Window Only) --------------------
minProm = 0.08 * (max(theta_final) - min(theta_final));
[pkVals, pkIdx] = findpeaks(theta_final, 'MinPeakProminence', minProm, 'MinPeakDistance', 5);

if numel(pkIdx) > 1
    T_mean = mean(diff(t_final(pkIdx(2:end)))); % Skipping transient peak
else
    T_mean = NaN;
end

%% ---- 5. Plot --------------------------------------------------------
figure('Name', 'Heaving Buoy – Pitch Angle', 'NumberTitle', 'off', ...
       'Color', 'white', 'Units', 'normalized', 'Position', [0.1 0.2 0.75 0.55]);

% Main data plot
plot(t_final, theta_final, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
     'MarkerFaceColor', 'r', 'DisplayName', 'Pitch \theta(t)'); 
hold on;

% Zero reference line
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'Zero line');

% Steady-State Peaks (to show where period is measured)
if ~isnan(T_mean)
    plot(t_final(pkIdx(2:end)), pkVals(2:end), 'g^', 'MarkerSize', 8, ...
         'LineWidth', 1.5, 'DisplayName', 'SS peaks (period)');
end

% Formatting
xlabel('Time  (s)', 'FontSize', 13);
ylabel('Pitch Angle  (°)', 'FontSize', 13);
title('Heaving Buoy – Pitch Angle (\theta) vs Time (s)', ...
      'FontSize', 15, 'FontWeight', 'bold');
grid on; grid minor;

% Clean Legend (Removed: crop start, flatline, transient peak)
legend('show', 'Location', 'northeast', 'FontSize', 10);

hold off;
xlim([0 t_final(end)]);

%% ---- 6. Summary to Command Window ----------------------------------
fprintf('\n--- Pitch Angle Analysis (t=0 at release) ---\n');
fprintf('  Final Data Point : X = %.5f, Y = %.5f\n', t_final(end), theta_final(end));
if ~isnan(T_mean)
    fprintf('  Measured Period T : %.4f s\n', T_mean);
    fprintf('  Measured Freq   f : %.4f Hz\n', 1/T_mean);
end