%% plot_heave_data.m
% Reads and plots the heave (y) data for a heaving buoy from HeaveData2.txt
%
% NOTE ON FILE FORMAT:
%   The file uses a comma as BOTH the decimal separator AND the field
%   delimiter (European locale convention).  Each row therefore looks like:
%       0,799,1,897,2,108
%   which actually represents three values: t=0.799, x=1.897, y=2.108
%   The script handles this by reading every row as 6 raw tokens and then
%   reassembling the three columns (t, x, y) from adjacent token pairs.

clear; clc; close all;

%% ---- 1. File path --------------------------------------------------
filename = 'HeaveData2.txt';   % change path if the file lives elsewhere

%% ---- 2. Read raw data -----------------------------------------------
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end

% Skip the two header lines  ("mass A,," and "t,x,y,")
fgetl(fid);
fgetl(fid);

raw = [];
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line), continue; end

    % Split on comma  ->  6 tokens per data row
    tokens = strsplit(line, ',');
    if numel(tokens) < 6, continue; end

    % Reassemble decimal numbers: token pairs (integer + fractional part)
    vals = zeros(1, 3);
    for col = 1:3
        intPart  = str2double(tokens{2*col - 1});
        fracPart = str2double(tokens{2*col});
        % Determine the number of decimal digits in the fractional token
        fracStr  = strtrim(tokens{2*col});
        nDigits  = numel(fracStr);
        vals(col) = intPart + fracPart / (10^nDigits);
    end
    raw(end+1, :) = vals; %#ok<AGROW>
end
fclose(fid);

%% ---- 3. Extract columns --------------------------------------------
t = raw(:, 1);   % time  [s]
% x = raw(:, 2); % surge / horizontal – not plotted here
y = raw(:, 3);   % heave [mm or units as recorded]

%% ---- 4. Period calculation via local peak detection ----------------
% The signal has a large transient at the start that biases the mean,
% making zero-crossing methods unreliable.  Instead we detect ALL local
% peaks (findpeaks) with a minimum prominence to avoid noise spikes,
% then compute peak-to-peak intervals.
%
% We skip the first transient peak (the large launch spike) and work on
% the subsequent oscillating tail which represents the natural period.

% Minimum peak prominence: 20% of overall range to filter noise
minProm = 0.20 * (max(y) - min(y));
% Minimum separation between peaks: ~10 samples (avoids double-detection)
minSep  = 10;

[pkVals, pkIdx] = findpeaks(y, 'MinPeakProminence', minProm, ...
                                'MinPeakDistance',   minSep);
t_peaks = t(pkIdx);

% Skip the first (transient) peak — it is the large launch spike
if numel(pkVals) > 1
    pkVals_ss  = pkVals(2:end);     % steady-state peaks
    pkIdx_ss   = pkIdx(2:end);
    t_peaks_ss = t_peaks(2:end);
else
    pkVals_ss  = pkVals;
    pkIdx_ss   = pkIdx;
    t_peaks_ss = t_peaks;
end

% Period from successive steady-state peaks
if numel(t_peaks_ss) >= 2
    periods   = diff(t_peaks_ss);
    T_mean    = mean(periods);
    T_std     = std(periods);
    freq_mean = 1 / T_mean;
else
    warning('Not enough peaks found to estimate period.');
    periods = NaN; T_mean = NaN; T_std = NaN; freq_mean = NaN;
end

%% ---- 5. Plot --------------------------------------------------------
figure('Name', 'Heaving Buoy – Heave Response', 'NumberTitle', 'off', ...
       'Color', 'white', 'Units', 'normalized', 'Position', [0.1 0.2 0.75 0.55]);

plot(t, y, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'b');

xlabel('Time  (s)',   'FontSize', 13);
ylabel('Heave  (mm)', 'FontSize', 13);
title('Heaving Buoy – Heave (y) vs Time', 'FontSize', 15, 'FontWeight', 'bold');
grid on;
grid minor;
hold on;

% Mark the transient (first) peak in red
[peakVal, peakIdx] = max(y);
plot(t(peakIdx), peakVal, 'rv', 'MarkerSize', 10, 'LineWidth', 2, ...
     'DisplayName', 'Transient peak');
text(t(peakIdx) + 0.05, peakVal, ...
     sprintf('Peak = %.2f mm\n@ t = %.3f s', peakVal, t(peakIdx)), ...
     'FontSize', 10, 'Color', 'r', 'VerticalAlignment', 'bottom');

% Mark steady-state peaks used for period estimation in green
plot(t_peaks_ss, pkVals_ss, 'g^', 'MarkerSize', 8, 'LineWidth', 1.5, ...
     'DisplayName', 'SS peaks (period)');

% Period annotation box
if ~isnan(T_mean)
    annotation('textbox', [0.13 0.72 0.22 0.16], ...
        'String', {sprintf('Mean period:  T = %.3f s', T_mean), ...
                   sprintf('Std of period:    %.3f s', T_std), ...
                   sprintf('Frequency:    f = %.3f Hz', freq_mean)}, ...
        'FitBoxToText', 'on', 'BackgroundColor', [0.95 0.95 0.8], ...
        'EdgeColor', [0.4 0.4 0], 'FontSize', 10);
end

hold off;
legend({'Heave  y(t)', 'Transient peak', 'SS peaks (period)'}, ...
       'Location', 'northeast', 'FontSize', 11);
xlim([t(1) t(end)]);

%% ---- 6. Quick summary to Command Window ----------------------------
fprintf('\n--- Heave Data Summary (mass A) ---\n');
fprintf('  Duration   : %.3f s  (%d samples)\n', t(end)-t(1), numel(t));
fprintf('  Min heave  : %.2f mm  @ t = %.3f s\n', min(y), t(find(y == min(y), 1)));
fprintf('  Max heave  : %.2f mm  @ t = %.3f s\n', peakVal, t(peakIdx));
fprintf('  Mean heave : %.2f mm\n', mean(y));
fprintf('  Std  heave : %.2f mm\n', std(y));
fprintf('\n--- Period Estimation (peak-to-peak, steady-state) ---\n');
fprintf('  SS peaks found       : %d\n', numel(t_peaks_ss));
fprintf('  Number of periods    : %d\n', numel(periods));
if ~isnan(T_mean)
    fprintf('  Individual periods   : ');
    fprintf('%.3f  ', periods);
    fprintf('s\n');
    fprintf('  Mean period  T       : %.4f s\n', T_mean);
    fprintf('  Std  period          : %.4f s\n', T_std);
    fprintf('  Frequency    f       : %.4f Hz\n\n', freq_mean);
end