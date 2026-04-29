% =========================================================================
% UFO_Buoy_Free_Decay_ULTIMATE.m
% =========================================================================
clear; clc; close all;

% 1. DATA INPUT (Capytaine Results)
w_cap = (0.2:0.302:15.0)';
A33_cap = [6.0445 6.0941 6.1580 6.2241 6.2837 6.3304 6.3591 6.3661 ...
           6.3492 6.3078 6.2423 6.1537 6.0447 5.9180 5.7766 5.6241 ...
           5.4637 5.2980 5.1306 4.9628 4.7973 4.6359 4.4792 4.3288 ...
           4.1854 4.0489 3.9202 3.7993 3.6861 3.5805 3.4826 3.3920 ...
           3.3086 3.2321 3.1622 3.0989 3.0417 2.9906 2.9452 2.9054 ...
           2.8712 2.8425 2.8194 2.8018 2.7900 2.7844 2.7856 2.7948 ...
           2.8133 2.8440]';
B33_cap = [0.0014 0.0226 0.0920 0.2358 0.4767 0.8323 1.3139 1.9264 ...
           2.6679 3.5299 4.4984 5.5550 6.6789 7.8476 9.0389 10.2314 ...
           11.4064 12.5468 13.6397 14.6737 15.6408 16.5354 17.3537 ...
           18.0949 18.7590 19.3458 19.8590 20.3009 20.6742 20.9830 ...
           21.2317 21.4242 21.5647 21.6575 21.7067 21.7170 21.6925 ...
           21.6376 21.5570 21.4554 21.3376 21.2082 21.0736 20.9404 ...
           20.8152 20.7096 20.6354 20.6085 20.6554 20.8122]';

% 2. PHYSICAL PARAMETERS
m = 1.4431;      
C33 = 619.1542;  
z0 = 0.030;      
rho = 1025;      
% Cd = 1.5;        
% Area = pi*0.14^2;

% 3. NATURAL FREQUENCY & MASS
% -------------------------------------------------------------------------
wn = 10;
for i = 1:20
    wn = sqrt(C33 / (m + interp1(w_cap, A33_cap, wn, 'pchip')));
end
Tn_theory = 2*pi/wn;

% A_inf is the infinite frequency limit. 3.5 is a realistic estimate
% given that A33 is trending down toward the end of your data.
A_inf = 3.5; 
M_total = m + A_inf;

% 4. RETARDATION KERNEL K(t) with TAIL EXTENSION
% -------------------------------------------------------------------------
dw = 0.01;
w_fine = (0:dw:60)'; % Integration up to 60 rad/s
B_fine = zeros(size(w_fine));

% Interpolate original data
idx_original = w_fine <= 15;
B_fine(idx_original) = interp1([0; w_cap], [0; B33_cap], w_fine(idx_original), 'pchip');

% Add exponential tail to bring B to 0 (prevents numerical ringing)
idx_tail = w_fine > 15;
B_end = B33_cap(end);
B_fine(idx_tail) = B_end * exp(-0.15 * (w_fine(idx_tail) - 15));

% Define K-time vector BEFORE the loop
dt_K = 0.002;
t_K = (0:dt_K:5)';
K = zeros(size(t_K));

% Inverse Cosine Transform
for i = 1:length(t_K)
    K(i) = (2/pi) * trapz(w_fine, B_fine .* cos(w_fine * t_K(i)));
end

% 5. TIME DOMAIN SIMULATION (Cummins Equation)
% -------------------------------------------------------------------------
dt = 0.002;
t = (0:dt:4)';
n = length(t);
z = zeros(n,1); v = zeros(n,1);
z(1) = z0;

fprintf('Running Simulation...\n');
for i = 1:n-1
    if i > 1
        % Determine how much history matches K
        n_hist = min(i, length(K));
        history = v(i:-1:i-n_hist+1);
        k_part  = K(1:n_hist);
        
        % Radiation force (Convolution)
        F_rad = dt * (sum(k_part .* history) - 0.5*k_part(1)*history(1) - 0.5*k_part(end)*history(end));
    else
        F_rad = 0;
    end
    
    % Viscous Force
    %F_visc = 0.5 * rho * Cd * Area * v(i) * abs(v(i));
  
    % Acceleration
    a = (-C33*z(i) - F_rad ) / M_total; %- F_visc
    
    % Integration
    v(i+1) = v(i) + dt*a;
    z(i+1) = z(i) + dt*v(i+1);
end

% 6. ANALYSIS & PLOTTING
% -------------------------------------------------------------------------
[pks, locs] = findpeaks(z, 'MinPeakDistance', round(0.25/dt), 'MinPeakHeight', 0.0005);
T_d = mean(diff(t(locs)));
log_dec = mean(log(pks(1:end-1) ./ pks(2:end)));
zeta = log_dec / sqrt(4*pi^2 + log_dec^2);

figure('Color','w','Position',[100 100 800 600]);
subplot(2,1,1);
plot(t, z*1000, 'b', 'LineWidth', 1.5); hold on;
plot(t(locs), pks*1000, 'ro'); grid on;
title(['Heave Free Decay: T_d = ', num2str(T_d,3), 's (Theory T_n = ', num2str(Tn_theory,3), 's)']);
ylabel('Heave (mm)'); xlabel('Time (s)');

subplot(2,1,2);
plot(t_K, K, 'r'); grid on; xlim([0 2]);
title('Retardation Function K(t) (Fixed Tail)');
ylabel('K(t)'); xlabel('Time (s)');

fprintf('--- FINAL CHECK ---\n');
fprintf('Theoretical Tn : %.4f s\n', Tn_theory);
fprintf('Simulated Td   : %.4f s\n', T_d);
fprintf('Damping Ratio  : %.2f %%\n', zeta*100);

% 7. ADDITIONAL PLOTS: Hydrodynamic Coefficients
% -------------------------------------------------------------------------
figure('Name','Hydrodynamic Coefficients','Color','w','Position',[150 150 800 600]);

% Plot Added Mass (A33)
subplot(2,1,1);
plot(w_cap, A33_cap, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
yline(A_inf, 'r--', ['Assumed A_\infty = ', num2str(A_inf)], 'LineWidth', 1.5);
grid on;
title('Heave Added Mass A_{33}(\omega)');
xlabel('Frequency \omega [rad/s]');
ylabel('Mass [kg]');
legend('Capytaine Data', 'A_{\infty} used in Sim');

% Plot Radiation Damping (B33) with the Fixed Tail
subplot(2,1,2);
plot(w_fine, B_fine, 'r', 'LineWidth', 2); hold on;
plot(w_cap, B33_cap, 'bo', 'MarkerSize', 4); 
grid on;
title('Heave Radiation Damping B_{33}(\omega) with Exponential Tail');
xlabel('Frequency \omega [rad/s]');
ylabel('Damping [N-s/m]');
legend('Interpolated + Tail', 'Original Capytaine Data');