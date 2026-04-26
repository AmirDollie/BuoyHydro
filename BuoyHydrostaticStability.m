% =========================================================================
% BUOY HYDROSTATIC STABILITY ANALYSIS
% Based on Fossen's GM_surfaced2submerged (MSS Toolbox)
%
% Buoy geometry: Two hemispherical bowls + polyethylene stability ring
%   - Total mass:      1406.85 g  (Physical Mass from new SW props)
%   - Total volume:    1,406,850.00 mm^3 (from Buoyancy Model)
%   - Total Height:    200 mm
%   - Hull diameter:   200 mm
%   - Ring diameter:   280 mm  (waterplane)
%   - KG (from keel):  88.08 mm  (SolidWorks CoM Y = 88.08 mm)
%   - KB (from keel):  56.40 mm  (Buoyancy Model CoM Y = 56.40 mm)
%
% SolidWorks Inertia (at CoM, principal axes) [g·mm^2]:
%   (Mapped: SW Lxx -> Marine Ixx, SW Lzz -> Marine Iyy, SW Lyy -> Marine Izz)
%   Ixx_cm = 6,653,303.96   Iyy_cm = 6,653,158.46   Izz_cm = 5,828,458.75
% SolidWorks Inertia (at keel/output coord origin) [g·mm^2]:
%   Ixx_k  = 17,571,343.43  Iyy_k  = 17,569,729.07  Izz_k  = 5,828,458.75
%
% Sign convention (Fossen): z positive DOWNWARD from waterline
%   z_from_WL = T - z_from_keel
%
% Toolbox: MSS (Fossen 2021) - ensure GM_surfaced2submerged is on your path
% =========================================================================
close all; clear; clc;

% -------------------------------------------------------------------------
% 1. GEOMETRY & MASS PROPERTIES  (all SI units: metres, kg)
% -------------------------------------------------------------------------
D_ring  = 0.280;            % Polyethylene ring / waterplane diameter [m]
D_hull  = 0.200;            % Hull outer diameter [m]
T       = 0.0866;           % Draft at equilibrium / height of buoyancy model [m]

% --- SolidWorks mass (as-designed, pre-ballast) ---
m_SW    = 1.40685;          % Physical Mass from SolidWorks [kg]
rho_sw  = 1025;             % Seawater density [kg/m^3]
g_acc   = 9.81;             % Gravity [m/s^2]

% Displaced volume from the new buoyancy model
nabla   = 1406850.00 * 1e-9; % Displaced volume at equilibrium [m^3]

m       = nabla * rho_sw;    % Total operational mass required [kg] (~1.442 kg)
delta_m = m - m_SW;          % Ballast required to hit draft [kg] (~0.035 kg)

% --- SolidWorks moments of inertia (convert g·mm^2 → kg·m^2) ---
%     At centre of mass (principal axes) - CORRECTED FOR Y/Z SWAP
I_xx_cm = 6653303.96 * 1e-9;   % [kg·m^2]  roll  axis (SW Lxx)
I_yy_cm = 6653158.46 * 1e-9;   % [kg·m^2]  pitch axis (SW Lzz)
I_zz_cm = 5828458.75 * 1e-9;   % [kg·m^2]  yaw   axis (SW Lyy)

% -------------------------------------------------------------------------
% 2. HYDROSTATIC PARAMETERS AT EQUILIBRIUM (from keel)
% -------------------------------------------------------------------------
KG = 0.08808;               % Centre of gravity above keel [m]
KB = 0.05640;               % Centre of buoyancy above keel [m]

% Convert to Fossen sign convention: z positive DOWN from waterline
z_g           = T - KG;     % CoG below waterline (positive down)
z_b_surface   = T - KB;     % CoB below waterline at surface condition
z_b_submerged = T / 2;      % CoB when fully submerged

% -------------------------------------------------------------------------
% 3. WATERPLANE MOMENT OF INERTIA
% -------------------------------------------------------------------------
I_T_ring = (pi * D_ring^4) / 64;   % [m^4]  Circular waterplane (ring)
I_T_hull = (pi * D_hull^4) / 64;   % [m^4]  Hull-only (for comparison)

fprintf('===== BUOY HYDROSTATIC SUMMARY =====\n');
fprintf('Draft T            = %.1f mm\n',  T*1e3);
fprintf('Displaced volume   = %.6f m^3\n', nabla);
fprintf('Displaced mass     = %.3f kg\n',  nabla*rho_sw);
fprintf('KG                 = %.2f mm\n',  KG*1e3);
fprintf('KB                 = %.2f mm\n',  KB*1e3);

fprintf('\n--- At equilibrium ---\n');
BM  = I_T_ring / nabla;
GM  = KB - KG + BM;

fprintf('I_T (ring Ø280mm)  = %.4e m^4\n', I_T_ring);
fprintf('BM                 = %.2f mm\n',  BM*1e3);
fprintf('GM_T               = %.2f mm  (%s)\n', GM*1e3, ...
    ternary(GM > 0, 'STABLE ✓', 'UNSTABLE ✗'));
fprintf('BG                 = %.2f mm\n', (KG-KB)*1e3);
fprintf('=====================================\n\n');

% -------------------------------------------------------------------------
% 4. SWEEP: GM vs submergence depth (surface → fully submerged)
% -------------------------------------------------------------------------
targetDepth = 0.4;                  % Plot to 400 mm depth
zn_values   = 0 : 0.002 : targetDepth;
N           = length(zn_values);
GM_T_ring   = zeros(1, N);
GM_T_hull   = zeros(1, N);
z_b_vals    = zeros(1, N);
BG_z_vals   = zeros(1, N);

for i = 1:N
    zn = zn_values(i);
    % Using ring waterplane (actual configuration)
    [GM_T_ring(i), ~, z_b_vals(i)] = GM_surfaced2submerged( ...
        I_T_ring, nabla, zn, T, z_b_surface, z_b_submerged, z_g);
    
    % Using hull-only waterplane (for comparison — no ring)
    [GM_T_hull(i), ~, ~] = GM_surfaced2submerged( ...
        I_T_hull, nabla, zn, T, z_b_surface, z_b_submerged, z_g);
    
    BG_z_vals(i) = z_g - z_b_vals(i);
end

% Find exact zn where GM crosses zero (neutral stability / capsize point)
idx_GM0 = find(GM_T_ring <= 0, 1);
if ~isempty(idx_GM0) && idx_GM0 > 1
    zn_GM_zero = interp1(GM_T_ring(idx_GM0-1:idx_GM0), zn_values(idx_GM0-1:idx_GM0), 0);
else
    zn_GM_zero = NaN; % If it never capsizes in the sweep range
end

% -------------------------------------------------------------------------
% 5. PLOT
% -------------------------------------------------------------------------

% 5. PLOT
% -------------------------------------------------------------------------
fig = figure('Name', 'Buoy Hydrostatic Stability', ...
             'Color', [0.97 0.97 0.97], ...
             'Position', [100 100 900 600]);

subplot(2,1,1);
hold on; grid on; box on;
plot(zn_values*1e3, GM_T_ring*1e3,  'r-',  'LineWidth', 2.5, ...
    'DisplayName', 'GM_T  — with ring (Ø280mm)');
plot(zn_values*1e3, GM_T_hull*1e3,  'b--', 'LineWidth', 1.5, ...
    'DisplayName', 'GM_T  — hull only (Ø200mm)');

yline(0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'GM = 0 (Neutral)');

% --- Clean Capsize Line + Custom Label ---
if ~isnan(zn_GM_zero)
    x_cap = zn_GM_zero * 1e3;

    % Vertical line (no built-in label)
    xline(x_cap, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Custom label (clean, readable placement)
    y_offset = 10; % mm above GM = 0 line

    text(x_cap + 5, y_offset, ...
        sprintf('Capsize (%.1f mm)', x_cap), ...
        'FontSize', 10, ...
        'Color', 'k', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.8 0.8 0.8], ...
        'Margin', 4, ...
        'VerticalAlignment', 'bottom');
end

ylabel('Height [mm]');
title('Transverse Metacentric Height GM_T vs. Submergence Depth');
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11);

subplot(2,1,2);
hold on; grid on; box on;
plot(zn_values*1e3, BG_z_vals*1e3, 'g-', 'LineWidth', 2.0, ...
    'DisplayName', 'BG_z = z_g - z_b  (lever arm)');
plot(zn_values*1e3, z_b_vals*1e3,  'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'z_b (CoB below WL)');
yline(z_g*1e3, 'c--', 'LineWidth', 1.5, ...
    'Label', sprintf('z_g = %.1f mm', z_g*1e3), ...
    'LabelHorizontalAlignment', 'right', ...
    'DisplayName', 'z_g (CoG below WL, constant)');

% Matching vertical line (no label)
if ~isnan(zn_GM_zero)
    xline(zn_GM_zero*1e3, 'k:', 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end

xlabel('Submergence depth z_n  [mm] (Additional depth beyond equilibrium)');
ylabel('Distance from waterline [mm]  (+ve down)');
title('Centre of Buoyancy and Gravity vs. Submergence Depth');
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11);

sgtitle('Buoy Hydrostatic Stability — Fossen MSS Toolbox', ...
    'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% 6. STIFFNESS MATRIX ELEMENT C33, C44, C55 (heave & roll restoring)
% -------------------------------------------------------------------------
fprintf('\n===== LINEARISED RESTORING COEFFICIENTS =====\n');
A_WP = pi/4 * D_ring^2;        % Waterplane area [m^2]
C33  = rho_sw * g_acc * A_WP;  % Heave restoring [N/m]
C44  = rho_sw * g_acc * nabla * GM;  % Roll restoring [N·m/rad]
C55  = C44;                     % Pitch = Roll (axisymmetric)

fprintf('Waterplane area A_WP = %.4f m^2\n', A_WP);
fprintf('C33 (heave stiffness) = %.2f N/m\n',     C33);
fprintf('C44 (roll  stiffness) = %.4f N·m/rad\n', C44);
fprintf('C55 (pitch stiffness) = %.4f N·m/rad\n', C55);

% Heave natural period (undamped, added mass ignored as first estimate)
T_heave = 2*pi * sqrt(m / C33);
fprintf('\nUndamped heave natural period T_n3 ≈ %.3f s\n', T_heave);

% Roll & Pitch natural periods — using EXACT SolidWorks I_xx at CoM
I_xx_total = I_xx_cm + delta_m * KG^2;   % [kg·m^2]  total roll inertia at CoM
I_yy_total = I_yy_cm + delta_m * KG^2;   % [kg·m^2]  pitch (same by symmetry)

T_roll  = 2*pi * sqrt(I_xx_total / C44);
T_pitch = 2*pi * sqrt(I_yy_total / C55);
T_yaw   = 2*pi * sqrt(I_zz_cm   / 1e-6); % yaw: no hydrostatic restoring, placeholder

fprintf('Roll  inertia I_xx (at CoM, incl. ballast) = %.4e kg·m^2\n', I_xx_total);
fprintf('Pitch inertia I_yy (at CoM, incl. ballast) = %.4e kg·m^2\n', I_yy_total);
fprintf('Yaw   inertia I_zz (spin, at CoM)          = %.4e kg·m^2\n', I_zz_cm);

fprintf('\nUndamped roll  natural period T_n4 ≈ %.3f s  (exact SolidWorks I_xx)\n', T_roll);
fprintf('Undamped pitch natural period T_n5 ≈ %.3f s  (exact SolidWorks I_yy)\n', T_pitch);
fprintf('(Yaw: no hydrostatic restoring — mooring dependent)\n');
fprintf('=============================================\n');

% =========================================================================
% Helper function (place at end of script for MATLAB R2016b+ / Octave)
% =========================================================================
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end