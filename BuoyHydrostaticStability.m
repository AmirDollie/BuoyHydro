% =========================================================================
% BUOY HYDROSTATIC STABILITY ANALYSIS
% Based on Fossen's GM_surfaced2submerged (MSS Toolbox)
%
% Buoy geometry: Two hemispherical bowls + polyethylene stability ring
%   - Total mass:    1443 g  (incl. 127.5 g ballast for 50% submergence)
%   - Total volume:  3,064,433 mm^3
%   - Hull diameter: 200 mm
%   - Ring diameter: 280 mm  (waterplane)
%   - KG (from keel): 89.71 mm
%   - KB (from keel): 72.19 mm
%
% Sign convention (Fossen): z positive DOWNWARD from waterline
%   z_from_WL = T - z_from_keel
%
% Author:  [Your name]
% Toolbox: MSS (Fossen 2021) - ensure GM_surfaced2submerged is on your path
% =========================================================================
close all; clear; clc;

% -------------------------------------------------------------------------
% 1. GEOMETRY & MASS PROPERTIES  (all SI units: metres, kg)
% -------------------------------------------------------------------------
D_ring  = 0.280;            % Polyethylene ring / waterplane diameter [m]
D_hull  = 0.200;            % Hull outer diameter [m]
T       = 0.1005;           % Draft at 50% submergence [m]
%   Total height ~ 201 mm, 50% => T = 100.5 mm

m       = 1.5705;           % Total mass incl. 127.5 g ballast [kg]
rho_sw  = 1025;             % Seawater density [kg/m^3]
g_acc   = 9.81;             % Gravity [m/s^2]

nabla   = m / rho_sw;       % Displaced volume at equilibrium [m^3]
%   = 1.5705 / 1025 = 0.001532 m^3  (matches V_total/2)

% -------------------------------------------------------------------------
% 2. HYDROSTATIC PARAMETERS AT EQUILIBRIUM (from keel)
% -------------------------------------------------------------------------
KG = 0.08971;               % Centre of gravity above keel [m]
KB = 0.07219;               % Centre of buoyancy above keel [m]

% Convert to Fossen sign convention: z positive DOWN from waterline
% z = T - z_keel
z_g           = T - KG;     % CoG below waterline (positive down)
z_b_surface   = T - KB;     % CoB below waterline at surface condition

% When fully submerged the CoB rises to the geometric centroid of the body.
% Approximating as a sphere/oblate body: z_b_submerged ~ T/2 from waterline
z_b_submerged = T / 2;

% -------------------------------------------------------------------------
% 3. WATERPLANE MOMENT OF INERTIA
%    Using the RING diameter (280 mm) — this is what sets BM
% -------------------------------------------------------------------------
I_T_ring = (pi * D_ring^4) / 64;   % [m^4]  Circular waterplane (ring)
I_T_hull = (pi * D_hull^4) / 64;   % [m^4]  Hull-only (for comparison)

fprintf('===== BUOY HYDROSTATIC SUMMARY =====\n');
fprintf('Draft T            = %.1f mm\n',  T*1e3);
fprintf('Displaced volume   = %.6f m^3\n', nabla);
fprintf('Displaced mass     = %.3f kg\n',  nabla*rho_sw);
fprintf('KG                 = %.2f mm\n',  KG*1e3);
fprintf('KB                 = %.2f mm\n',  KB*1e3);
fprintf('\n--- At equilibrium (50%% submerged) ---\n');
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
BM_T_ring   = zeros(1, N);
z_b_vals    = zeros(1, N);
BG_z_vals   = zeros(1, N);

for i = 1:N
    zn = zn_values(i);

    % Using ring waterplane (actual configuration)
    [GM_T_ring(i), BM_T_ring(i), z_b_vals(i)] = GM_surfaced2submerged( ...
        I_T_ring, nabla, zn, T, z_b_surface, z_b_submerged, z_g);

    % Using hull-only waterplane (for comparison — no ring)
    [GM_T_hull(i), ~, ~] = GM_surfaced2submerged( ...
        I_T_hull, nabla, zn, T, z_b_surface, z_b_submerged, z_g);

    BG_z_vals(i) = z_g - z_b_vals(i);
end

% -------------------------------------------------------------------------
% 5. PLOT
% -------------------------------------------------------------------------
fig = figure('Name', 'Buoy Hydrostatic Stability', ...
             'Color', [0.97 0.97 0.97], ...
             'Position', [100 100 900 600]);

% Mark equilibrium draft
xline_val = T;  % equilibrium draft

subplot(2,1,1);
hold on; grid on; box on;
plot(zn_values*1e3, GM_T_ring*1e3,  'r-',  'LineWidth', 2.5, ...
    'DisplayName', 'GM_T  — with ring (Ø280mm)');
plot(zn_values*1e3, GM_T_hull*1e3,  'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'GM_T  — hull only (Ø200mm)');
plot(zn_values*1e3, BM_T_ring*1e3,  'b-',  'LineWidth', 2.0, ...
    'DisplayName', 'BM_T  — with ring');
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'GM = 0 (neutral)');
xline(T*1e3, 'm:', 'LineWidth', 2, ...
    'Label', sprintf('50%% submergence\nT=%.1fmm', T*1e3), ...
    'LabelHorizontalAlignment', 'left');
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
xline(T*1e3, 'm:', 'LineWidth', 2, ...
    'Label', sprintf('50%% subm. T=%.1fmm', T*1e3), ...
    'LabelHorizontalAlignment', 'left');
xlabel('Submergence depth z_n  [mm]');
ylabel('Distance from waterline [mm]  (+ve down)');
title('Centre of Buoyancy and Gravity vs. Submergence Depth');
legend('Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11);

sgtitle('Buoy Hydrostatic Stability — Fossen MSS Toolbox', ...
    'FontSize', 13, 'FontWeight', 'bold');
% (per-line widths already set above in each plot call)

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

% Roll natural period (rough estimate, need I_xx)
% I_xx approx for thin disk: m*R^2/4
R_hull = D_hull/2;
I_xx_est = m * R_hull^2 / 4;
T_roll = 2*pi * sqrt(I_xx_est / C44);
fprintf('Undamped roll  natural period T_n4 ≈ %.3f s  (rough, needs CAD I_xx)\n', T_roll);
fprintf('=============================================\n');

% =========================================================================
% Helper function (place at end of script for MATLAB R2016b+ / Octave)
% =========================================================================
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

