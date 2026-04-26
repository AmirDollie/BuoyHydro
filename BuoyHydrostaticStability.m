% =========================================================================
% BUOY HYDROSTATIC & DYNAMIC MATRIX ANALYSIS
% Combines: Hydrostatic sweep (Fossen GM_surfaced2submerged) +
%           Rigid-body matrices (rbody, Gmtrx)
%
% Buoy geometry: Two hemispherical bowls + polyethylene stability ring
%   - Total Height:    200 mm  |  Hull Ø: 200 mm  |  Ring Ø: 280 mm
%   - Draft (T):        86.6 mm  (buoyancy model height)
%
% SOURCE — SolidWorks Mass Properties:
%   Assembly3_MoreSpace (full buoy, Y_upCoords):
%     Mass              = 1443.14 g
%     Volume            = 3,064,876.89 mm^3
%     KG (Y of CoM)     = 89.71 mm above keel
%     Ixx_cm = 6,901,795.53   Iyy_cm = 6,905,469.92   Izz_cm = 9,784,262.59  [g·mm^2]
%     Ixx_k  = 18,520,213.57  Iyy_k  = 9,784,262.59   Izz_k  = 18,516,539.19 [g·mm^2]
%
%   BuoyancyModelPETG (X_upCoordskeel):
%     Volume (nabla)          = 1,464,525.51 mm^3
%     KB (X of CoM from keel) = 56.52 mm
%     Mass (rho_sw * nabla)   = 1464.53 g  → buoy floats with +0.569 N reserve
%
% Toolbox: MSS (Fossen 2021) — rbody, Gmtrx, GM_surfaced2submerged
% =========================================================================
close all; clear; clc;

% =========================================================================
% SECTION 1 — GEOMETRY & MASS
% =========================================================================
D_ring  = 0.280;               % Ring / waterplane diameter      [m]
D_hull  = 0.200;               % Hull outer diameter             [m]
T       = 0.0866;              % Equilibrium draft               [m]

m       = 1.44314;             % As-built assembly mass (SW)     [kg]
rho_sw  = 1025;                % Seawater density                [kg/m^3]
g_acc   = 9.81;                % Gravity                         [m/s^2]

nabla   = 1464525.51e-9;       % Displaced volume (buoyancy mdl) [m^3]

% =========================================================================
% SECTION 2 — MOMENTS OF INERTIA  (g·mm^2 → kg·m^2)
%   SW Y_upCoords axis mapping → Marine body-fixed:
%     SW Px (Lxx) → Ixx (roll)
%     SW Py (Lyy) → Iyy (pitch)   [NOTE: Pz=Izz_cm is the spin/yaw axis]
%     SW Pz (Lzz) → Izz (yaw)
% =========================================================================
I_xx_cm = 6901795.53e-9;       % Roll  axis  [kg·m^2]
I_yy_cm = 6905469.92e-9;       % Pitch axis  [kg·m^2]
I_zz_cm = 9784262.59e-9;       % Yaw   axis  [kg·m^2]

% At keel (for parallel-axis verification only)
I_xx_k  = 18520213.57e-9;
I_yy_k  = 9784262.59e-9;
I_zz_k  = 18516539.19e-9;

KG      = 0.08971;             % CoG above keel  [m]

% Radii of gyration at CoM (no ballast — use as-built mass directly)
k44 = sqrt(I_xx_cm / m);      % Roll  [m]
k55 = sqrt(I_yy_cm / m);      % Pitch [m]
k66 = sqrt(I_zz_cm / m);      % Yaw   [m]

% =========================================================================
% SECTION 3 — HYDROSTATICS AT EQUILIBRIUM
% =========================================================================
KB      = 0.05652;             % CoB above keel  [m]

I_T_ring = (pi * D_ring^4) / 64;   % 2nd moment of waterplane area (ring) [m^4]
I_T_hull = (pi * D_hull^4) / 64;   % Hull-only (comparison)               [m^4]
A_WP     = pi/4 * D_ring^2;        % Waterplane area                      [m^2]

BM  = I_T_ring / nabla;            % Metacentric radius    [m]
GM  = KB - KG + BM;                % Transverse GM         [m]
GML = GM;                          % = GM_T (axisymmetric)

% Fossen sign convention: z positive DOWN from waterline
z_g           = T - KG;
z_b_surface   = T - KB;
z_b_submerged = T / 2;

% =========================================================================
% SECTION 4 — PRINT HYDROSTATIC SUMMARY
% =========================================================================
fprintf('================================================================\n');
fprintf('  BUOY HYDROSTATIC SUMMARY\n');
fprintf('================================================================\n');
fprintf('  Assembly mass m           = %8.4f kg  (%7.2f g)\n', m, m*1e3);
fprintf('  Displaced volume nabla    = %8.6f m^3\n', nabla);
fprintf('  Buoyancy force            = %8.4f N\n', rho_sw*g_acc*nabla);
fprintf('  Weight                    = %8.4f N\n', m*g_acc);
fprintf('  Net upward force          = %+8.4f N  (%s)\n', ...
    (rho_sw*nabla - m)*g_acc, ...
    ternary((rho_sw*nabla - m) > 0, 'floats ✓', 'sinks ✗'));
fprintf('  Draft T                   = %8.2f mm\n',  T*1e3);
fprintf('  KG                        = %8.3f mm\n',  KG*1e3);
fprintf('  KB                        = %8.3f mm\n',  KB*1e3);
fprintf('  BG  (= KG - KB)           = %8.3f mm\n',  (KG-KB)*1e3);
fprintf('  BM  (= I_T / nabla)       = %8.3f mm\n',  BM*1e3);
fprintf('  GM_T                      = %8.3f mm  [%s]\n', GM*1e3, ...
    ternary(GM > 0, 'STABLE ✓', 'UNSTABLE ✗'));
fprintf('  A_WP (ring Ø280 mm)       = %8.6f m^2\n', A_WP);
fprintf('  I_T  (ring Ø280 mm)       = %8.4e m^4\n', I_T_ring);
fprintf('----------------------------------------------------------------\n');
fprintf('  Radii of gyration (at CoM, as-built SW inertias):\n');
fprintf('  k44 (roll)  = %7.4f m  (%6.3f mm)\n', k44, k44*1e3);
fprintf('  k55 (pitch) = %7.4f m  (%6.3f mm)\n', k55, k55*1e3);
fprintf('  k66 (yaw)   = %7.4f m  (%6.3f mm)\n', k66, k66*1e3);
fprintf('================================================================\n\n');

% Parallel-axis verification
fprintf('  Parallel-axis check (I_cm + m*KG^2 should match SW keel value):\n');
I_xx_PA = I_xx_cm + m * KG^2;
fprintf('  Computed I_xx_keel  = %12.2f g·mm^2\n', I_xx_PA*1e9);
fprintf('  SolidWorks I_xx_k   = %12.2f g·mm^2\n', I_xx_k *1e9);
fprintf('  Difference          = %+.2f g·mm^2\n\n', (I_xx_PA - I_xx_k)*1e9);

% =========================================================================
% SECTION 5 — LINEARISED RESTORING COEFFICIENTS & NATURAL PERIODS
% =========================================================================
C33 = rho_sw * g_acc * A_WP;          % Heave stiffness    [N/m]
C44 = rho_sw * g_acc * nabla * GM;    % Roll  stiffness    [N·m/rad]
C55 = C44;                            % Pitch (axisymmetric)

T_heave = 2*pi * sqrt(m       / C33);
T_roll  = 2*pi * sqrt(I_xx_cm / C44);
T_pitch = 2*pi * sqrt(I_yy_cm / C55);

fprintf('================================================================\n');
fprintf('  LINEARISED RESTORING COEFFICIENTS\n');
fprintf('================================================================\n');
fprintf('  C33 (heave) = %10.4f N/m\n',      C33);
fprintf('  C44 (roll)  = %10.6f N·m/rad\n',  C44);
fprintf('  C55 (pitch) = %10.6f N·m/rad\n',  C55);
fprintf('  C66 (yaw)   =   0  (mooring dependent)\n');
fprintf('----------------------------------------------------------------\n');
fprintf('  Undamped natural periods (no added mass — lower bound):\n');
fprintf('  T_n3 heave  = %.4f s\n', T_heave);
fprintf('  T_n4 roll   = %.4f s\n', T_roll);
fprintf('  T_n5 pitch  = %.4f s\n', T_pitch);
fprintf('  T_n6 yaw    = N/A  (no hydrostatic restoring)\n');
fprintf('================================================================\n\n');

% =========================================================================
% SECTION 6 — RIGID-BODY & HYDROSTATIC MATRICES  (rbody + Gmtrx)
% =========================================================================
% rbody: r_G from CO (keel), z-up body-fixed → CoG above keel = +KG
r_G   = [0; 0; KG];
r_bP  = [0; 0; 0];            % Reference point = keel
nu2   = [0; 0; 0];            % At rest
x_F   = 0;                    % Centre of flotation (axisymmetric)

[MRB, CRB] = rbody(m, k44, k55, k66, nu2, r_G);
G = Gmtrx(nabla, A_WP, GM, GML, x_F, r_bP);

% =========================================================================
% SECTION 7 — PRINT MATRICES
% =========================================================================
dof = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};

fprintf('================================================================\n');
fprintf('  BUOY SYSTEM MATRICES  (Reference point: Keel)\n');
fprintf('================================================================\n\n');

fprintf('  1. Hydrostatic Stiffness Matrix G  [N/m, N·m/rad]\n');
fprintf('     (non-zero: G33=heave, G44=roll, G55=pitch)\n');
disp(array2table(round(G,4), 'VariableNames', dof, 'RowNames', dof));

fprintf('  2. Rigid-Body Mass Matrix MRB  [kg, kg·m^2]\n');
disp(array2table(round(MRB,5), 'VariableNames', dof, 'RowNames', dof));

fprintf('  3. Coriolis Matrix CRB at rest (nu=[0,0,0]) — should be zero\n');
disp(CRB);
fprintf('================================================================\n\n');

% =========================================================================
% SECTION 8 — HYDROSTATIC DEPTH SWEEP (GM vs submergence)
% =========================================================================
targetDepth = 0.40;
zn_values   = 0 : 0.001 : targetDepth;
N           = length(zn_values);

GM_T_ring = zeros(1,N);  GM_T_hull = zeros(1,N);
BM_T_ring = zeros(1,N);
z_b_vals  = zeros(1,N);  BG_z_vals = zeros(1,N);

for i = 1:N
    zn = zn_values(i);
    [GM_T_ring(i), BM_T_ring(i), z_b_vals(i)] = GM_surfaced2submerged( ...
        I_T_ring, nabla, zn, T, z_b_surface, z_b_submerged, z_g);
    [GM_T_hull(i), ~, ~] = GM_surfaced2submerged( ...
        I_T_hull, nabla, zn, T, z_b_surface, z_b_submerged, z_g);
    BG_z_vals(i) = z_g - z_b_vals(i);
end

% Capsize depth (GM = 0 crossing)
idx_GM0 = find(GM_T_ring <= 0, 1);
if ~isempty(idx_GM0) && idx_GM0 > 1
    zn_GM_zero = interp1(GM_T_ring(idx_GM0-1:idx_GM0), ...
                         zn_values(idx_GM0-1:idx_GM0), 0);
    fprintf('  WARNING: GM crosses zero at zn = %.1f mm\n\n', zn_GM_zero*1e3);
else
    zn_GM_zero = NaN;
    fprintf('  GM remains positive across full sweep range.\n\n');
end

% =========================================================================
% SECTION 9 — PLOTS
% =========================================================================
figure('Name','Buoy Hydrostatic Stability', ...
       'Color',[0.97 0.97 0.97], 'Position',[80 60 960 700]);

subplot(2,1,1);
hold on; grid on; box on;
plot(zn_values*1e3, GM_T_ring*1e3, 'r-',  'LineWidth', 2.5, ...
    'DisplayName', 'GM_T — with ring (Ø280 mm)');
plot(zn_values*1e3, GM_T_hull*1e3, 'b--', 'LineWidth', 1.5, ...
    'DisplayName', 'GM_T — hull only (Ø200 mm)');
plot(zn_values*1e3, BM_T_ring*1e3, 'b-',  'LineWidth', 1.5, ...
    'DisplayName', 'BM_T — with ring');
yline(0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'GM = 0  (neutral)');
if ~isnan(zn_GM_zero)
    xline(zn_GM_zero*1e3, 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');
    text(zn_GM_zero*1e3 + 4, max(GM_T_ring)*1e3*0.12, ...
        sprintf('Capsize\n%.1f mm', zn_GM_zero*1e3), ...
        'FontSize', 9, 'Color','k', 'BackgroundColor','w', ...
        'EdgeColor',[0.7 0.7 0.7], 'Margin', 3);
end
ylabel('Height [mm]');
title('Transverse Metacentric Height GM_T vs. Additional Submergence Depth');
legend('Location','northeast','FontSize',10);
set(gca,'FontSize',11);

subplot(2,1,2);
hold on; grid on; box on;
plot(zn_values*1e3, BG_z_vals*1e3, 'g-',  'LineWidth', 2.0, ...
    'DisplayName', 'BG_z = z_g - z_b  (righting lever)');
plot(zn_values*1e3, z_b_vals*1e3,  'k-',  'LineWidth', 1.5, ...
    'DisplayName', 'z_b  (CoB below waterline)');
yline(z_g*1e3, 'c--', 'LineWidth', 1.5, ...
    'Label', sprintf('z_g = %.2f mm', z_g*1e3), ...
    'LabelHorizontalAlignment','right', ...
    'DisplayName', 'z_g  (CoG below waterline, const.)');
if ~isnan(zn_GM_zero)
    xline(zn_GM_zero*1e3, 'k:', 'LineWidth', 1.5, 'HandleVisibility','off');
end
xlabel('Additional submergence depth z_n  [mm]');
ylabel('Distance below waterline [mm]  (+ve down)');
title('Centre of Buoyancy z_b and Righting Lever BG_z vs. Submergence');
legend('Location','northeast','FontSize',10);
set(gca,'FontSize',11);

sgtitle('Buoy Hydrostatic Stability — Fossen MSS Toolbox', ...
    'FontSize',13,'FontWeight','bold');

% =========================================================================
% LOCAL HELPER
% =========================================================================
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end