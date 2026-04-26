% =========================================================================
% BUOY DYNAMIC MATRICES (Rigid Body + Hydrostatics)
% =========================================================================

% 1. PHYSICAL PROPERTIES (From previous SolidWorks/Hydro outputs)
m         = 1.4431;      % Total mass (including ballast) [kg]
KG        = 0.08971;     % CoG above keel [m] (Preserved precision)
k44       = 0.0692;      % Roll radius of gyration @ CG [m]
k55       = 0.0692;      % Pitch radius of gyration @ CG [m]
k66       = 0.0823;      % Yaw radius of gyration @ CG [m]
nabla_val = 0.001465;    % Displaced volume [m^3]
A_wp_val  = 0.06158;     % Waterplane area [m^2]
GMT_val   = 0.17283;     % Transverse GM [m]
GML_val   = 0.17283;     % Longitudinal GM (axisymmetric) [m]

% 2. COORDINATE VECTORS
% In Fossen's convention, z is positive DOWN. 
% Origin (CO) is at the keel. Since CG is physically above the keel, 
% its z-coordinate is negative.
r_G  = [0; 0; -KG];     
r_bP = [0; 0; 0];       % Point P is the Keel (CO)
nu2  = [0; 0; 0];       % Angular velocity [p, q, r] at rest
x_F  = 0;               % Center of flotation (center of buoy)

% 3. COMPUTE MATRICES
% Compute Rigid-Body Mass Matrix (MRB) and Coriolis Matrix (CRB)
[MRB, CRB] = rbody(m, k44, k55, k66, nu2, r_G);

% Compute Hydrostatic Stiffness Matrix (G)
G = Gmtrx(nabla_val, A_wp_val, GMT_val, GML_val, x_F, r_bP);

% 4. DISPLAY OUTPUTS
clc;
fprintf('<strong>=== BUOY SYSTEM MATRICES (Reference: Keel, Z-Down) ===</strong>\n');

fprintf('\n<strong>1. Hydrostatic Stiffness Matrix (G) [N/m, Nm/rad]</strong>\n');
% Describes the "buoyancy springs" in Heave, Roll, and Pitch
disp(array2table(round(G, 4), 'VariableNames', {'Surge','Sway','Heave','Roll','Pitch','Yaw'}));

fprintf('\n<strong>2. Rigid-Body Mass Matrix (MRB) [kg, kg·m^2]</strong>\n');
% Includes mass and inertia shifted from CG to Keel
disp(array2table(round(MRB, 5), 'VariableNames', {'Surge','Sway','Heave','Roll','Pitch','Yaw'}));

fprintf('\n<strong>3. Coriolis Matrix (CRB) at nu=[0,0,0]</strong>\n');
disp(CRB);