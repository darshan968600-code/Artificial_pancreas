function [mpc_obj, plant_ss_d, p] = design_mpc_controller()
% DESIGN_MPC_CONTROLLER Linearizes the Hovorka plant at basal steady-state
% and designs the Model Predictive Controller (MPC) with full state feedback
% (Insulin-on-Board IOB tracking) to achieve:
%   - Peak Glucose strictly < 180 mg/dL with Meal Preview
%   - Minimum Glucose strictly > 70 mg/dL (Zero Hypoglycemia)
%
% States x (11x1):
%   [S1; S2; I; x1; x2; x3; Q1; Q2; D1; D2; G_meas]

%% 1. Load Parameters and Basal Steady-State
p = hovorka_params(70, 100);

x0 = p.ss.x0;
u0_Uh = p.ss.u_basal_Uh;      % ~0.392 U/h
d0_gmin = 0;                  % g/min

%% 2. Compute Numerical Jacobian (Linearization)
nx = 11;

u_Uh_to_mumin = p.UH_TO_MUMIN;
carbs_gmin_to_mmolmin = 1000 / p.MOL_WT_GLUCOSE;

f_eval = @(x, u_Uh, d_gmin) hovorka_ode(0, x, ...
    u_Uh * u_Uh_to_mumin, ...
    d_gmin * carbs_gmin_to_mmolmin, p);

h_x = 1e-6 * max(1, abs(x0));
h_u = 1e-6;
h_d = 1e-6;

% A matrix (df/dx)
A = zeros(nx, nx);
for j = 1:nx
    x_pert_plus  = x0; x_pert_plus(j)  = x_pert_plus(j) + h_x(j);
    x_pert_minus = x0; x_pert_minus(j) = x_pert_minus(j) - h_x(j);
    f_plus  = f_eval(x_pert_plus,  u0_Uh, d0_gmin);
    f_minus = f_eval(x_pert_minus, u0_Uh, d0_gmin);
    A(:, j) = (f_plus - f_minus) / (2 * h_x(j));
end

% B_u matrix (df/du in U/h)
f_u_plus  = f_eval(x0, u0_Uh + h_u, d0_gmin);
f_u_minus = f_eval(x0, u0_Uh - h_u, d0_gmin);
B_u = (f_u_plus - f_u_minus) / (2 * h_u);

% B_d matrix (df/dd in g/min)
f_d_plus  = f_eval(x0, u0_Uh, d0_gmin + h_d);
f_d_minus = f_eval(x0, u0_Uh, d0_gmin - h_d);
B_d = (f_d_plus - f_d_minus) / (2 * h_d);

% Output matrix C: y = G_meas * 18.01557 [mg/dL]
C = zeros(1, nx);
C(1, 11) = p.MMOLL_TO_MGDL;
D = [0, 0];

%% 3. Continuous and Discrete State-Space Models
plant_ss_c = ss(A, [B_u, B_d], C, D);
plant_ss_c.InputName  = {'Insulin_Infusion', 'Meal_Carb_Rate'};
plant_ss_c.InputUnit  = {'U/h', 'g/min'};
plant_ss_c.OutputName = {'Glucose'};
plant_ss_c.OutputUnit = {'mg/dL'};
plant_ss_c.InputGroup.ManipulatedVariables = 1;
plant_ss_c.InputGroup.MeasuredDisturbances = 2;

Ts = 5; % minutes
plant_ss_d = c2d(plant_ss_c, Ts, 'zoh');

%% 4. Configure MATLAB MPC Object
Np = 48; % 4-hour prediction horizon captures full insulin absorption & clearance
Nc = 4;  % 20-min control horizon

mpc_obj = mpc(plant_ss_d, Ts, Np, Nc);

% Set nominal operating values
mpc_obj.Model.Nominal.X = x0;
mpc_obj.Model.Nominal.U = [u0_Uh; 0];
mpc_obj.Model.Nominal.Y = 100;

% Set custom state estimator (enables full state feedback & IOB tracking)
setEstimator(mpc_obj, 'custom');

% Actuator Constraints (Allows fast prandial bolus while respecting pump maximum)
mpc_obj.ManipulatedVariables(1).Min = 0.0;     % Minimum infusion: 0 U/h
mpc_obj.ManipulatedVariables(1).Max = 4.0;     % Prandial delivery ceiling: 4.0 U/h
mpc_obj.ManipulatedVariables(1).RateMin = -2.0; % Slew rate limit
mpc_obj.ManipulatedVariables(1).RateMax = 2.0;

% Tuning Weights: State-Feedback Anti-Hypo Tuning
mpc_obj.Weights.OutputVariables = 1.0;          % Glucose tracking penalty
mpc_obj.Weights.ManipulatedVariables = 0.5;     % Balanced insulin penalty
mpc_obj.Weights.ManipulatedVariablesRate = 5.0; % Delivery smoothness penalty

fprintf('Safe State-Feedback MPC Controller configured!\n');
fprintf('  Prediction Horizon Np = %d steps (%d min / %.1f hours)\n', Np, Np*Ts, (Np*Ts)/60);
fprintf('  Control Horizon Nc    = %d steps (%d min)\n', Nc, Nc*Ts);
fprintf('  Sample Time Ts        = %d min\n', Ts);
fprintf('  Actuator Bounds:        [%.1f, %.1f] U/h\n', ...
    mpc_obj.ManipulatedVariables(1).Min, mpc_obj.ManipulatedVariables(1).Max);

end
