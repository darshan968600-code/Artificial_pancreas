function [dxdt, G, G_mgdl, UI, UG] = hovorka_ode(t, x, u_val, D_meal_val, p)
% HOVORKA_ODE Computes state derivatives for the Hovorka (2004) T1D model.
%
% States x (11x1):
%   x(1)  - S1: Subcutaneous insulin depot 1 [mU]
%   x(2)  - S2: Subcutaneous insulin depot 2 [mU]
%   x(3)  - I:  Plasma insulin concentration [mU/L]
%   x(4)  - x1: Remote insulin action on transport [min^-1]
%   x(5)  - x2: Remote insulin action on disposal [min^-1]
%   x(6)  - x3: Remote insulin action on EGP suppression [-]
%   x(7)  - Q1: Accessible compartment glucose mass [mmol]
%   x(8)  - Q2: Non-accessible compartment glucose mass [mmol]
%   x(9)  - D1: Gut carbohydrate compartment 1 [mmol]
%   x(10) - D2: Gut carbohydrate compartment 2 [mmol]
%   x(11) - G_meas: Subcutaneous CGM glucose concentration [mmol/L]
%
% Inputs:
%   u_val      - Insulin infusion rate [mU/min] (can be function handle or scalar)
%   D_meal_val - Meal ingestion rate [mmol/min] (can be function handle or scalar)
%   p          - Parameters struct from hovorka_params.m

% Handle scalar or function-handle inputs
if isa(u_val, 'function_handle')
    u = u_val(t);
else
    u = u_val;
end

if isa(D_meal_val, 'function_handle')
    D_meal = D_meal_val(t);
else
    D_meal = D_meal_val;
end

% Prevent negative inputs
u = max(0, u);
D_meal = max(0, D_meal);

% Unpack states
S1     = max(0, x(1));
S2     = max(0, x(2));
I      = max(0, x(3));
x1     = max(0, x(4));
x2     = max(0, x(5));
x3     = max(0, x(6));
Q1     = max(0.01, x(7));
Q2     = max(0.01, x(8));
D1     = max(0, x(9));
D2     = max(0, x(10));
G_meas = max(0.01, x(11));

%% 1. Subcutaneous Insulin Absorption Subsystem
UI = S2 / p.tau_S; % Insulin appearance rate [mU/min]
dS1 = u - S1 / p.tau_S;
dS2 = (S1 - S2) / p.tau_S;

%% 2. Plasma Insulin Kinetics Subsystem
dI = UI / p.V_I - p.ke * I;

%% 3. Insulin Action Subsystem
dx1 = -p.ka1 * x1 + p.kb1 * I;
dx2 = -p.ka2 * x2 + p.kb2 * I;
dx3 = -p.ka3 * x3 + p.kb3 * I;

%% 4. Gut Glucose Absorption Subsystem
UG = D2 / p.tau_D; % Glucose appearance rate [mmol/min]
dD1 = D_meal * p.A_G - D1 / p.tau_D;
dD2 = (D1 - D2) / p.tau_D;

%% 5. Glucose Kinetics Subsystem
G = Q1 / p.V_G; % Plasma glucose [mmol/L]
G_mgdl = G * p.MMOLL_TO_MGDL; % [mg/dL]

% Non-insulin dependent glucose flux (CNS):
if G >= 4.5
    F01_c = p.F01;
else
    F01_c = p.F01 * (G / 4.5);
end

% Renal glucose clearance:
if G >= 9.0
    FR = 0.003 * (G - 9.0) * p.V_G;
else
    FR = 0;
end

% Endogenous glucose production:
EGP = max(0, p.EGP0 * (1 - x3));

% Compartment mass balances:
dQ1 = -(F01_c + x1 * Q1) + p.k12 * Q2 - FR + EGP + UG;
dQ2 = x1 * Q1 - (p.k12 + x2) * Q2;

%% 6. CGM Sensor Delay
dG_meas = (G - G_meas) / p.tau_cgm;

%% Pack derivatives
dxdt = [dS1; dS2; dI; dx1; dx2; dx3; dQ1; dQ2; dD1; dD2; dG_meas];

end
