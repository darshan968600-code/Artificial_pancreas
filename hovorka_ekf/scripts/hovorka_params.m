function p = hovorka_params(BW, G_target_mgdl)
% HOVORKA_PARAMS Defines physiological parameters and analytical basal
% steady-state for the Hovorka (2004) Type 1 Diabetes model.
%
% References:
% Hovorka et al., "Nonlinear model predictive control of glucose concentration 
% in subjects with type 1 diabetes", Physiol. Meas. 25 (2004) 905-920.
%
% Inputs:
%   BW             - Patient body weight in kg (default: 70 kg)
%   G_target_mgdl  - Basal target glucose in mg/dL (default: 100 mg/dL)
%
% Output:
%   p              - Struct containing all parameters, unit conversions, 
%                    and analytical steady-state values.

if nargin < 1 || isempty(BW)
    BW = 70; % kg (standard adult)
end
if nargin < 2 || isempty(G_target_mgdl)
    G_target_mgdl = 100; % mg/dL
end

%% 1. Basic conversions
p.BW = BW;                                    % Body weight [kg]
p.MOL_WT_GLUCOSE = 180.1557;                 % g/mol
p.MGDL_TO_MMOLL  = 1 / 18.01557;             % 1 mg/dL = 0.0555 mmol/L
p.MMOLL_TO_MGDL  = 18.01557;                 % 1 mmol/L = 18.01557 mg/dL
p.U_TO_MU        = 1000;                     % 1 U = 1000 mU
p.UH_TO_MUMIN    = 1000 / 60;                % 1 U/h = 16.6667 mU/min
p.MUMIN_TO_UH    = 60 / 1000;                % 1 mU/min = 0.06 U/h

p.G_target_mgdl  = G_target_mgdl;            % mg/dL
p.G_target_mmoll = G_target_mgdl * p.MGDL_TO_MMOLL; % mmol/L

%% 2. Glucose Subsystem Parameters
p.V_G = 0.16 * p.BW;                         % Glucose distribution volume [L]
p.k12 = 0.066;                               % Transfer rate Q1 <-> Q2 [min^-1]
p.F01 = 0.0097 * p.BW;                       % Non-insulin-dependent glucose flux [mmol/min]
p.EGP0 = 0.0161 * p.BW;                      % Endogenous glucose production at 0 insulin [mmol/min]

%% 3. Insulin Action Subsystem Parameters
% Deactivation rate constants [min^-1]
p.ka1 = 0.006;
p.ka2 = 0.06;
p.ka3 = 0.03;

% Insulin sensitivities (Hovorka 2004 nominal values)
p.S_IT = 51.2e-4;                            % Transport sensitivity [min^-1 / (mU/L)]
p.S_ID = 8.2e-4;                             % Disposal sensitivity [min^-1 / (mU/L)]
p.S_IE = 520e-4;                             % EGP suppression sensitivity [L / mU]

% Activation rate constants kb = ka * S [min^-1 / (mU/L)]
p.kb1 = p.ka1 * p.S_IT;
p.kb2 = p.ka2 * p.S_ID;
p.kb3 = p.ka3 * p.S_IE;

%% 4. Insulin Kinetics & Subcutaneous Absorption Parameters
p.V_I = 0.12 * p.BW;                         % Insulin distribution volume [L]
p.ke  = 0.138;                               % Insulin elimination rate [min^-1]
p.tau_S = 55;                                % Subcutaneous insulin absorption time constant [min]

%% 5. Gut Glucose Absorption Parameters
p.tau_D = 40;                                % Meal absorption time constant [min]
p.A_G   = 0.80;                              % Carbohydrate bioavailability [-]

%% 6. Continuous Glucose Monitor (CGM) Sensor Delay
p.tau_cgm = 10;                              % Subcutaneous sensing delay [min]

%% 7. Analytical Basal Steady-State Calculation
% At basal state (no meal: U_G = 0), glucose is constant at G_target:
% G0 = G_target_mmoll (>= 4.5 mmol/L, so F01_c = F01; < 9 mmol/L, so F_R = 0)
% Q1_0 = G0 * V_G
% x1_0 = S_IT * I0
% x2_0 = S_ID * I0
% x3_0 = S_IE * I0
%
% Steady state equation:
% EGP0*(1 - S_IE*I0) - F01 - (S_IT*S_ID*Q1_0 * I0^2) / (k12 + S_ID*I0) = 0
% Multiplied by (k12 + S_ID*I0), this is quadratic in I0: a*I0^2 + b*I0 + c = 0

G0 = p.G_target_mmoll;
Q1_0 = G0 * p.V_G;
c0 = p.EGP0 - p.F01;

% Quadratic coefficients:
% a * I0^2 + b * I0 + c = 0
a_quad = -(p.EGP0 * p.S_IE * p.S_ID + p.S_IT * p.S_ID * Q1_0);
b_quad = c0 * p.S_ID - p.k12 * p.EGP0 * p.S_IE;
c_quad = c0 * p.k12;

% Solve quadratic equation:
roots_I0 = roots([a_quad, b_quad, c_quad]);
% Pick the positive real root
I0 = roots_I0(roots_I0 > 0 & imag(roots_I0) == 0);
if isempty(I0)
    error('Could not find a valid positive steady-state plasma insulin concentration.');
end
I0 = I0(1); % Plasma insulin [mU/L]

% Remote insulin actions at steady state:
x1_0 = p.S_IT * I0;
x2_0 = p.S_ID * I0;
x3_0 = p.S_IE * I0;

% Non-accessible glucose compartment Q2:
Q2_0 = (x1_0 * Q1_0) / (p.k12 + x2_0);

% Insulin appearance rate at steady state [mU/min]:
UI_0 = p.V_I * p.ke * I0;

% Basal insulin infusion rate:
u_basal_mumin = UI_0;                        % [mU/min]
u_basal_Uh    = u_basal_mumin * p.MUMIN_TO_UH; % [U/h]

% Subcutaneous insulin depots:
S1_0 = u_basal_mumin * p.tau_S;              % [mU]
S2_0 = u_basal_mumin * p.tau_S;              % [mU]

% Meal absorption compartments at steady state (no meal):
D1_0 = 0;
D2_0 = 0;

% Sensor glucose:
G_meas_0 = G0;

%% 8. Pack steady-state values into p
p.ss.I0            = I0;
p.ss.x1_0          = x1_0;
p.ss.x2_0          = x2_0;
p.ss.x3_0          = x3_0;
p.ss.Q1_0          = Q1_0;
p.ss.Q2_0          = Q2_0;
p.ss.G0            = G0;
p.ss.G0_mgdl       = G0 * p.MMOLL_TO_MGDL;
p.ss.S1_0          = S1_0;
p.ss.S2_0          = S2_0;
p.ss.D1_0          = D1_0;
p.ss.D2_0          = D2_0;
p.ss.G_meas_0      = G_meas_0;
p.ss.u_basal_mumin = u_basal_mumin;
p.ss.u_basal_Uh    = u_basal_Uh;

% Full steady-state state vector:
% x = [S1, S2, I, x1, x2, x3, Q1, Q2, D1, D2, G_meas]'
p.ss.x0 = [S1_0; S2_0; I0; x1_0; x2_0; x3_0; Q1_0; Q2_0; D1_0; D2_0; G_meas_0];

% State names for documentation
p.state_names = { ...
    'S1 (Subcut insulin depot 1 [mU])', ...
    'S2 (Subcut insulin depot 2 [mU])', ...
    'I (Plasma insulin [mU/L])', ...
    'x1 (Transport insulin action [min^-1])', ...
    'x2 (Disposal insulin action [min^-1])', ...
    'x3 (EGP suppression action [-])', ...
    'Q1 (Accessible glucose [mmol])', ...
    'Q2 (Non-accessible glucose [mmol])', ...
    'D1 (Gut meal compartment 1 [mmol])', ...
    'D2 (Gut meal compartment 2 [mmol])', ...
    'G_meas (CGM measured glucose [mmol/L])' ...
};

end
