function params = bergman_params()
% BERGMAN_PARAMS Returns parameters for the Bergman Minimal Model
% Based on standard literature values and user's original implementation.

params = struct();

% Basal Steady-State Values
params.Gb = 90.0;          % Basal blood glucose [mg/dL]
params.Ib = 7.0;           % Basal plasma insulin [uU/mL]

% Kinetic Parameters (per minute)
params.p1 = 0.0287;        % Glucose effectiveness [1/min]
params.p2 = 0.0142;        % Insulin action clearance rate [1/min]
params.p3 = 1.5e-5;        % Insulin sensitivity index [mL/(uU min^2)]
params.n  = 0.2814;        % Insulin disappearance rate [1/min]

% Initial State Vector at Basal Equilibrium
% x = [G; X; I]
params.x0 = [params.Gb; 0.0; params.Ib];

% Nominal Basal Infusion Rate to maintain Ib
params.u_basal = 0.0;      % In the differential equation dI/dt = -n*(I-Ib) + u, at equilibrium I = Ib requires u = 0 (or u represents deviation from basal)

% Simulation Settings
params.Ts_sec = 60;        % 1-minute sampling interval (in seconds)
params.Ts_min = 1;         % 1-minute sampling interval (in minutes)

end
