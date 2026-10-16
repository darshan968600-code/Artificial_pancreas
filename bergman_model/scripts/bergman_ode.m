function dxdt = bergman_ode(t, x, u, D, params)
% BERGMAN_ODE Differential equations for the 3-state Bergman Minimal Model
%
% States:
%   x(1) = G: Blood glucose [mg/dL]
%   x(2) = X: Remote insulin action in interstitial compartment [1/min]
%   x(3) = I: Plasma insulin [uU/mL]
%
% Inputs:
%   u: Exogenous insulin infusion rate
%   D: Gut glucose appearance rate [mg/dL/min]

G = x(1);
X = x(2);
I = x(3);

Gb = params.Gb;
Ib = params.Ib;
p1 = params.p1;
p2 = params.p2;
p3 = params.p3;
n  = params.n;

% Rates of change (per minute)
dG = -p1 * (G - Gb) - X * G + D;
dX = -p2 * X + p3 * (I - Ib);
dI = -n  * (I - Ib) + u;

dxdt = [dG; dX; dI];

end
