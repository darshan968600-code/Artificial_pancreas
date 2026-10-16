% VERIFY_STEADY_STATE Verifies basal equilibrium of Hovorka plant
clc; clear; close all;

fprintf('=========================================================\n');
fprintf('  HOVORKA MODEL: BASAL STEADY-STATE VERIFICATION\n');
fprintf('=========================================================\n');

% Robust project root and path resolution
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts')
    project_root = fileparts(this_dir);
else
    project_root = this_dir;
end
addpath(fullfile(project_root, 'scripts'));

% Load parameters for 70 kg adult, target 100 mg/dL
p = hovorka_params(70, 100);

fprintf('Patient Weight:        %.1f kg\n', p.BW);
fprintf('Target Glucose:        %.2f mg/dL (%.2f mmol/L)\n', p.ss.G0_mgdl, p.ss.G0);
fprintf('Basal Insulin Rate:    %.3f mU/min (%.3f U/h)\n', p.ss.u_basal_mumin, p.ss.u_basal_Uh);
fprintf('Basal Plasma Insulin:  %.2f mU/L\n', p.ss.I0);
fprintf('Accessible Glucose Q1: %.2f mmol\n', p.ss.Q1_0);
fprintf('Non-access Glucose Q2: %.2f mmol\n', p.ss.Q2_0);
fprintf('Insulin Action x1:     %.4e min^-1\n', p.ss.x1_0);
fprintf('Insulin Action x2:     %.4e min^-1\n', p.ss.x2_0);
fprintf('Insulin Action x3:     %.4f\n', p.ss.x3_0);

% Check derivative at steady state:
dxdt_0 = hovorka_ode(0, p.ss.x0, p.ss.u_basal_mumin, 0, p);
max_deriv_err = max(abs(dxdt_0));

fprintf('---------------------------------------------------------\n');
fprintf('Maximum |dx/dt| at t=0: %.3e\n', max_deriv_err);
if max_deriv_err < 1e-10
    fprintf('[PASS] Analytical steady-state matches ODE derivatives to machine precision!\n');
else
    fprintf('[FAIL] State derivative is non-zero (%.3e)!\n', max_deriv_err);
    return;
end

% Run 24-hour simulation (1440 minutes) with basal infusion and no meal
tspan = [0, 1440]; % 24 hours in minutes
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, x] = ode15s(@(t, x) hovorka_ode(t, x, p.ss.u_basal_mumin, 0, p), tspan, p.ss.x0, opts);

G_mgdl_sim = (x(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
max_drift = max(abs(G_mgdl_sim - 100));

fprintf('24-Hour Simulation Maximum Glucose Drift: %.4e mg/dL\n', max_drift);
if max_drift < 1e-5
    fprintf('[PASS] Glucose remained rock-solid flat at 100 mg/dL for 24 hours!\n');
else
    fprintf('[FAIL] Glucose drifted during basal run by %.4e mg/dL!\n', max_drift);
end
fprintf('=========================================================\n');
