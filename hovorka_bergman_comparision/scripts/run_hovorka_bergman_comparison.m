function [metrics_table, results] = run_hovorka_bergman_comparison(meal_carbs_g, meal_dur, meal_start_time, t_end)
% RUN_HOVORKA_BERGMAN_COMPARISON Master comparative benchmark evaluating
% the Bergman Minimal Model (3-state) vs the Hovorka Model (11-state)
% under an identical, physiologically calibrated meal carbohydrate disturbance
% with a commercial Hybrid Closed-Loop (HCL) prandial bolus (150 U/h speed).
%
% Compares 6 distinct closed-loop cases:
%   1. Bergman Open-Loop (Calibrated Bolus + Fixed Basal)
%   2. Bergman MPC Without Meal Preview (Calibrated Bolus + Reactive Basal)
%   3. Bergman MPC With Meal Preview (Calibrated Bolus + Proactive Basal)
%   4. Hovorka Open-Loop (Prandial Bolus + Fixed Basal)
%   5. Hovorka MPC Without Meal Preview (Prandial Bolus + Reactive Basal)
%   6. Hovorka MPC With Meal Preview (Prandial Bolus + Proactive Basal)
%
% Usage:
%   run_hovorka_bergman_comparison()             -> Interactive prompt (defaults: 75g, 15 min)
%   run_hovorka_bergman_comparison(75)           -> 75g CHO over 15 min
%   run_hovorka_bergman_comparison(100, 20)      -> 100g CHO over 20 min
%   run_hovorka_bergman_comparison(50, 15, 60, 360) -> 50g CHO, 6h timeline

clc;
fprintf('========================================================================\n');
fprintf('  ARTIFICIAL PANCREAS BENCHMARK: BERGMAN (3-STATE) vs HOVORKA (11-STATE)\n');
fprintf('  Physiologically Calibrated Disturbance + 150 U/h Prandial Bolus\n');
fprintf('========================================================================\n');

%% 1. Interactive / Default Inputs
if nargin < 1 || isempty(meal_carbs_g)
    try
        user_val = input('Enter meal carbohydrate amount in grams [Default: 75]: ', 's');
        if isempty(strtrim(user_val))
            meal_carbs_g = 75;
        else
            meal_carbs_g = str2double(user_val);
            if isnan(meal_carbs_g) || meal_carbs_g <= 0
                fprintf('  Invalid input. Defaulting to 75g CHO.\n');
                meal_carbs_g = 75;
            end
        end
    catch
        meal_carbs_g = 75;
    end
end

if nargin < 2 || isempty(meal_dur)
    try
        dur_val = input('Enter meal ingestion duration in minutes [Default: 15]: ', 's');
        if isempty(strtrim(dur_val))
            meal_dur = 15;
        else
            meal_dur = str2double(dur_val);
            if isnan(meal_dur) || meal_dur <= 0
                fprintf('  Invalid input. Defaulting to 15 min duration.\n');
                meal_dur = 15;
            end
        end
    catch
        meal_dur = 15;
    end
end

if nargin < 3 || isempty(meal_start_time)
    meal_start_time = 60; % min (1-hour basal equilibrium baseline)
end

if nargin < 4 || isempty(t_end)
    t_end = 360; % min (6 hours total simulation time)
end

%% 2. Shared Prandial Bolus Parameters (Commercial HCL Standard)
% Patient ICR: 24.0 g/U for standard 70kg adult
ICR = 24.0; % g/U
bolus_fraction = 0.85; % 85% upfront bolus
bolus_U = (meal_carbs_g / ICR) * bolus_fraction;

% Bolus motor delivery speed: 150 U/h
bolus_speed_Uh = 150.0; % U/h
% Bolus delivery duration in minutes:
t_bolus_min = (bolus_U / bolus_speed_Uh) * 60; % duration in minutes

fprintf('\nSimulation Benchmark Scenario:\n');
fprintf('  - Meal Carbohydrates:   %.1f grams\n', meal_carbs_g);
fprintf('  - Meal Ingestion Time:  %.1f min (at t = %.0f min)\n', meal_dur, meal_start_time);
fprintf('  - Prandial Bolus Dose:  %.2f Units (85%% of %.1fg @ ICR=%.0f)\n', bolus_U, meal_carbs_g, ICR);
fprintf('  - Bolus Delivery Speed: %.1f U/h (Duration: %.1f seconds)\n', bolus_speed_Uh, t_bolus_min * 60);
fprintf('  - Total Duration:       %.1f min (%.1f hours)\n', t_end, t_end/60);

%% 3. Setup Paths & Directories
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts')
    project_root = fileparts(this_dir);
else
    project_root = this_dir;
end
workspace_root = fileparts(project_root);

bergman_scripts = fullfile(workspace_root, 'bergman_model', 'scripts');
bergman_models  = fullfile(workspace_root, 'bergman_model', 'models');
hovorka_scripts = fullfile(workspace_root, 'hovorka_phase1', 'scripts');

addpath(fullfile(project_root, 'scripts'));
addpath(bergman_scripts);
addpath(bergman_models);
addpath(hovorka_scripts);

fig_dir  = fullfile(project_root, 'results', 'figures');
data_dir = fullfile(project_root, 'results', 'data');
if ~exist(fig_dir, 'dir'),  mkdir(fig_dir);  end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

save_clean = @(fh, fp) clean_save_helper(fh, fp);
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

%% =====================================================================
%% PART A: BERGMAN MINIMAL MODEL (3-STATE) SIMULATION (CALIBRATED)
%% =====================================================================
fprintf('\n------------------------------------------------------------------------\n');
fprintf('  [PART 1/2] SIMULATING BERGMAN MINIMAL MODEL (PHYSIOLOGICALLY CALIBRATED)\n');
fprintf('------------------------------------------------------------------------\n');

p_berg = bergman_params();
Gb_berg = p_berg.Gb; % 90 mg/dL
Ib_berg = p_berg.Ib; % 7 uU/mL
x0_berg = [Gb_berg; 0.0; Ib_berg];

Ts_berg = 1; % min
N_berg = t_end / Ts_berg;
time_berg = (0:N_berg)' * Ts_berg;

% Calibrated Meal Disturbance for Bergman:
% Total mass = meal_carbs_g * 0.80 bioavailable into V_G = 120 dL
% Integral(D) = D_peak * e * T_peak = delta_G_total
D_peak_berg = (meal_carbs_g * 0.80 * 1000) / (120 * exp(1) * 30);
T_peak_berg = 30; % min
D_berg = zeros(N_berg + 1, 1);
for k = 1:(N_berg + 1)
    t_curr = time_berg(k);
    if t_curr >= meal_start_time
        tau = t_curr - meal_start_time;
        D_berg(k) = D_peak_berg * (tau / T_peak_berg) * exp(1 - tau / T_peak_berg);
    end
end

% Unit conversion for Bergman infusion:
% 1 U in V_I (8,400 mL) = 119.05 uU/mL
% 1 U/h = 1.984 uU/mL/min
u_conv_berg = 1.984;
bolus_speed_berg = bolus_speed_Uh * u_conv_berg; % ~297.6 uU/mL/min

% Load Bergman MPC Controllers
load(fullfile(bergman_models, 'mpc_baseline.mat'), 'mpcobj');
load(fullfile(bergman_models, 'mpc_meal_preview.mat'), 'mpcobj_md');

% --- 1. Bergman Open-Loop (Bolus + Basal) ---
fprintf('  Simulating Bergman Case 1: Open-Loop Baseline (Bolus + Basal)...\n');
x_berg_OL = zeros(N_berg + 1, 3);
u_berg_OL = zeros(N_berg, 1);
x_berg_OL(1, :) = x0_berg';

for k = 1:N_berg
    t_k = time_berg(k);
    t_next = time_berg(k+1);
    D_k = D_berg(k);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode45(@(t, x) bergman_ode(t, x, bolus_speed_berg, D_k, p_berg), [t_k, t_split], x_berg_OL(k, :)', ode_opts);
        [~, s2] = ode45(@(t, x) bergman_ode(t, x, 0.0, D_k, p_berg), [t_split, t_next], s1(end, :)', ode_opts);
        x_berg_OL(k+1, :) = s2(end, :);
        u_berg_OL(k) = (bolus_U / (Ts_berg / 60));
    else
        [~, s] = ode45(@(t, x) bergman_ode(t, x, 0.0, D_k, p_berg), [t_k, t_next], x_berg_OL(k, :)', ode_opts);
        x_berg_OL(k+1, :) = s(end, :);
        u_berg_OL(k) = 0.0;
    end
end

% --- 2. Bergman MPC Without Preview (Bolus + Reactive Basal) ---
fprintf('  Simulating Bergman Case 2: MPC Without Meal Preview (Reactive)...\n');
x_berg_NP = zeros(N_berg + 1, 3);
u_berg_NP = zeros(N_berg, 1);
x_berg_NP(1, :) = x0_berg';
xmpc_berg_NP = mpcstate(mpcobj);

for k = 1:N_berg
    t_k = time_berg(k);
    t_next = time_berg(k+1);
    y_meas = x_berg_NP(k, 1);
    u_opt = mpcmove(mpcobj, xmpc_berg_NP, y_meas, Gb_berg);
    D_k = D_berg(k);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode45(@(t, x) bergman_ode(t, x, u_opt + bolus_speed_berg, D_k, p_berg), [t_k, t_split], x_berg_NP(k, :)', ode_opts);
        [~, s2] = ode45(@(t, x) bergman_ode(t, x, u_opt, D_k, p_berg), [t_split, t_next], s1(end, :)', ode_opts);
        x_berg_NP(k+1, :) = s2(end, :);
        u_berg_NP(k) = (u_opt / u_conv_berg) + (bolus_U / (Ts_berg / 60));
    else
        [~, s] = ode45(@(t, x) bergman_ode(t, x, u_opt, D_k, p_berg), [t_k, t_next], x_berg_NP(k, :)', ode_opts);
        x_berg_NP(k+1, :) = s(end, :);
        u_berg_NP(k) = u_opt / u_conv_berg;
    end
end

% --- 3. Bergman MPC With Meal Preview (Bolus + Proactive Basal) ---
fprintf('  Simulating Bergman Case 3: MPC With Meal Preview (Proactive)...\n');
x_berg_MP = zeros(N_berg + 1, 3);
u_berg_MP = zeros(N_berg, 1);
x_berg_MP(1, :) = x0_berg';
xmpc_berg_MP = mpcstate(mpcobj_md);
Np_berg = mpcobj_md.PredictionHorizon;

for k = 1:N_berg
    t_k = time_berg(k);
    t_next = time_berg(k+1);
    y_meas = x_berg_MP(k, 1);
    v_preview = zeros(Np_berg + 1, 1);
    for h = 0:Np_berg
        v_preview(h + 1) = D_berg(min(k + h, N_berg + 1));
    end
    u_opt = mpcmove(mpcobj_md, xmpc_berg_MP, y_meas, Gb_berg, v_preview);
    D_k = D_berg(k);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode45(@(t, x) bergman_ode(t, x, u_opt + bolus_speed_berg, D_k, p_berg), [t_k, t_split], x_berg_MP(k, :)', ode_opts);
        [~, s2] = ode45(@(t, x) bergman_ode(t, x, u_opt, D_k, p_berg), [t_split, t_next], s1(end, :)', ode_opts);
        x_berg_MP(k+1, :) = s2(end, :);
        u_berg_MP(k) = (u_opt / u_conv_berg) + (bolus_U / (Ts_berg / 60));
    else
        [~, s] = ode45(@(t, x) bergman_ode(t, x, u_opt, D_k, p_berg), [t_k, t_next], x_berg_MP(k, :)', ode_opts);
        x_berg_MP(k+1, :) = s(end, :);
        u_berg_MP(k) = u_opt / u_conv_berg;
    end
end

G_berg_OL = x_berg_OL(:, 1);
G_berg_NP = x_berg_NP(:, 1);
G_berg_MP = x_berg_MP(:, 1);

%% =====================================================================
%% PART B: HOVORKA NONLINEAR MODEL (11-STATE) SIMULATION
%% =====================================================================
fprintf('\n------------------------------------------------------------------------\n');
fprintf('  [PART 2/2] SIMULATING HOVORKA MODEL (PHYSIOLOGICAL PLANT + BOLUS)\n');
fprintf('------------------------------------------------------------------------\n');

p_hov = hovorka_params(70, 100);
[mpc_hov, plant_hov_d] = design_mpc_controller();

Ts_hov = 5; % min
N_hov = t_end / Ts_hov;
time_hov = (0:N_hov)' * Ts_hov;

meal_rate_gmin = meal_carbs_g / meal_dur;
d_hov_gmin = zeros(N_hov + 1, 1);
for k = 1:(N_hov + 1)
    t_curr = time_hov(k);
    if t_curr >= meal_start_time && t_curr < (meal_start_time + meal_dur)
        d_hov_gmin(k) = meal_rate_gmin;
    end
end

Np_hov = mpc_hov.PredictionHorizon;
target_hov_mgdl = 100;

% --- 4. Hovorka Open-Loop (Bolus + Basal) ---
fprintf('  Simulating Hovorka Case 1: Open-Loop Baseline (Bolus + Basal)...\n');
x_hov_OL = zeros(N_hov + 1, 11);
u_hov_OL = zeros(N_hov, 1);
x_hov_OL(1, :) = p_hov.ss.x0';

for k = 1:N_hov
    t_k = time_hov(k);
    t_next = time_hov(k+1);
    d_curr = d_hov_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p_hov.MOL_WT_GLUCOSE);
    u_basal = p_hov.ss.u_basal_Uh;
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_basal + bolus_speed_Uh)*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_k, t_split], x_hov_OL(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_basal*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_hov_OL(k+1, :) = s2(end, :);
        u_hov_OL(k) = u_basal + (bolus_U / (Ts_hov / 60));
    else
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_basal*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_k, t_next], x_hov_OL(k, :)', ode_opts);
        x_hov_OL(k+1, :) = s(end, :);
        u_hov_OL(k) = u_basal;
    end
end

% --- 5. Hovorka MPC Without Preview (Bolus + Reactive Basal) ---
fprintf('  Simulating Hovorka Case 2: MPC Without Meal Preview...\n');
x_hov_NP = zeros(N_hov + 1, 11);
u_hov_NP = zeros(N_hov, 1);
x_hov_NP(1, :) = p_hov.ss.x0';
xmpc_hov_NP = mpcstate(mpc_hov);

for k = 1:N_hov
    t_k = time_hov(k);
    t_next = time_hov(k+1);
    y_meas = x_hov_NP(k, 11) * p_hov.MMOLL_TO_MGDL;
    d_curr = d_hov_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p_hov.MOL_WT_GLUCOSE);
    
    xmpc_hov_NP.Plant = x_hov_NP(k, :)';
    u_opt = mpcmove(mpc_hov, xmpc_hov_NP, y_meas, target_hov_mgdl, d_curr);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_k, t_split], x_hov_NP(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_hov_NP(k+1, :) = s2(end, :);
        u_hov_NP(k) = u_opt + (bolus_U / (Ts_hov / 60));
    else
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_k, t_next], x_hov_NP(k, :)', ode_opts);
        x_hov_NP(k+1, :) = s(end, :);
        u_hov_NP(k) = u_opt;
    end
end

% --- 6. Hovorka MPC With Meal Preview (Bolus + Proactive Basal) ---
fprintf('  Simulating Hovorka Case 3: MPC With Meal Preview...\n');
x_hov_MP = zeros(N_hov + 1, 11);
u_hov_MP = zeros(N_hov, 1);
x_hov_MP(1, :) = p_hov.ss.x0';
xmpc_hov_MP = mpcstate(mpc_hov);

for k = 1:N_hov
    t_k = time_hov(k);
    t_next = time_hov(k+1);
    y_meas = x_hov_MP(k, 11) * p_hov.MMOLL_TO_MGDL;
    v_preview = zeros(Np_hov + 1, 1);
    for h = 0:Np_hov
        v_preview(h + 1) = d_hov_gmin(min(k + h, N_hov + 1));
    end
    xmpc_hov_MP.Plant = x_hov_MP(k, :)';
    u_opt = mpcmove(mpc_hov, xmpc_hov_MP, y_meas, target_hov_mgdl, v_preview);
    d_curr = d_hov_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p_hov.MOL_WT_GLUCOSE);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_k, t_split], x_hov_MP(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_hov_MP(k+1, :) = s2(end, :);
        u_hov_MP(k) = u_opt + (bolus_U / (Ts_hov / 60));
    else
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p_hov.UH_TO_MUMIN, d_k_mmolmin, p_hov), ...
            [t_k, t_next], x_hov_MP(k, :)', ode_opts);
        x_hov_MP(k+1, :) = s(end, :);
        u_hov_MP(k) = u_opt;
    end
end

G_hov_OL = (x_hov_OL(:, 7) ./ p_hov.V_G) .* p_hov.MMOLL_TO_MGDL;
G_hov_NP = (x_hov_NP(:, 7) ./ p_hov.V_G) .* p_hov.MMOLL_TO_MGDL;
G_hov_MP = (x_hov_MP(:, 7) ./ p_hov.V_G) .* p_hov.MMOLL_TO_MGDL;

%% =====================================================================
%% 4. COMPUTE CLINICAL METRICS FOR ALL 6 CASES
%% =====================================================================
calc_met = @(G, t, u, Ts_m, G0) struct( ...
    'Baseline', G0, ...
    'Peak', max(G), ...
    't_Peak', t(find(G == max(G), 1)), ...
    'Excursion', max(G) - G0, ...
    'Min', min(G), ...
    'TIR', 100 * mean(G >= 70 & G <= 180), ...
    'TAR', 100 * mean(G > 180), ...
    'TBR', 100 * mean(G < 70), ...
    'TotalU', sum(u) * (Ts_m / 60) ...
);

m_b_ol = calc_met(G_berg_OL, time_berg, u_berg_OL, Ts_berg, Gb_berg);
m_b_np = calc_met(G_berg_NP, time_berg, u_berg_NP, Ts_berg, Gb_berg);
m_b_mp = calc_met(G_berg_MP, time_berg, u_berg_MP, Ts_berg, Gb_berg);

m_h_ol = calc_met(G_hov_OL, time_hov, u_hov_OL, Ts_hov, 100);
m_h_np = calc_met(G_hov_NP, time_hov, u_hov_NP, Ts_hov, 100);
m_h_mp = calc_met(G_hov_MP, time_hov, u_hov_MP, Ts_hov, 100);

red_b_mgdl = m_b_ol.Peak - m_b_mp.Peak;
red_b_pct  = (red_b_mgdl / m_b_ol.Peak) * 100;
red_h_mgdl = m_h_ol.Peak - m_h_mp.Peak;
red_h_pct  = (red_h_mgdl / m_h_ol.Peak) * 100;

%% 5. DISPLAY COMPREHENSIVE SIDE-BY-SIDE BENCHMARK TABLE
fprintf('\n===============================================================================================\n');
fprintf('  SIDE-BY-SIDE CLINICAL BENCHMARK: %.1fg CHO MEAL (150 U/h Bolus Speed)\n', meal_carbs_g);
fprintf('===============================================================================================\n');
fprintf('%-26s | %-32s | %-32s\n', 'Metric', 'BERGMAN MINIMAL MODEL (3-State)', 'HOVORKA MODEL (11-State)');
fprintf('%-26s | %-9s %-10s %-10s | %-9s %-10s %-10s\n', '', 'Open-Loop', 'No-Preview', 'Preview', 'Open-Loop', 'No-Preview', 'Preview');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Baseline Glucose (mg/dL)', m_b_ol.Baseline, m_b_np.Baseline, m_b_mp.Baseline, m_h_ol.Baseline, m_h_np.Baseline, m_h_mp.Baseline);
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Peak Glucose (mg/dL)', m_b_ol.Peak, m_b_np.Peak, m_b_mp.Peak, m_h_ol.Peak, m_h_np.Peak, m_h_mp.Peak);
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Time to Peak (min)', m_b_ol.t_Peak, m_b_np.t_Peak, m_b_mp.t_Peak, m_h_ol.t_Peak, m_h_np.t_Peak, m_h_mp.t_Peak);
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Max Excursion Delta (mg/dL)', m_b_ol.Excursion, m_b_np.Excursion, m_b_mp.Excursion, m_h_ol.Excursion, m_h_np.Excursion, m_h_mp.Excursion);
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Minimum Glucose (mg/dL)', m_b_ol.Min, m_b_np.Min, m_b_mp.Min, m_h_ol.Min, m_h_np.Min, m_h_mp.Min);
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Time in Range 70-180 (%)', m_b_ol.TIR, m_b_np.TIR, m_b_mp.TIR, m_h_ol.TIR, m_h_np.TIR, m_h_mp.TIR);
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Time Above 180 (%)', m_b_ol.TAR, m_b_np.TAR, m_b_mp.TAR, m_h_ol.TAR, m_h_np.TAR, m_h_mp.TAR);
fprintf('%-26s | %-9.1f %-10.1f %-10.1f | %-9.1f %-10.1f %-10.1f\n', 'Time Below 70 (%)', m_b_ol.TBR, m_b_np.TBR, m_b_mp.TBR, m_h_ol.TBR, m_h_np.TBR, m_h_mp.TBR);
fprintf('%-26s | %-9.2f %-10.2f %-10.2f | %-9.2f %-10.2f %-10.2f\n', 'Prandial Bolus (U)', bolus_U, bolus_U, bolus_U, bolus_U, bolus_U, bolus_U);
fprintf('%-26s | %-9.2f %-10.2f %-10.2f | %-9.2f %-10.2f %-10.2f\n', 'Total Insulin (U)', m_b_ol.TotalU, m_b_np.TotalU, m_b_mp.TotalU, m_h_ol.TotalU, m_h_np.TotalU, m_h_mp.TotalU);
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('Preview Peak Reduction:       | %.1f mg/dL (%.1f%%) vs OL             | %.1f mg/dL (%.1f%%) vs OL\n', ...
    red_b_mgdl, red_b_pct, red_h_mgdl, red_h_pct);
fprintf('===============================================================================================\n');

%% =====================================================================
%% 6. GENERATE HIGH-RESOLUTION COMPARISON FIGURES
%% =====================================================================
max_y_plot = max(260, ceil(max([m_b_ol.Peak, m_h_ol.Peak, 220])/20)*20 + 20);

% --- FIGURE 1: Glucose Comparison (Bergman Top, Hovorka Bottom) ---
fig1 = figure('Color', 'w', 'Position', [80, 60, 1100, 850], 'Visible', 'on');

subplot(2, 1, 1);
plot(time_berg, G_berg_OL, 'r--', 'LineWidth', 2.0, 'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL, TIR: %.1f%%)', m_b_ol.Peak, m_b_ol.TIR));
hold on;
plot(time_berg, G_berg_NP, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No-Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_b_np.Peak, m_b_np.TIR));
plot(time_berg, G_berg_MP, 'g-',  'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal-Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_b_mp.Peak, m_b_mp.TIR));
yline(90,  'k:', 'LineWidth', 1.0, 'DisplayName', 'Target (90 mg/dL)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180 mg/dL)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70 mg/dL)');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Target Range [70-180]');
ylabel('Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Bergman Minimal Model (Calibrated 3-State Plant) - %.1fg CHO Meal (%.2f U Bolus)', ...
    meal_carbs_g, bolus_U), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, t_end]); ylim([50, max_y_plot]);

subplot(2, 1, 2);
plot(time_hov, G_hov_OL, 'r--', 'LineWidth', 2.0, 'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL, TIR: %.1f%%)', m_h_ol.Peak, m_h_ol.TIR));
hold on;
plot(time_hov, G_hov_NP, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No-Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_h_np.Peak, m_h_np.TIR));
plot(time_hov, G_hov_MP, 'g-',  'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal-Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_h_mp.Peak, m_h_mp.TIR));
yline(100, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Target (100 mg/dL)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180 mg/dL)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70 mg/dL)');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Target Range [70-180]');
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Hovorka Physiological Model (11-State Plant) - %.1fg CHO Meal (%.2f U Bolus)', ...
    meal_carbs_g, bolus_U), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, t_end]); ylim([50, max_y_plot]);

fig1_path = fullfile(fig_dir, 'hovorka_vs_bergman_glucose_comparison.png');
save_clean(fig1, fig1_path);

% --- FIGURE 2: Direct Glucose Excursion Delta Overlay ---
fig2 = figure('Color', 'w', 'Position', [100, 80, 1050, 650], 'Visible', 'off');
plot(time_berg, G_berg_OL - Gb_berg, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Bergman Open-Loop (\\Delta: +%.1f mg/dL)', m_b_ol.Excursion));
hold on;
plot(time_berg, G_berg_MP - Gb_berg, 'm-',  'LineWidth', 2.0, 'DisplayName', sprintf('Bergman MPC Preview (\\Delta: +%.1f mg/dL)', m_b_mp.Excursion));
plot(time_hov,  G_hov_OL  - 100,     'k--', 'LineWidth', 1.8, 'DisplayName', sprintf('Hovorka Open-Loop (\\Delta: +%.1f mg/dL)', m_h_ol.Excursion));
plot(time_hov,  G_hov_MP  - 100,     'g-',  'LineWidth', 2.2, 'DisplayName', sprintf('Hovorka MPC Preview (\\Delta: +%.1f mg/dL)', m_h_mp.Excursion));
yline(0, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Basal Baseline (\Delta = 0)');
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('\Delta Glucose Excursion from Baseline (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Net Postprandial Glucose Excursion (\\DeltaG): Bergman vs Hovorka (%.1fg CHO)', meal_carbs_g), ...
    'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on; xlim([0, t_end]);
save_clean(fig2, fullfile(fig_dir, 'hovorka_vs_bergman_excursion_overlay.png'));
close(fig2);

% --- FIGURE 3: Clinical Metrics Bar Chart Comparison ---
fig3 = figure('Color', 'w', 'Position', [120, 100, 1000, 550], 'Visible', 'off');
categories = {'Open-Loop', 'MPC No-Preview', 'MPC Meal-Preview'};
peak_data = [m_b_ol.Peak, m_h_ol.Peak; m_b_np.Peak, m_h_np.Peak; m_b_mp.Peak, m_h_mp.Peak];
tir_data  = [m_b_ol.TIR,  m_h_ol.TIR;  m_b_np.TIR,  m_h_np.TIR;  m_b_mp.TIR,  m_h_mp.TIR];

subplot(1, 2, 1);
b1 = bar(categorical(categories, categories), peak_data);
b1(1).FaceColor = [0.2 0.5 0.8]; % Bergman blue
b1(2).FaceColor = [0.8 0.3 0.2]; % Hovorka red
ylabel('Peak Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Peak Glucose Comparison (%.1fg CHO)', meal_carbs_g), 'FontSize', 12, 'FontWeight', 'bold');
legend({'Bergman (3-State)', 'Hovorka (11-State)'}, 'Location', 'northeast');
grid on; yline(180, 'r--', 'Hyperglycemia Limit (180)', 'LineWidth', 1.2);

subplot(1, 2, 2);
b2 = bar(categorical(categories, categories), tir_data);
b2(1).FaceColor = [0.2 0.5 0.8];
b2(2).FaceColor = [0.8 0.3 0.2];
ylabel('Time in Range 70-180 (%)', 'FontSize', 11, 'FontWeight', 'bold');
title('Time in Range (TIR %) Comparison', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Bergman (3-State)', 'Hovorka (11-State)'}, 'Location', 'southeast');
grid on; ylim([0, 110]); yline(70, 'g--', 'Clinical Goal (70%)', 'LineWidth', 1.2);

save_clean(fig3, fullfile(fig_dir, 'hovorka_vs_bergman_metrics_barchart.png'));
close(fig3);

% Save full comparison dataset
mat_path = fullfile(data_dir, 'hovorka_vs_bergman_data.mat');
results = struct();
results.meal_carbs_g = meal_carbs_g;
results.meal_dur     = meal_dur;
results.meal_start_time = meal_start_time;
results.t_end        = t_end;
results.bolus_U      = bolus_U;
results.bolus_speed_Uh = bolus_speed_Uh;
results.bergman = struct('time', time_berg, 'G_OL', G_berg_OL, 'G_NP', G_berg_NP, 'G_MP', G_berg_MP, ...
    'u_OL', u_berg_OL, 'u_NP', u_berg_NP, 'u_MP', u_berg_MP, 'm_OL', m_b_ol, 'm_NP', m_b_np, 'm_MP', m_b_mp);
results.hovorka = struct('time', time_hov, 'G_OL', G_hov_OL, 'G_NP', G_hov_NP, 'G_MP', G_hov_MP, ...
    'u_OL', u_hov_OL, 'u_NP', u_hov_NP, 'u_MP', u_hov_MP, 'm_OL', m_h_ol, 'm_NP', m_h_np, 'm_MP', m_h_mp);

save(mat_path, 'results');
fprintf('\nComparison dataset successfully saved to: %s\n', mat_path);
fprintf('===============================================================================================\n');

metrics_table = struct('Bergman', struct('OL', m_b_ol, 'NP', m_b_np, 'MP', m_b_mp), ...
                       'Hovorka', struct('OL', m_h_ol, 'NP', m_h_np, 'MP', m_h_mp));
end

%% Helper to cleanly overwrite figures
function clean_save_helper(fig_handle, file_path)
    if exist(file_path, 'file')
        try, delete(file_path); catch, end
    end
    try
        exportgraphics(fig_handle, file_path, 'Resolution', 300);
    catch
        saveas(fig_handle, file_path);
    end
    fprintf('  Saved Figure: %s\n', file_path);
end
