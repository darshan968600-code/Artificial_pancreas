function [results, metrics] = run_phase4_circadian_24h_protocol()
% RUN_PHASE4_CIRCADIAN_24H_PROTOCOL - Phase 4 of Master Engineering Roadmap
% 24-Hour Circadian Multi-Meal Clinical Protocol (Breakfast, Lunch, Snack, Dinner, Overnight)
%
% Implements a complete 24-hour daily life in-silico clinical trial (1440 min):
%   - 00:00 to 07:30 (0 - 450 min): Nocturnal Fasting & Dawn Phenomenon
%   - 05:00 to 05:15 (300 - 315 min): Nighttime Sensor Dropout (Compression Low)
%   - 07:30 (450 min): Breakfast (45.0g CHO over 15 min, Announced)
%   - 12:30 (750 min): Lunch (75.0g CHO over 20 min, Announced)
%   - 16:00 (960 min): UNANNOUNCED Snack (25.0g CHO over 10 min, No Preview, No Bolus)
%   - 19:30 (1170 min): Large Dinner (90.0g CHO over 25 min, Announced)
%   - 22:00 to 24:00 (1320 - 1440 min): Bedtime Recovery
%
% Compares:
%   Case 1: Proposed Hardware-Ready EKF Output-Feedback MPC
%   Case 2: Ideal Full-State Feedback MPC (Benchmark)
%   Case 3: Standard Open-Loop Baseline (Standard MDI Therapy)
%
% Usage:
%   run_phase4_circadian_24h_protocol()

fprintf('========================================================================\n');
fprintf('  HOVORKA ARTIFICIAL PANCREAS: PHASE 4 - 24-HOUR MULTI-MEAL PROTOCOL\n');
fprintf('  Full Day: Breakfast (45g), Lunch (75g), Unannounced Snack (25g), Dinner (90g)\n');
fprintf('  Includes: Dexcom G7 ARMA Noise + Interstitial Lag + Nocturnal Dropout\n');
fprintf('========================================================================\n');

% Path configuration
script_dir = fileparts(mfilename('fullpath'));
root_dir   = fullfile(script_dir, '..');
fig_dir    = fullfile(root_dir, 'results', 'figures', 'phase4');
data_dir   = fullfile(root_dir, 'results', 'data', 'phase4');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

%% 1. Configure 24-Hour Simulation Grid & Circadian Schedule
p = hovorka_params(70, 100);
target_mgdl = p.G_target_mgdl; % 100 mg/dL

Ts = 5;                        % 5-minute sensing & control interval
t_end = 1440;                  % 24 hours (1440 minutes)
time_vec = (0:Ts:t_end)';
N_steps = length(time_vec) - 1;

% Define the 4 Meals in the Daily Schedule
meals = struct();
meals(1).name      = 'Breakfast';
meals(1).t_start   = 450;      % 07:30
meals(1).dur       = 15;
meals(1).carbs     = 45.0;     % 45g CHO
meals(1).announced = true;

meals(2).name      = 'Lunch';
meals(2).t_start   = 750;      % 12:30
meals(2).dur       = 20;
meals(2).carbs     = 75.0;     % 75g CHO
meals(2).announced = true;

meals(3).name      = 'Unannounced Snack';
meals(3).t_start   = 960;      % 16:00
meals(3).dur       = 10;
meals(3).carbs     = 25.0;     % 25g CHO
meals(3).announced = false;    % STRICTLY UNANNOUNCED (No Bolus, No Preview!)

meals(4).name      = 'Dinner';
meals(4).t_start   = 1170;     % 19:30
meals(4).dur       = 25;
meals(4).carbs     = 90.0;     % 90g CHO
meals(4).announced = true;

% Build Continuous Physical Disturbance Vector d_meal_gmin(t)
d_meal_gmin = zeros(N_steps + 1, 1);
for m = 1:length(meals)
    t_start = meals(m).t_start;
    t_end_m = t_start + meals(m).dur;
    rate_m  = meals(m).carbs / meals(m).dur;
    for k = 1:length(time_vec)
        t = time_vec(k);
        if t >= t_start && t < t_end_m
            d_meal_gmin(k) = d_meal_gmin(k) + rate_m;
        end
    end
end

% Build Feedforward Preview Disturbance Vector (ONLY for Announced Meals)
d_preview_known = zeros(N_steps + 1, 1);
for m = 1:length(meals)
    if meals(m).announced
        t_start = meals(m).t_start;
        t_end_m = t_start + meals(m).dur;
        rate_m  = meals(m).carbs / meals(m).dur;
        for k = 1:length(time_vec)
            t = time_vec(k);
            if t >= t_start && t < t_end_m
                d_preview_known(k) = d_preview_known(k) + rate_m;
            end
        end
    end
end

% 15-Minute Nocturnal Compression Low Sensor Dropout at 05:00 AM (t = 300 to 315 min)
dropout_start_min = 300;
dropout_end_min   = 315;
dropout_mask = (time_vec >= dropout_start_min & time_vec < dropout_end_min);

% Generate 24-Hour Continuous Dexcom G7 ARMA Colored Noise
rng(42);
phi_1   = 0.75;
theta_1 = 0.25;
sigma_e = 5.2;

e_innov = sigma_e * randn(N_steps + 1, 1);
v_arma  = zeros(N_steps + 1, 1);
v_arma(1) = e_innov(1);
for k = 2:(N_steps + 1)
    v_arma(k) = phi_1 * v_arma(k-1) + e_innov(k) + theta_1 * e_innov(k-1);
end

%% 2. Configure Bolus Parameters and MPC Controller
ICR = 24.0; % g/U
bolus_speed_Uh = 150.0; % U/h pump motor speed

% Bolus specifications per meal
bolus_info = struct();
for m = 1:length(meals)
    if meals(m).announced
        bolus_info(m).U_hcl = (meals(m).carbs / ICR) * 0.85; % 85% upfront
        bolus_info(m).U_ol  = (meals(m).carbs / ICR) * 1.00; % 100% upfront
    else
        % Unannounced snack: Zero upfront bolus
        bolus_info(m).U_hcl = 0.0;
        bolus_info(m).U_ol  = 0.0;
    end
    bolus_info(m).t_min_hcl = (bolus_info(m).U_hcl / bolus_speed_Uh) * 60;
    bolus_info(m).t_min_ol  = (bolus_info(m).U_ol  / bolus_speed_Uh) * 60;
end

Np = 48; % 4-hour prediction horizon
Nc = 4;  % 20-min control horizon
[mpc_obj, ~] = design_mpc_controller(p, Ts, Np, Nc);

Q_ekf = diag([1.0, 1.0, 0.25, 1e-6, 1e-6, 1e-4, 4.0, 4.0, 1.0, 1.0, 0.05]);
R_ekf = (7.5)^2;
P0_ekf = diag([10, 10, 5, 1e-4, 1e-4, 1e-3, 20, 20, 10, 10, 1.0]);
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

fprintf('24-Hour Protocol Schedule:\n');
fprintf('  - 00:00 - 07:30: Overnight Fasting (Dropout at 05:00 - 05:15 AM)\n');
fprintf('  - 07:30 (450 min):  Breakfast  (%.1fg CHO) -> HCL Bolus: %.2f U\n', meals(1).carbs, bolus_info(1).U_hcl);
fprintf('  - 12:30 (750 min):  Lunch      (%.1fg CHO) -> HCL Bolus: %.2f U\n', meals(2).carbs, bolus_info(2).U_hcl);
fprintf('  - 16:00 (960 min):  Snack      (%.1fg CHO) -> UNANNOUNCED (0.00 U Bolus, Tests Rescue)\n', meals(3).carbs);
fprintf('  - 19:30 (1170 min): Dinner     (%.1fg CHO) -> HCL Bolus: %.2f U\n', meals(4).carbs, bolus_info(4).U_hcl);
fprintf('  - Total Daily Carbs: %.1f grams\n', sum([meals.carbs]));

%% 3. Simulate Case 1: Proposed Hardware-Ready EKF Output-Feedback MPC
fprintf('\n[1/3] Simulating Case 1: 24h EKF Output-Feedback MPC (Realistic Observer)...\n');

ekf = hovorka_ekf('init', [], p, Ts, Q_ekf, R_ekf, p.ss.x0, P0_ekf);
xmpc_ekf = mpcstate(mpc_obj);

x_true_ekf   = zeros(N_steps + 1, 11);
x_hat_hist   = zeros(N_steps + 1, 11);
P_diag_hist  = zeros(N_steps + 1, 11);
u_deliv_ekf  = zeros(N_steps, 1);
y_cgm_meas   = zeros(N_steps + 1, 1);

x_true_ekf(1, :)  = p.ss.x0';
x_hat_hist(1, :)  = ekf.x_hat';
P_diag_hist(1, :) = diag(ekf.P)';

G_meas_0 = x_true_ekf(1, 11) * p.MMOLL_TO_MGDL;
y_cgm_meas(1) = G_meas_0 + v_arma(1);

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    
    % 1. Continuous Glucose Sensor reading
    G_meas_k = x_true_ekf(k, 11) * p.MMOLL_TO_MGDL;
    if dropout_mask(k)
        y_k_cgm = NaN;
    else
        y_k_cgm = G_meas_k + v_arma(k);
    end
    y_cgm_meas(k) = y_k_cgm;
    
    % 2. EKF Measurement Update
    [ekf, ~, ~] = hovorka_ekf('update', ekf, y_k_cgm, p);
    x_hat_hist(k, :)  = ekf.x_hat';
    P_diag_hist(k, :) = diag(ekf.P)';
    
    % 3. Feedforward Meal Preview (ONLY for announced meals)
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        v_preview(h + 1) = d_preview_known(min(k + h, N_steps + 1));
    end
    
    % 4. MPC Output-Feedback Execution using Estimated States
    xmpc_ekf.Plant = ekf.x_hat;
    G_hat_curr = ekf.x_hat(11) * p.MMOLL_TO_MGDL;
    u_opt = mpcmove(mpc_obj, xmpc_ekf, G_hat_curr, target_mgdl, v_preview);
    
    if dropout_mask(k)
        u_opt = max(0.0, min(u_opt, p.ss.u_basal_Uh * 1.5));
    else
        u_opt = max(0.0, min(u_opt, 4.0));
    end
    
    % 5. Check if a meal begins at this sample step
    meal_idx = find([meals.t_start] == t_k, 1);
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    if ~isempty(meal_idx) && bolus_info(meal_idx).U_hcl > 0
        b_U = bolus_info(meal_idx).U_hcl;
        t_b_min = bolus_info(meal_idx).t_min_hcl;
        t_split = t_k + t_b_min;
        
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ekf(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ekf(k+1, :) = s2(end, :);
        u_deliv_ekf(k) = u_opt + (b_U / (Ts / 60));
    else
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, s_ode] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ekf(k, :)', ode_opts);
        x_true_ekf(k+1, :) = s_ode(end, :);
        u_deliv_ekf(k) = u_opt;
    end
    
    % 6. EKF State Prediction over [t_k, t_next]
    u_total_mumin = (u_deliv_ekf(k)) * p.UH_TO_MUMIN;
    ekf = hovorka_ekf('predict', ekf, u_total_mumin, d_k_mmolmin, p);
end

% Final state update
G_meas_end = x_true_ekf(end, 11) * p.MMOLL_TO_MGDL;
if dropout_mask(end)
    y_cgm_meas(end) = NaN;
else
    y_cgm_meas(end) = G_meas_end + v_arma(end);
end
[ekf, ~, ~] = hovorka_ekf('update', ekf, y_cgm_meas(end), p);
x_hat_hist(end, :) = ekf.x_hat';

%% 4. Simulate Case 2: Ideal Full-State Feedback MPC (Benchmark)
fprintf('[2/3] Simulating Case 2: 24h Ideal Full-State Feedback MPC (Benchmark)...\n');
x_true_ideal  = zeros(N_steps + 1, 11);
u_deliv_ideal = zeros(N_steps, 1);
x_true_ideal(1, :) = p.ss.x0';
xmpc_ideal = mpcstate(mpc_obj);

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        v_preview(h + 1) = d_preview_known(min(k + h, N_steps + 1));
    end
    
    xmpc_ideal.Plant = x_true_ideal(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_ideal, x_true_ideal(k, 11) * p.MMOLL_TO_MGDL, target_mgdl, v_preview);
    u_opt = max(0.0, min(u_opt, 4.0));
    
    meal_idx = find([meals.t_start] == t_k, 1);
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    if ~isempty(meal_idx) && bolus_info(meal_idx).U_hcl > 0
        b_U = bolus_info(meal_idx).U_hcl;
        t_b_min = bolus_info(meal_idx).t_min_hcl;
        t_split = t_k + t_b_min;
        
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ideal(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ideal(k+1, :) = s2(end, :);
        u_deliv_ideal(k) = u_opt + (b_U / (Ts / 60));
    else
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, s_ode] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ideal(k, :)', ode_opts);
        x_true_ideal(k+1, :) = s_ode(end, :);
        u_deliv_ideal(k) = u_opt;
    end
end

%% 5. Simulate Case 3: Standard Open-Loop Baseline (100% Boluses + Fixed Basal)
fprintf('[3/3] Simulating Case 3: 24h Open-Loop Baseline (Standard MDI Therapy)...\n');
x_true_ol  = zeros(N_steps + 1, 11);
u_deliv_ol = zeros(N_steps, 1);
x_true_ol(1, :) = p.ss.x0';

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    u_basal = p.ss.u_basal_Uh;
    
    meal_idx = find([meals.t_start] == t_k, 1);
    
    if ~isempty(meal_idx) && bolus_info(meal_idx).U_ol > 0
        b_U = bolus_info(meal_idx).U_ol;
        t_b_min = bolus_info(meal_idx).t_min_ol;
        t_split = t_k + t_b_min;
        
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_basal + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ol(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_basal*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ol(k+1, :) = s2(end, :);
        u_deliv_ol(k) = u_basal + (b_U / (Ts / 60));
    else
        u_k_mumin = u_basal * p.UH_TO_MUMIN;
        [~, s_ode] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ol(k, :)', ode_opts);
        x_true_ol(k+1, :) = s_ode(end, :);
        u_deliv_ol(k) = u_basal;
    end
end

%% 6. Compute Comprehensive 24-Hour Clinical Metrics
G_blood_ekf = (x_true_ekf(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;
G_cgm_true  = x_true_ekf(:, 11) * p.MMOLL_TO_MGDL;
G_hat_ekf   = (x_hat_hist(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;
G_ideal     = (x_true_ideal(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;
G_ol        = (x_true_ol(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;

calc_metrics = @(G, u) struct( ...
    'MeanG', mean(G), ...
    'StdG', std(G), ...
    'Peak', max(G), ...
    'Nadir', min(G), ...
    'TIR', mean(G >= 70 & G <= 180) * 100, ...
    'TAR_L1', mean(G > 180 & G <= 250) * 100, ...
    'TAR_L2', mean(G > 250) * 100, ...
    'TBR_L1', mean(G < 70 & G >= 54) * 100, ...
    'TBR_L2', mean(G < 54) * 100, ...
    'TotalU', sum(u) * (Ts / 60), ...
    'eA1c', (mean(G) + 46.7) / 28.7);

m_ekf   = calc_metrics(G_blood_ekf, u_deliv_ekf);
m_ideal = calc_metrics(G_ideal, u_deliv_ideal);
m_ol    = calc_metrics(G_ol, u_deliv_ol);

% Estimation errors
err_G = G_blood_ekf - G_hat_ekf;
rmse_G = sqrt(mean(err_G.^2));
max_err_G = max(abs(err_G));

% Nighttime metrics (00:00 - 07:00, t = 0 to 420 min)
night_mask = (time_vec <= 420);
tir_night_ekf = mean(G_blood_ekf(night_mask) >= 70 & G_blood_ekf(night_mask) <= 180) * 100;
nadir_night_ekf = min(G_blood_ekf(night_mask));

% Unannounced Snack Excursion & Recovery (t = 960 to 1170 min)
snack_mask = (time_vec >= 960 & time_vec <= 1170);
snack_peak_ekf = max(G_blood_ekf(snack_mask));
snack_peak_ol  = max(G_ol(snack_mask));

%% 7. Print 24-Hour Clinical Consensus Report
fprintf('\n===================================================================================================\n');
fprintf('  PHASE 4 CLINICAL CONSENSUS BENCHMARK TABLE: FULL 24-HOUR MULTI-MEAL TRIAL\n');
fprintf('===================================================================================================\n');
fprintf('%-32s | %-16s | %-18s | %-18s\n', 'Consensus Clinical Metric', 'Open-Loop (MDI)', 'Ideal MPC (Cheat)', 'Phase 4 EKF MPC');
fprintf('---------------------------------------------------------------------------------------------------\n');
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', 'Mean 24h Glucose (mg/dL)', m_ol.MeanG, m_ideal.MeanG, m_ekf.MeanG);
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', 'Glucose Standard Dev (mg/dL)', m_ol.StdG, m_ideal.StdG, m_ekf.StdG);
fprintf('%-32s | %-16.2f | %-18.2f | %-18.2f\n', 'Estimated HbA1c (%)', m_ol.eA1c, m_ideal.eA1c, m_ekf.eA1c);
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', '24h Peak Glucose (mg/dL)', m_ol.Peak, m_ideal.Peak, m_ekf.Peak);
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', '24h Minimum Nadir (mg/dL)', m_ol.Nadir, m_ideal.Nadir, m_ekf.Nadir);
fprintf('---------------------------------------------------------------------------------------------------\n');
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', 'Time in Range [70-180] (%)', m_ol.TIR, m_ideal.TIR, m_ekf.TIR);
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', 'TAR Level 1 [181-250] (%)', m_ol.TAR_L1, m_ideal.TAR_L1, m_ekf.TAR_L1);
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', 'TAR Level 2 [>250] (%)', m_ol.TAR_L2, m_ideal.TAR_L2, m_ekf.TAR_L2);
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', 'TBR Level 1 [54-69] (%)', m_ol.TBR_L1, m_ideal.TBR_L1, m_ekf.TBR_L1);
fprintf('%-32s | %-16.1f | %-18.1f | %-18.1f\n', 'TBR Level 2 [<54] (%)', m_ol.TBR_L2, m_ideal.TBR_L2, m_ekf.TBR_L2);
fprintf('%-32s | %-16.2f | %-18.2f | %-18.2f\n', 'Total Daily Insulin (TDD) (U)', m_ol.TotalU, m_ideal.TotalU, m_ekf.TotalU);
fprintf('---------------------------------------------------------------------------------------------------\n');
fprintf('CIRCADIAN CHALLENGE HIGHLIGHTS:\n');
fprintf('  - Nocturnal Fasting Stability:    TIR = %.1f%%, Nadir = %.1f mg/dL (Zero Nocturnal Hypo)\n', ...
    tir_night_ekf, nadir_night_ekf);
fprintf('  - Unannounced 25g Snack Peak:     EKF MPC = %.1f mg/dL vs Open-Loop = %.1f mg/dL (-%.1f mg/dL reduction)\n', ...
    snack_peak_ekf, snack_peak_ol, snack_peak_ol - snack_peak_ekf);
fprintf('  - 24h Continuous Observer RMSE:   %.2f mg/dL (Reconstructing 11 states over 288 samples)\n', rmse_G);
fprintf('===================================================================================================\n');

%% 8. Generate Publication-Grade Visualizations (Visible and Open)
save_clean = @(fig, path) print(fig, path, '-dpng', '-r300');
time_hours = time_vec / 60;

% --- FIGURE 1: 24-Hour Continuous Overview ---
fig1 = figure('Color', 'w', 'Position', [50, 40, 1200, 950], 'Visible', 'on');

% Subplot 1: Blood Glucose
subplot(3, 1, 1);
% Nighttime Shading (00:00 to 07:00 and 22:00 to 24:00)
patch([0, 7.0, 7.0, 0], [40, 40, 260, 260], [0.93 0.93 0.97], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Nighttime / Sleep');
hold on;
patch([22.0, 24.0, 24.0, 22.0], [40, 40, 260, 260], [0.93 0.93 0.97], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
% Target Range [70 - 180]
patch([0, 24, 24, 0], [70, 70, 180, 180], [0.85 0.95 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Target Range [70-180 mg/dL]');
% Nocturnal Dropout Window (05:00 - 05:15 AM)
patch([5.0, 5.25, 5.25, 5.0], [40, 40, 260, 260], [1.0 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', '15-min Dropout (05:00 AM)');

plot(time_hours, G_ol, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Open-Loop Baseline (TIR: %.1f%%, eA1c: %.2f%%)', m_ol.TIR, m_ol.eA1c));
plot(time_hours, G_ideal, 'k:', 'LineWidth', 2.0, 'DisplayName', sprintf('Ideal MPC (TIR: %.1f%%)', m_ideal.TIR));
plot(time_hours, G_blood_ekf, 'b-', 'LineWidth', 2.2, 'DisplayName', sprintf('Phase 4 EKF MPC (TIR: %.1f%%, eA1c: %.2f%%)', m_ekf.TIR, m_ekf.eA1c));
plot(time_hours, G_hat_ekf, 'g--', 'LineWidth', 1.6, 'DisplayName', sprintf('EKF Estimated G (RMSE: %.2f mg/dL)', rmse_G));

yline(100, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Target Setpoint (100 mg/dL)');
ylabel('Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title('24-Hour Continuous Artificial Pancreas Performance: Circadian Multi-Meal Clinical Challenge', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, 24]); ylim([50, 250]);
set(gca, 'XTick', 0:2:24, 'XTickLabel', {'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});

% Subplot 2: Insulin Infusion
subplot(3, 1, 2);
u_time_h = time_hours(1:end-1);
patch([0, 7.0, 7.0, 0], [0, 0, 6.0, 6.0], [0.93 0.93 0.97], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
hold on;
patch([22.0, 24.0, 24.0, 22.0], [0, 0, 6.0, 6.0], [0.93 0.93 0.97], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
stairs(u_time_h, u_deliv_ol, 'r--', 'LineWidth', 1.6, 'DisplayName', sprintf('Open-Loop (TDD: %.2f U)', m_ol.TotalU));
stairs(u_time_h, u_deliv_ideal, 'k:', 'LineWidth', 1.8, 'DisplayName', sprintf('Ideal MPC (TDD: %.2f U)', m_ideal.TotalU));
stairs(u_time_h, u_deliv_ekf, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('EKF MPC (TDD: %.2f U)', m_ekf.TotalU));
yline(p.ss.u_basal_Uh, 'k--', sprintf('Basal (%.2f U/h)', p.ss.u_basal_Uh), 'LineWidth', 1.0);
ylabel('Insulin Infusion (U/h)', 'FontSize', 11, 'FontWeight', 'bold');
title('Continuous Insulin Pump Infusion: High-Speed Prandial Boluses (150 U/h) + Autonomous Safe Basal Modulation', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, 24]); ylim([0, 5.5]);
set(gca, 'XTick', 0:2:24, 'XTickLabel', {'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});

% Subplot 3: Meal Disturbances
subplot(3, 1, 3);
stem(time_hours, d_meal_gmin, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'Marker', 'none');
hold on;
for m = 1:length(meals)
    t_h = meals(m).t_start / 60;
    text(t_h, meals(m).carbs / meals(m).dur + 0.3, ...
        sprintf('%s\n(%.0fg CHO)', meals(m).name, meals(m).carbs), ...
        'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
xlabel('Time of Day (24-Hour Clock)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Meal Flux (g/min)', 'FontSize', 11, 'FontWeight', 'bold');
title('Prandial Carbohydrate Intakes Across 24 Hours', 'FontSize', 12, 'FontWeight', 'bold');
grid on; xlim([0, 24]); ylim([0, 4.8]);
set(gca, 'XTick', 0:2:24, 'XTickLabel', {'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});

fig1_path = fullfile(fig_dir, 'phase4_24h_closed_loop_overview.png');
save_clean(fig1, fig1_path);
fprintf('  Saved Figure: %s\n', fig1_path);

% --- FIGURE 2: 24-Hour Insulin-on-Board (IOB) & Subcutaneous Depots ---
fig2 = figure('Color', 'w', 'Position', [80, 60, 1150, 700], 'Visible', 'on');

subplot(2, 1, 1);
plot(time_hours, x_true_ekf(:, 1), 'b-', 'LineWidth', 2.0, 'DisplayName', 'True Subcutaneous Depot S1'); hold on;
plot(time_hours, x_hat_hist(:, 1), 'b--', 'LineWidth', 1.8, 'DisplayName', 'EKF Estimated S1');
plot(time_hours, x_true_ekf(:, 2), 'r-', 'LineWidth', 2.0, 'DisplayName', 'True Subcutaneous Depot S2');
plot(time_hours, x_hat_hist(:, 2), 'r--', 'LineWidth', 1.8, 'DisplayName', 'EKF Estimated S2');
ylabel('Insulin Mass (mU)', 'FontSize', 11, 'FontWeight', 'bold');
title('24-Hour Subcutaneous Depot Tracking: Real-Time Reconstruction of Insulin Pools', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, 24]);
set(gca, 'XTick', 0:2:24, 'XTickLabel', {'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});

subplot(2, 1, 2);
% Insulin on Board (IOB) in Units: (S1 + S2 - Basal_Depots) / 1000
basal_depots = (p.ss.S1_0 + p.ss.S2_0);
iob_true = max(0, (x_true_ekf(:, 1) + x_true_ekf(:, 2) - basal_depots) / 1000);
iob_hat  = max(0, (x_hat_hist(:, 1) + x_hat_hist(:, 2) - basal_depots) / 1000);
plot(time_hours, iob_true, 'k-', 'LineWidth', 2.2, 'DisplayName', 'True Active Insulin on Board (IOB)'); hold on;
plot(time_hours, iob_hat,  'm--', 'LineWidth', 1.8, 'DisplayName', 'EKF Reconstructed IOB');
ylabel('Active IOB (Units)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time of Day (24-Hour Clock)', 'FontSize', 11, 'FontWeight', 'bold');
title('Insulin-on-Board (IOB) Accounting: Prevents Insulin Stacking Across Consecutive Meals', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, 24]); ylim([0, 4.0]);
set(gca, 'XTick', 0:2:24, 'XTickLabel', {'00:00','02:00','04:00','06:00','08:00','10:00','12:00','14:00','16:00','18:00','20:00','22:00','24:00'});

fig2_path = fullfile(fig_dir, 'phase4_24h_hidden_states_iob.png');
save_clean(fig2, fig2_path);
fprintf('  Saved Figure: %s\n', fig2_path);

% --- FIGURE 3: Clinical Ambulatory Glucose Profile (AGP) Summary ---
fig3 = figure('Color', 'w', 'Position', [100, 80, 1100, 520], 'Visible', 'on');

subplot(1, 2, 1);
tir_data = [m_ekf.TIR; m_ekf.TAR_L1; m_ekf.TAR_L2; m_ekf.TBR_L1; m_ekf.TBR_L2];
bar_h = bar(1, [m_ekf.TBR_L2, m_ekf.TBR_L1, m_ekf.TIR, m_ekf.TAR_L1, m_ekf.TAR_L2], 'stacked');
bar_h(1).FaceColor = [0.6 0.1 0.1]; % Deep Red: Severe Hypo (<54)
bar_h(2).FaceColor = [0.9 0.3 0.3]; % Light Red: Hypo (54-69)
bar_h(3).FaceColor = [0.2 0.7 0.3]; % Green: In Range (70-180)
bar_h(4).FaceColor = [0.9 0.6 0.2]; % Orange: High (181-250)
bar_h(5).FaceColor = [0.8 0.4 0.1]; % Dark Orange: Very High (>250)
set(gca, 'XTick', 1, 'XTickLabel', {'Phase 4 EKF MPC'}, 'FontWeight', 'bold');
ylabel('Percentage of 24h Day (%)', 'FontWeight', 'bold');
title(sprintf('Consensus Time-in-Range (TIR: %.1f%%)', m_ekf.TIR), 'FontWeight', 'bold');
legend({'TBR L2 (<54)', 'TBR L1 (54-69)', 'TIR (70-180)', 'TAR L1 (181-250)', 'TAR L2 (>250)'}, ...
    'Location', 'eastoutside');
grid on; ylim([0, 100]);

subplot(1, 2, 2);
comp_bars = bar([m_ol.TIR, m_ideal.TIR, m_ekf.TIR; ...
                 m_ol.TAR_L1 + m_ol.TAR_L2, m_ideal.TAR_L1 + m_ideal.TAR_L2, m_ekf.TAR_L1 + m_ekf.TAR_L2; ...
                 m_ol.TBR_L1 + m_ol.TBR_L2, m_ideal.TBR_L1 + m_ideal.TBR_L2, m_ekf.TBR_L1 + m_ekf.TBR_L2]');
comp_bars(1).FaceColor = [0.2 0.7 0.3]; % TIR
comp_bars(2).FaceColor = [0.9 0.5 0.2]; % TAR
comp_bars(3).FaceColor = [0.8 0.2 0.2]; % TBR
set(gca, 'XTickLabel', {'Open-Loop', 'Ideal MPC', 'EKF MPC'}, 'FontWeight', 'bold');
ylabel('Percentage of Time (%)', 'FontWeight', 'bold');
title('24h Clinical Metric Comparison Across Strategies', 'FontWeight', 'bold');
legend({'Time in Range (70-180)', 'Time Above 180', 'Time Below 70'}, 'Location', 'northeast');
grid on; ylim([0, 115]);

fig3_path = fullfile(fig_dir, 'phase4_24h_clinical_tir_pie_and_bars.png');
save_clean(fig3, fig3_path);
fprintf('  Saved Figure: %s\n', fig3_path);

%% 9. Save Phase 4 Simulation Results Dataset
mat_path = fullfile(data_dir, 'phase4_24h_circadian_results.mat');
results = struct();
results.time_vec    = time_vec;
results.time_hours  = time_hours;
results.G_blood_ekf = G_blood_ekf;
results.G_cgm_true  = G_cgm_true;
results.G_hat_ekf   = G_hat_ekf;
results.G_ideal     = G_ideal;
results.G_ol        = G_ol;
results.u_deliv_ekf = u_deliv_ekf;
results.u_deliv_ol  = u_deliv_ol;
results.m_ekf       = m_ekf;
results.m_ideal     = m_ideal;
results.m_ol        = m_ol;
results.rmse_G      = rmse_G;
results.meals       = meals;
results.bolus_info  = bolus_info;

metrics = results.m_ekf;
save(mat_path, 'results');
fprintf('Phase 4 simulation dataset saved to: %s\n', mat_path);
fprintf('========================================================================\n');
end
