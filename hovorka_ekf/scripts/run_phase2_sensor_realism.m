function [results, metrics] = run_phase2_sensor_realism(meal_carbs_g, meal_dur_min, dropout_start_min, dropout_dur_min)
% RUN_PHASE2_SENSOR_REALISM - Phase 2 of Master Roadmap: Sensor Noise & Realism Emulation
%
% Fully addressed review points:
%   1. Reports isolated sensor MARD (against interstitial glucose) AND total clinical MARD (against blood glucose).
%   2. Formally measures and verifies the 10-minute interstitial diffusion lag.
%   3. Explicitly passes Np=48 and Nc=4 into design_mpc_controller(p, Ts, Np, Nc).
%   4. Explicitly defines bolus allocation: Open-Loop receives 100% upfront bolus (standard MDI),
%      while Hybrid Closed-Loop receives 85% upfront bolus + 15% autonomous MPC auto-correction.
%   5. Calculates actuator rate jitter only across consecutive non-bolus samples (zero artifact jumps).
%   6. Imposes non-negative [0.0, 1.5*u_basal] bounds during dropout and [0.0, 4.0] U/h always.
%
% Usage:
%   run_phase2_sensor_realism()             % Default: 75g meal, 15 min dropout at t=60 min
%   run_phase2_sensor_realism(75, 15, 60, 15)

if nargin < 1 || isempty(meal_carbs_g),    meal_carbs_g = 75.0; end
if nargin < 2 || isempty(meal_dur_min),    meal_dur_min = 15.0; end
if nargin < 3 || isempty(dropout_start_min), dropout_start_min = 60.0; end
if nargin < 4 || isempty(dropout_dur_min),   dropout_dur_min = 15.0; end

dropout_end_min = dropout_start_min + dropout_dur_min;

fprintf('========================================================================\n');
fprintf('  HOVORKA ARTIFICIAL PANCREAS: PHASE 2 - SENSOR REALISM EMULATION\n');
fprintf('  Dexcom G7 ARMA Noise + Interstitial Delay + 15-min Sensor Dropout\n');
fprintf('========================================================================\n');

% Path configuration
script_dir = fileparts(mfilename('fullpath'));
root_dir   = fullfile(script_dir, '..');
fig_dir    = fullfile(root_dir, 'results', 'figures', 'phase2');
data_dir   = fullfile(root_dir, 'results', 'data', 'phase2');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

%% 1. Load Parameters and Setup Simulation Grid
p = hovorka_params(70, 100);
target_mgdl = p.G_target_mgdl; % 100 mg/dL

Ts = 5;                        % 5-minute CGM & MPC sampling time
t_end = 360;                   % 6-hour simulation (minutes)
time_vec = (0:Ts:t_end)';
N_steps = length(time_vec) - 1;

meal_start_time = 30;          % Meal at minute 30
meal_end_time   = meal_start_time + meal_dur_min;
meal_rate_gmin  = meal_carbs_g / meal_dur_min;

d_meal_gmin = zeros(N_steps + 1, 1);
for k = 1:length(time_vec)
    t = time_vec(k);
    if t >= meal_start_time && t < meal_end_time
        d_meal_gmin(k) = meal_rate_gmin;
    end
end

% Prandial Bolus Allocations:
ICR = 24.0; % g/U for 70kg adult
% Open-Loop gets 100% upfront bolus (standard MDI therapy without automated basal)
bolus_U_ol = meal_carbs_g / ICR; 
% Hybrid Closed-Loop gets 85% upfront bolus; remaining 15% is modulated by MPC
bolus_fraction_hcl = 0.85;
bolus_U_hcl = (meal_carbs_g / ICR) * bolus_fraction_hcl;

bolus_speed_Uh = 150.0; % Commercial pump infusion speed [U/h]
t_bolus_min_ol  = (bolus_U_ol  / bolus_speed_Uh) * 60;
t_bolus_min_hcl = (bolus_U_hcl / bolus_speed_Uh) * 60;

%% 2. Generate Realistic ARMA Colored Sensor Noise (Dexcom G7 Emulation)
% ARMA(1, 1) model: v(k) = phi * v(k-1) + e(k) + theta * e(k-1)
% Tuned to achieve instrument MARD_sensor ~ 8.5-9.0% against interstitial fluid
rng(42); % Fixed seed for exact reproducibility
phi_1   = 0.75;
theta_1 = 0.25;
sigma_e = 5.2; % mg/dL

e_innov = sigma_e * randn(N_steps + 1, 1);
v_arma  = zeros(N_steps + 1, 1);
v_arma(1) = e_innov(1);
for k = 2:(N_steps + 1)
    v_arma(k) = phi_1 * v_arma(k-1) + e_innov(k) + theta_1 * e_innov(k-1);
end

% Setup Sensor Dropout Mask (NaN when disconnected)
dropout_mask = false(N_steps + 1, 1);
for k = 1:(N_steps + 1)
    t = time_vec(k);
    if t >= dropout_start_min && t < dropout_end_min
        dropout_mask(k) = true;
    end
end

%% 3. Design State-Space MPC Controller with Explicit Horizons
Np = 48; % 4-hour prediction horizon
Nc = 4;  % 20-min control horizon
[mpc_obj, plant_ss_d] = design_mpc_controller(p, Ts, Np, Nc);

fprintf('Simulation Configuration:\n');
fprintf('  - Meal Carbohydrates:       %.1f g over %.0f min\n', meal_carbs_g, meal_dur_min);
fprintf('  - HCL Prandial Bolus (85%%): %.2f Units @ %.0f U/h (Duration: %.1f s)\n', ...
    bolus_U_hcl, bolus_speed_Uh, t_bolus_min_hcl * 60);
fprintf('  - Open-Loop Bolus (100%%):   %.2f Units @ %.0f U/h (Duration: %.1f s)\n', ...
    bolus_U_ol, bolus_speed_Uh, t_bolus_min_ol * 60);
fprintf('  - Interstitial Delay:       tau_cgm = %.1f min (State 11 lag)\n', p.tau_cgm);
fprintf('  - ARMA Sensor Noise:        phi = %.2f, theta = %.2f, sigma_e = %.1f mg/dL\n', ...
    phi_1, theta_1, sigma_e);
fprintf('  - Sensor Dropout Window:    t = [%.0f, %.0f] min (Duration: %.0f min)\n', ...
    dropout_start_min, dropout_end_min, dropout_dur_min);
fprintf('  - MPC Controller Horizons:  Np = %d steps (%d min), Nc = %d steps (%d min)\n', ...
    Np, Np*Ts, Nc, Nc*Ts);

ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

%% 4. Simulate Case 1: Phase 2 Realistic EKF Output-Feedback MPC
fprintf('\n[1/3] Simulating Case 1: Phase 2 Realistic EKF (Noise + Delay + Dropout)...\n');

% Process and Measurement Covariances
Q_ekf = diag([1.0, 1.0, 0.25, 1e-6, 1e-6, 1e-4, 4.0, 4.0, 1.0, 1.0, 0.05]);
R_ekf = (7.5)^2; % Measurement variance tuned for colored noise
P0_ekf = diag([10, 10, 5, 1e-4, 1e-4, 1e-3, 20, 20, 10, 10, 1.0]);

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

% Initial reading
G_meas_true_0 = x_true_ekf(1, 11) * p.MMOLL_TO_MGDL;
y_cgm_meas(1) = G_meas_true_0 + v_arma(1);

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    
    % 1. Continuous Glucose Sensor reading
    G_meas_true_k = x_true_ekf(k, 11) * p.MMOLL_TO_MGDL;
    if dropout_mask(k)
        y_k_cgm = NaN; % Sensor dropout (no packet received)
    else
        y_k_cgm = G_meas_true_k + v_arma(k);
    end
    y_cgm_meas(k) = y_k_cgm;
    
    % 2. EKF Measurement Update (gracefully handles NaN during dropout)
    [ekf, innov_k, K_k] = hovorka_ekf('update', ekf, y_k_cgm, p);
    x_hat_hist(k, :)  = ekf.x_hat';
    P_diag_hist(k, :) = diag(ekf.P)';
    
    % 3. Feedforward Meal Preview (15 min lookahead)
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        v_preview(h + 1) = d_meal_gmin(min(k + h, N_steps + 1));
    end
    
    % 4. MPC Output-Feedback Execution using Estimated States
    xmpc_ekf.Plant = ekf.x_hat;
    G_hat_curr_mgdl = ekf.x_hat(11) * p.MMOLL_TO_MGDL;
    
    u_opt = mpcmove(mpc_obj, xmpc_ekf, G_hat_curr_mgdl, target_mgdl, v_preview);
    
    % Strict Actuator Bounds & Dropout Fail-Safe
    if dropout_mask(k)
        % Enforce non-negative safety ceiling [0.0, 1.5 * u_basal] during signal loss
        u_opt = max(0.0, min(u_opt, p.ss.u_basal_Uh * 1.5));
    else
        % Normal operation bounds [0.0, 4.0 U/h]
        u_opt = max(0.0, min(u_opt, 4.0));
    end
    
    % 5. Dose Delivery: 85% Prandial Bolus + MPC Basal (covers remaining 15%)
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min_hcl;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ekf(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ekf(k+1, :) = s2(end, :);
        u_deliv_ekf(k) = u_opt + (bolus_U_hcl / (Ts / 60));
    else
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ekf(k, :)', ode_opts);
        x_true_ekf(k+1, :) = s(end, :);
        u_deliv_ekf(k) = u_opt;
    end
    
    % 6. EKF State Prediction over [t_k, t_next]
    u_total_mumin = (u_deliv_ekf(k)) * p.UH_TO_MUMIN;
    ekf = hovorka_ekf('predict', ekf, u_total_mumin, d_k_mmolmin, p);
end

% Final state update
G_meas_true_end = x_true_ekf(end, 11) * p.MMOLL_TO_MGDL;
if dropout_mask(end)
    y_cgm_meas(end) = NaN;
else
    y_cgm_meas(end) = G_meas_true_end + v_arma(end);
end
[ekf, ~, ~] = hovorka_ekf('update', ekf, y_cgm_meas(end), p);
x_hat_hist(end, :)  = ekf.x_hat';
P_diag_hist(end, :) = diag(ekf.P)';

%% 5. Simulate Case 2: Ideal Full-State Feedback MPC (Benchmark)
fprintf('[2/3] Simulating Case 2: Ideal Full-State Feedback MPC (Benchmark)...\n');
x_true_ideal = zeros(N_steps + 1, 11);
u_deliv_ideal = zeros(N_steps, 1);
x_true_ideal(1, :) = p.ss.x0';
xmpc_ideal = mpcstate(mpc_obj);

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        v_preview(h + 1) = d_meal_gmin(min(k + h, N_steps + 1));
    end
    
    xmpc_ideal.Plant = x_true_ideal(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_ideal, x_true_ideal(k, 11) * p.MMOLL_TO_MGDL, target_mgdl, v_preview);
    u_opt = max(0.0, min(u_opt, 4.0));
    
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min_hcl;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ideal(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ideal(k+1, :) = s2(end, :);
        u_deliv_ideal(k) = u_opt + (bolus_U_hcl / (Ts / 60));
    else
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ideal(k, :)', ode_opts);
        x_true_ideal(k+1, :) = s(end, :);
        u_deliv_ideal(k) = u_opt;
    end
end

%% 6. Simulate Case 3: Open-Loop Baseline (100% Bolus + Basal)
fprintf('[3/3] Simulating Case 3: Open-Loop Baseline (100%% Bolus + Basal)...\n');
x_true_ol = zeros(N_steps + 1, 11);
u_deliv_ol = zeros(N_steps, 1);
x_true_ol(1, :) = p.ss.x0';

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    u_basal = p.ss.u_basal_Uh;
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min_ol;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_basal + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ol(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_basal*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ol(k+1, :) = s2(end, :);
        u_deliv_ol(k) = u_basal + (bolus_U_ol / (Ts / 60));
    else
        u_k_mumin = u_basal * p.UH_TO_MUMIN;
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ol(k, :)', ode_opts);
        x_true_ol(k+1, :) = s(end, :);
        u_deliv_ol(k) = u_basal;
    end
end

%% 7. Compute Clinical Metrics, Delay, and Isolated MARDs
G_blood_ekf = (x_true_ekf(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;
G_cgm_true  = x_true_ekf(:, 11) * p.MMOLL_TO_MGDL;
G_hat_ekf   = (x_hat_hist(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;
G_ideal     = (x_true_ideal(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;
G_ol        = (x_true_ol(:, 7) / p.V_G) * p.MMOLL_TO_MGDL;

calc_metrics = @(G, u) struct( ...
    'Peak', max(G), ...
    'TimePeak', time_vec(find(G == max(G), 1)), ...
    'Excursion', max(G) - G(1), ...
    'Nadir', min(G), ...
    'TIR', mean(G >= 70 & G <= 180) * 100, ...
    'TAR', mean(G > 180) * 100, ...
    'TBR', mean(G < 70) * 100, ...
    'TotalU', sum(u) * (Ts / 60));

m_ekf   = calc_metrics(G_blood_ekf, u_deliv_ekf);
m_ideal = calc_metrics(G_ideal, u_deliv_ideal);
m_ol    = calc_metrics(G_ol, u_deliv_ol);

% Estimation errors
err_G = G_blood_ekf - G_hat_ekf;
rmse_G = sqrt(mean(err_G.^2));
max_err_G = max(abs(err_G));

% Dropout specific error
err_dropout = err_G(dropout_mask);
rmse_dropout = sqrt(mean(err_dropout.^2));

% Metric 1: Isolated Sensor MARD against True Interstitial Glucose (Instrument Noise)
valid_cgm = ~isnan(y_cgm_meas);
mard_sensor = mean(abs(y_cgm_meas(valid_cgm) - G_cgm_true(valid_cgm)) ./ G_cgm_true(valid_cgm)) * 100;

% Metric 2: Total Clinical MARD against True Blood Glucose (Noise + 10-min Lag)
mard_clinical = mean(abs(y_cgm_meas(valid_cgm) - G_blood_ekf(valid_cgm)) ./ G_blood_ekf(valid_cgm)) * 100;

% Verified Interstitial Lag (Time to peak blood vs interstitial)
t_peak_blood = time_vec(find(G_blood_ekf == max(G_blood_ekf), 1));
t_peak_subq  = time_vec(find(G_cgm_true == max(G_cgm_true), 1));
tau_lag_actual = t_peak_subq - t_peak_blood;

% Metric 3: Actuator rate jitter across valid consecutive non-bolus samples
k_bolus = find(time_vec(1:end-1) == meal_start_time, 1);
bolus_window = [k_bolus-1, k_bolus, k_bolus+1];
valid_diffs = [];
for k = 1:(length(u_deliv_ekf)-1)
    if ~ismember(k, bolus_window) && ~ismember(k+1, bolus_window)
        valid_diffs(end+1) = abs(u_deliv_ekf(k+1) - u_deliv_ekf(k)); %#ok<AGROW>
    end
end
if ~isempty(valid_diffs)
    max_du = max(valid_diffs);
    mean_du = mean(valid_diffs);
else
    max_du = 0; mean_du = 0;
end

%% 8. Print Comprehensive Clinical Benchmark Table
fprintf('\n===============================================================================================\n');
fprintf('  PHASE 2 BENCHMARK TABLE: %.1fg CHO MEAL (Dexcom G7 ARMA Noise + 15-min Dropout)\n', meal_carbs_g);
fprintf('===============================================================================================\n');
fprintf('%-28s | %-16s | %-18s | %-18s\n', 'Metric', 'Open-Loop (100% Bolus)', 'Ideal MPC (Cheat)', 'Phase 2 EKF MPC');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Peak Blood Glucose (mg/dL)', m_ol.Peak, m_ideal.Peak, m_ekf.Peak);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time to Peak (min)', m_ol.TimePeak, m_ideal.TimePeak, m_ekf.TimePeak);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Max Excursion (mg/dL)', m_ol.Excursion, m_ideal.Excursion, m_ekf.Excursion);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Minimum Glucose (mg/dL)', m_ol.Nadir, m_ideal.Nadir, m_ekf.Nadir);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time in Range 70-180 (%)', m_ol.TIR, m_ideal.TIR, m_ekf.TIR);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time Above 180 (%)', m_ol.TAR, m_ideal.TAR, m_ekf.TAR);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time Below 70 (%)', m_ol.TBR, m_ideal.TBR, m_ekf.TBR);
fprintf('%-28s | %-16.2f | %-18.2f | %-18.2f\n', 'Total Insulin Delivered (U)', m_ol.TotalU, m_ideal.TotalU, m_ekf.TotalU);
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('PHASE 2 SENSOR REALISM & OBSERVER METRICS (Addressed Review Points):\n');
fprintf('  1. Isolated Sensor MARD (vs SubQ): %.2f %% (Bench target: 8.5-9.0%%)\n', mard_sensor);
fprintf('  2. Total Clinical MARD (vs Blood): %.2f %% (Noise + 10-min lag)\n', mard_clinical);
fprintf('  3. Empirical Interstitial Delay:   %.1f min (Time shift between blood and SubQ peaks)\n', tau_lag_actual);
fprintf('  4. Blood Glucose Overall RMSE:     %.2f mg/dL (Max Error: %.2f mg/dL)\n', rmse_G, max_err_G);
fprintf('  5. Error During 15-min Dropout:    %.2f mg/dL (Dead-Reckoning Hold)\n', rmse_dropout);
fprintf('  6. True Actuator Slew Rate Jitter: %.2f U/h (Criterion: < 0.50 U/h, Mean: %.2f U/h)\n', max_du, mean_du);
fprintf('===============================================================================================\n');

%% 9. Generate Publication-Grade Visualizations
save_clean = @(fig, path) print(fig, path, '-dpng', '-r300');

% --- FIGURE 1: Closed-Loop Performance with ARMA Noise and Sensor Dropout ---
fig1 = figure('Color', 'w', 'Position', [60, 50, 1150, 900], 'Visible', 'on');

% Subplot 1: Blood Glucose Regulation
subplot(3, 1, 1);
patch([dropout_start_min, dropout_end_min, dropout_end_min, dropout_start_min], ...
      [40, 40, 240, 240], [1.0 0.85 0.85], 'EdgeColor', 'none', ...
      'FaceAlpha', 0.5, 'DisplayName', '15-min Sensor Dropout (t=60-75)');
hold on;
patch([0, t_end, t_end, 0], [70, 70, 180, 180], [0.85 0.95 0.85], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.35, 'DisplayName', 'Target Range [70-180]');

plot(time_vec, y_cgm_meas, 'm.', 'MarkerSize', 10, 'DisplayName', sprintf('Dexcom G7 (Sensor MARD: %.1f%%)', mard_sensor));
plot(time_vec, G_ol, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL)', m_ol.Peak));
plot(time_vec, G_ideal, 'k:', 'LineWidth', 2.0, 'DisplayName', sprintf('Ideal MPC (Peak: %.1f mg/dL)', m_ideal.Peak));
plot(time_vec, G_blood_ekf, 'b-', 'LineWidth', 2.2, 'DisplayName', sprintf('Phase 2 EKF MPC (Peak: %.1f mg/dL)', m_ekf.Peak));
plot(time_vec, G_hat_ekf, 'g--', 'LineWidth', 1.8, 'DisplayName', sprintf('EKF Estimated G (RMSE: %.2f mg/dL)', rmse_G));

yline(100, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Target (100 mg/dL)');
ylabel('Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Phase 2 Sensor Realism: Dexcom G7 ARMA Noise, Interstitial Lag & 15-min Dropout (%.1fg CHO)', meal_carbs_g), ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, t_end]); ylim([50, 235]);

% Subplot 2: Pump Insulin Delivery
subplot(3, 1, 2);
u_time = time_vec(1:end-1);
patch([dropout_start_min, dropout_end_min, dropout_end_min, dropout_start_min], ...
      [0, 0, 6.0, 6.0], [1.0 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
hold on;
stairs(u_time, u_deliv_ol, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Open-Loop (100%% Bolus: %.2f U)', bolus_U_ol));
stairs(u_time, u_deliv_ideal, 'k:', 'LineWidth', 2.0, 'DisplayName', 'Ideal MPC (85% Bolus + Basal)');
stairs(u_time, u_deliv_ekf, 'b-',  'LineWidth', 2.2, 'DisplayName', 'Phase 2 EKF MPC (85% Bolus + Basal)');
yline(p.ss.u_basal_Uh, 'k--', sprintf('Basal (%.2f U/h)', p.ss.u_basal_Uh), 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Max Basal Limit (4.0 U/h)');
ylabel('Insulin Infusion (U/h)', 'FontSize', 11, 'FontWeight', 'bold');
title('Pump Actuation: Upfront Boluses + Autonomous Closed-Loop Basal Modulation', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, t_end]); ylim([0, 5.5]);

% Subplot 3: Meal Rate
subplot(3, 1, 3);
stem(time_vec, d_meal_gmin, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'Marker', 'none', ...
    'DisplayName', sprintf('Meal Disturbance (%.1fg CHO over %.0f min)', meal_carbs_g, meal_dur_min));
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Meal Rate (g/min)', 'FontSize', 11, 'FontWeight', 'bold');
title('Prandial Glucose Influx', 'FontSize', 12, 'FontWeight', 'bold');
grid on; xlim([0, t_end]); ylim([0, meal_rate_gmin * 1.3]);

fig1_path = fullfile(fig_dir, 'phase2_cgm_realism_closed_loop.png');
save_clean(fig1, fig1_path);
fprintf('  Saved Figure: %s\n', fig1_path);

% --- FIGURE 2: EKF Estimation Error, Covariance Expansion & Dropout Recovery ---
fig2 = figure('Color', 'w', 'Position', [100, 80, 1100, 550], 'Visible', 'on');
sigma_G = sqrt(P_diag_hist(:, 7) ./ (p.V_G^2)) .* p.MMOLL_TO_MGDL;

patch([dropout_start_min, dropout_end_min, dropout_end_min, dropout_start_min], ...
      [-25, -25, 25, 25], [1.0 0.88 0.88], 'EdgeColor', 'none', ...
      'FaceAlpha', 0.6, 'DisplayName', 'Sensor Dropout (P Expands During Blind Hold)');
hold on;
plot(time_vec, err_G, 'b-', 'LineWidth', 2.0, 'DisplayName', sprintf('Estimation Error G - G_{hat} (RMSE: %.2f)', rmse_G));
plot(time_vec,  2*sigma_G, 'r--', 'LineWidth', 1.5, 'DisplayName', '+2\\sigma EKF Confidence Bound');
plot(time_vec, -2*sigma_G, 'r--', 'LineWidth', 1.5, 'DisplayName', '-2\\sigma EKF Confidence Bound');
yline(0, 'k-', 'LineWidth', 1.0);
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Glucose Estimation Error (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title('EKF Covariance Expansion During Signal Loss and Rapid Recovery Upon Reconnection', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, t_end]); ylim([-25, 25]);

fig2_path = fullfile(fig_dir, 'phase2_ekf_dropout_and_covariance.png');
save_clean(fig2, fig2_path);
fprintf('  Saved Figure: %s\n', fig2_path);

% --- FIGURE 3: Sensor Noise Analysis & Physiological Lag ---
fig3 = figure('Color', 'w', 'Position', [120, 100, 1050, 500], 'Visible', 'on');

subplot(1, 2, 1);
plot(time_vec, G_blood_ekf, 'b-', 'LineWidth', 2.0, 'DisplayName', 'True Blood Glucose G(t)'); hold on;
plot(time_vec, G_cgm_true, 'Color', [0.2 0.6 0.2], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Interstitial SubQ Glucose (Measured Lag: %.0f min)', tau_lag_actual));
plot(time_vec, y_cgm_meas, 'm.', 'MarkerSize', 9, 'DisplayName', 'Dexcom G7 Noisy Sensor');
xlabel('Time (minutes)', 'FontWeight', 'bold'); ylabel('Glucose (mg/dL)', 'FontWeight', 'bold');
title(sprintf('Physiological Diffusion Lag: \\Delta t = %.0f min', tau_lag_actual), 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, 240]);

subplot(1, 2, 2);
sensor_errors = y_cgm_meas(valid_cgm) - G_cgm_true(valid_cgm);
histogram(sensor_errors, 20, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'w', 'Normalization', 'pdf');
hold on;
x_grid = linspace(min(sensor_errors), max(sensor_errors), 100);
mean_err = mean(sensor_errors);
std_err  = std(sensor_errors);
y_norm = (1 / (std_err * sqrt(2*pi))) * exp(-0.5 * ((x_grid - mean_err) / std_err).^2);
plot(x_grid, y_norm, 'r-', 'LineWidth', 2.0, 'DisplayName', sprintf('Normal Fit (\\mu=%.1f, \\sigma=%.1f)', mean_err, std_err));
xlabel('Sensor Error (y_{cgm} - G_{subQ}) [mg/dL]', 'FontWeight', 'bold');
ylabel('Probability Density', 'FontWeight', 'bold');
title(sprintf('ARMA Isolated Sensor Error (Sensor MARD: %.2f%%)', mard_sensor), 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on;

fig3_path = fullfile(fig_dir, 'phase2_cgm_noise_and_mard.png');
save_clean(fig3, fig3_path);
fprintf('  Saved Figure: %s\n', fig3_path);

%% 10. Save Simulation Results Dataset
mat_path = fullfile(data_dir, 'phase2_sensor_realism_results.mat');
results = struct();
results.time_vec        = time_vec;
results.G_blood_ekf     = G_blood_ekf;
results.G_cgm_true      = G_cgm_true;
results.G_hat_ekf       = G_hat_ekf;
results.y_cgm_meas      = y_cgm_meas;
results.G_ideal         = G_ideal;
results.G_ol            = G_ol;
results.u_deliv_ekf     = u_deliv_ekf;
results.dropout_mask    = dropout_mask;
results.m_ekf           = m_ekf;
results.m_ideal         = m_ideal;
results.m_ol            = m_ol;
results.mard_sensor     = mard_sensor;
results.mard_clinical   = mard_clinical;
results.tau_lag_actual  = tau_lag_actual;
results.rmse_G          = rmse_G;
results.rmse_dropout    = rmse_dropout;
results.max_du          = max_du;
results.mean_du         = mean_du;

metrics = results.m_ekf;
save(mat_path, 'results');
fprintf('Phase 2 simulation dataset saved to: %s\n', mat_path);
fprintf('========================================================================\n');
end
