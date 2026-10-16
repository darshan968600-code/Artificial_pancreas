function [metrics, results] = run_ekf_closed_loop(meal_carbs_g, meal_dur, meal_start_time, t_end, cgm_noise_sigma)
% RUN_EKF_CLOSED_LOOP Hardware-Ready Realistic Artificial Pancreas Simulation.
% Implements an Extended Kalman Filter (EKF) to reconstruct all 11 unmeasured
% physiological states (subcutaneous depots S1, S2, plasma insulin I, and
% remote receptor actions x1, x2, x3) from a single NOISY Continuous Glucose
% Monitor (CGM) sensor signal, driving the Model Predictive Controller (MPC).
%
% Compares:
%   1. EKF Output-Feedback MPC (Realistic: Noisy CGM -> EKF Observer -> MPC)
%   2. Ideal Full-State Feedback MPC ("Omniscient" theoretical benchmark)
%   3. Open-Loop Baseline
%
% Usage:
%   run_ekf_closed_loop()                    -> Interactive prompts (defaults: 75g, 15 min)
%   run_ekf_closed_loop(75)                  -> 75g CHO over 15 min at t=60min
%   run_ekf_closed_loop(100, 20)             -> 100g CHO over 20 min
%   run_ekf_closed_loop(75, 15, 60, 360, 5)  -> Custom timeline & noise std = 5 mg/dL

clc;
fprintf('========================================================================\n');
fprintf('  HOVORKA ARTIFICIAL PANCREAS: EXTENDED KALMAN FILTER (EKF) CLOSED LOOP\n');
fprintf('  Realistic Observer-Based Output-Feedback Control (Hardware Ready)\n');
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
    meal_start_time = 60; % min
end

if nargin < 4 || isempty(t_end)
    t_end = 360; % min (6 hours)
end

if nargin < 5 || isempty(cgm_noise_sigma)
    cgm_noise_sigma = 5.0; % mg/dL (realistic Dexcom G7 / Abbott Libre 3 noise)
end

%% 2. Setup Paths & Directories
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts')
    project_root = fileparts(this_dir);
else
    project_root = this_dir;
end
addpath(fullfile(project_root, 'scripts'));

fig_dir  = fullfile(project_root, 'results', 'figures', 'phase1');
data_dir = fullfile(project_root, 'results', 'data', 'phase1');
if ~exist(fig_dir, 'dir'),  mkdir(fig_dir);  end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

save_clean = @(fh, fp) clean_save_helper(fh, fp);
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

%% 3. Load Physiology & Configure MPC + EKF
p = hovorka_params(70, 100);
[mpc_obj, plant_ss_d] = design_mpc_controller();

Ts = 5; % 5-minute control & sensing interval
N_steps = t_end / Ts;
time_vec = (0:N_steps)' * Ts;

% Disturbance array
meal_rate_gmin = meal_carbs_g / meal_dur;
d_meal_gmin = zeros(N_steps + 1, 1);
for k = 1:(N_steps + 1)
    t_curr = time_vec(k);
    if t_curr >= meal_start_time && t_curr < (meal_start_time + meal_dur)
        d_meal_gmin(k) = meal_rate_gmin;
    end
end

% Commercial HCL Prandial Bolus
ICR = 24.0; % g/U
bolus_fraction = 0.85;
bolus_U = (meal_carbs_g / ICR) * bolus_fraction;
bolus_speed_Uh = 150.0; % U/h motor speed
t_bolus_min = (bolus_U / bolus_speed_Uh) * 60;

% Initialize EKF Observer
R_cov = cgm_noise_sigma^2; % Measurement noise covariance
ekf = hovorka_ekf('init', [], p, Ts, [], R_cov, p.ss.x0);

fprintf('\nSimulation Scenario:\n');
fprintf('  - Meal Carbohydrates:   %.1f grams (Duration: %.0f min)\n', meal_carbs_g, meal_dur);
fprintf('  - Prandial Bolus:       %.2f Units @ %.1f U/h (Duration: %.1f s)\n', bolus_U, bolus_speed_Uh, t_bolus_min*60);
fprintf('  - Sensor Noise Level:   sigma = %.1f mg/dL (R = %.1f)\n', cgm_noise_sigma, R_cov);
fprintf('  - MPC Control Interval: Ts = %d min\n', Ts);

%% 4. Simulate Case 1: EKF-Driven Output-Feedback MPC (Realistic Closed-Loop)
fprintf('\n[1/3] Simulating Case 1: Realistic EKF Output-Feedback MPC...\n');

x_true_ekf  = zeros(N_steps + 1, 11);
x_hat_hist  = zeros(N_steps + 1, 11);
P_diag_hist = zeros(N_steps + 1, 11);
u_deliv_ekf = zeros(N_steps, 1);
y_cgm_noisy = zeros(N_steps + 1, 1);

x_true_ekf(1, :) = p.ss.x0';
x_hat_hist(1, :) = ekf.x_hat';
P_diag_hist(1, :) = diag(ekf.P)';

% Seed reproducible sensor noise
rng(42);
v_noise = cgm_noise_sigma * randn(N_steps + 1, 1);
y_cgm_noisy(1) = x_true_ekf(1, 11) * p.MMOLL_TO_MGDL + v_noise(1);

xmpc_ekf = mpcstate(mpc_obj);
Np = mpc_obj.PredictionHorizon;
target_mgdl = 100;

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    
    % 1. Sensor Measurement arrives from CGM
    y_k_cgm = y_cgm_noisy(k);
    
    % 2. EKF Measurement Update Step: correct state estimate using noisy CGM
    [ekf, innov_k, K_k] = hovorka_ekf('update', ekf, y_k_cgm, p);
    x_hat_hist(k, :) = ekf.x_hat';
    P_diag_hist(k, :) = diag(ekf.P)';
    
    % 3. Build Preview Sequence
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        v_preview(h + 1) = d_meal_gmin(min(k + h, N_steps + 1));
    end
    
    % 4. MPC computes optimal command using ONLY EKF ESTIMATED STATES
    xmpc_ekf.Plant = ekf.x_hat; % Feeds estimated states, NOT true states!
    u_opt = mpcmove(mpc_obj, xmpc_ekf, ekf.x_hat(11) * p.MMOLL_TO_MGDL, target_mgdl, v_preview);
    
    % 5. Apply Insulin & Prandial Bolus to TRUE Physical Plant
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ekf(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ekf(k+1, :) = s2(end, :);
        u_deliv_ekf(k) = u_opt + (bolus_U / (Ts / 60));
        u_step_avg_mumin = (u_deliv_ekf(k)) * p.UH_TO_MUMIN;
    else
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ekf(k, :)', ode_opts);
        x_true_ekf(k+1, :) = s(end, :);
        u_deliv_ekf(k) = u_opt;
        u_step_avg_mumin = u_k_mumin;
    end
    
    % Generate next noisy CGM reading
    y_cgm_noisy(k+1) = x_true_ekf(k+1, 11) * p.MMOLL_TO_MGDL + v_noise(k+1);
    
    % 6. EKF State Prediction Step forward to t_{k+1}
    ekf = hovorka_ekf('predict', ekf, u_step_avg_mumin, d_k_mmolmin, p);
end
x_hat_hist(end, :) = ekf.x_hat';

%% 5. Simulate Case 2: Ideal Full-State Feedback MPC ("Omniscient" Benchmark)
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
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_opt + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ideal(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_opt*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ideal(k+1, :) = s2(end, :);
        u_deliv_ideal(k) = u_opt + (bolus_U / (Ts / 60));
    else
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ideal(k, :)', ode_opts);
        x_true_ideal(k+1, :) = s(end, :);
        u_deliv_ideal(k) = u_opt;
    end
end

%% 6. Simulate Case 3: Open-Loop Baseline
fprintf('[3/3] Simulating Case 3: Open-Loop Baseline (Bolus + Basal)...\n');

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
        t_split = t_k + t_bolus_min;
        [~, s1] = ode15s(@(t, x) hovorka_ode(t, x, (u_basal + bolus_speed_Uh)*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_split], x_true_ol(k, :)', ode_opts);
        [~, s2] = ode15s(@(t, x) hovorka_ode(t, x, u_basal*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_split, t_next], s1(end, :)', ode_opts);
        x_true_ol(k+1, :) = s2(end, :);
        u_deliv_ol(k) = u_basal + (bolus_U / (Ts / 60));
    else
        [~, s] = ode15s(@(t, x) hovorka_ode(t, x, u_basal*p.UH_TO_MUMIN, d_k_mmolmin, p), ...
            [t_k, t_next], x_true_ol(k, :)', ode_opts);
        x_true_ol(k+1, :) = s(end, :);
        u_deliv_ol(k) = u_basal;
    end
end

%% 7. Extract Signals & Calculate Clinical & Estimation Metrics
G_ekf_true  = (x_true_ekf(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_ekf_hat   = (x_hat_hist(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_ideal     = (x_true_ideal(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_ol        = (x_true_ol(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;

% Estimation Error Metrics
G_err = G_ekf_true - G_ekf_hat;
rmse_G = sqrt(mean(G_err.^2));
max_err_G = max(abs(G_err));

% Subcutaneous S1, S2 Error
S1_err = x_true_ekf(:, 1) - x_hat_hist(:, 1);
S2_err = x_true_ekf(:, 2) - x_hat_hist(:, 2);
I_err  = x_true_ekf(:, 3) - x_hat_hist(:, 3);
rmse_S1 = sqrt(mean(S1_err.^2));
rmse_S2 = sqrt(mean(S2_err.^2));
rmse_I  = sqrt(mean(I_err.^2));

calc_met = @(G, u) struct( ...
    'Peak', max(G), ...
    't_Peak', time_vec(find(G == max(G), 1)), ...
    'Excursion', max(G) - 100, ...
    'Min', min(G), ...
    'TIR', 100 * mean(G >= 70 & G <= 180), ...
    'TAR', 100 * mean(G > 180), ...
    'TBR', 100 * mean(G < 70), ...
    'TotalU', sum(u) * (Ts / 60) ...
);

m_ekf   = calc_met(G_ekf_true, u_deliv_ekf);
m_ideal = calc_met(G_ideal, u_deliv_ideal);
m_ol    = calc_met(G_ol, u_deliv_ol);

%% 8. Display Results Summary Table
fprintf('\n===============================================================================================\n');
fprintf('  HOVORKA EKF CLINICAL BENCHMARK TABLE: %.1fg CHO MEAL (CGM Noise \\sigma = %.1f mg/dL)\n', ...
    meal_carbs_g, cgm_noise_sigma);
fprintf('===============================================================================================\n');
fprintf('%-28s | %-16s | %-18s | %-18s\n', 'Metric', 'Open-Loop', 'Ideal MPC (Cheat)', 'EKF MPC (Realistic)');
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Peak Blood Glucose (mg/dL)', m_ol.Peak, m_ideal.Peak, m_ekf.Peak);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time to Peak (min)', m_ol.t_Peak, m_ideal.t_Peak, m_ekf.t_Peak);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Max Excursion (mg/dL)', m_ol.Excursion, m_ideal.Excursion, m_ekf.Excursion);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Minimum Glucose (mg/dL)', m_ol.Min, m_ideal.Min, m_ekf.Min);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time in Range 70-180 (%)', m_ol.TIR, m_ideal.TIR, m_ekf.TIR);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time Above 180 (%)', m_ol.TAR, m_ideal.TAR, m_ekf.TAR);
fprintf('%-28s | %-16.1f | %-18.1f | %-18.1f\n', 'Time Below 70 (%)', m_ol.TBR, m_ideal.TBR, m_ekf.TBR);
fprintf('%-28s | %-16.2f | %-18.2f | %-18.2f\n', 'Total Insulin Delivered (U)', m_ol.TotalU, m_ideal.TotalU, m_ekf.TotalU);
fprintf('-----------------------------------------------------------------------------------------------\n');
fprintf('EKF ESTIMATION ACCURACY METRICS (Reconstructing 11 Hidden States from Noisy CGM):\n');
fprintf('  - Blood Glucose Estimation RMSE:   %.2f mg/dL (Max Abs Error: %.2f mg/dL)\n', rmse_G, max_err_G);
fprintf('  - Subcutaneous Depot S1 RMSE:      %.2f mU\n', rmse_S1);
fprintf('  - Subcutaneous Depot S2 RMSE:      %.2f mU\n', rmse_S2);
fprintf('  - Blood Plasma Insulin I RMSE:     %.2f mU/L\n', rmse_I);
fprintf('===============================================================================================\n');

%% =====================================================================
%% 9. GENERATE PUBLICATION-GRADE EKF FIGURES
%% =====================================================================

% --- FIGURE 1: Closed-Loop Performance (True vs Estimated vs Ideal vs Noisy CGM) ---
fig1 = figure('Color', 'w', 'Position', [80, 60, 1100, 850], 'Visible', 'on');

% Subplot 1: Glucose Regulation
subplot(3, 1, 1);
plot(time_vec, y_cgm_noisy, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 8, 'DisplayName', 'Noisy CGM Sensor (\sigma=5 mg/dL)');
hold on;
plot(time_vec, G_ol, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL)', m_ol.Peak));
plot(time_vec, G_ideal, 'k:', 'LineWidth', 2.0, 'DisplayName', sprintf('Ideal Full-State MPC (Peak: %.1f mg/dL)', m_ideal.Peak));
plot(time_vec, G_ekf_true, 'b-', 'LineWidth', 2.2, 'DisplayName', sprintf('EKF MPC True Glucose (Peak: %.1f mg/dL)', m_ekf.Peak));
plot(time_vec, G_ekf_hat, 'g--', 'LineWidth', 1.8, 'DisplayName', sprintf('EKF Estimated G (RMSE: %.2f mg/dL)', rmse_G));
yline(100, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Target (100 mg/dL)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180 mg/dL)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70 mg/dL)');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Target Range [70-180]');
ylabel('Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Realistic Output-Feedback Control: Hovorka EKF vs Ideal MPC (%.1fg CHO Meal)', meal_carbs_g), ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, t_end]); ylim([50, max(240, ceil(m_ol.Peak/20)*20 + 20)]);

% Subplot 2: Insulin Infusion Rate
subplot(3, 1, 2);
u_time = time_vec(1:end-1);
stairs(u_time, u_deliv_ol, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Open-Loop');
hold on;
stairs(u_time, u_deliv_ideal, 'k:', 'LineWidth', 2.0, 'DisplayName', 'Ideal Full-State MPC');
stairs(u_time, u_deliv_ekf, 'b-',  'LineWidth', 2.2, 'DisplayName', 'EKF Output-Feedback MPC');
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'MPC Basal Limit (4.0 U/h)');
ylabel('Insulin Infusion (U/h)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Insulin Pump Infusion (Prandial Bolus: %.2f U @ 150 U/h + Safe MPC Basal)', bolus_U), ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on; xlim([0, t_end]); ylim([0, 6.0]);

% Subplot 3: Meal Disturbance
subplot(3, 1, 3);
stem(time_vec, d_meal_gmin, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'Marker', 'none', ...
    'DisplayName', sprintf('Custom Meal (%.1fg CHO over %.0f min)', meal_carbs_g, meal_dur));
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Meal Rate (g/min)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Prandial Glucose Disturbance (%.1fg CHO at %.2f g/min)', meal_carbs_g, meal_rate_gmin), ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on; xlim([0, t_end]);

fig1_path = fullfile(fig_dir, 'ekf_closed_loop_performance.png');
save_clean(fig1, fig1_path);

% --- FIGURE 2: Hidden State Tracking (Proof that EKF reconstructs unmeasured biology) ---
fig2 = figure('Color', 'w', 'Position', [100, 80, 1100, 800], 'Visible', 'on');

% Subplot 1: Subcutaneous Depot S1 & S2
subplot(2, 2, 1);
plot(time_vec, x_true_ekf(:, 1), 'b-', 'LineWidth', 2.0, 'DisplayName', 'True S1'); hold on;
plot(time_vec, x_hat_hist(:, 1), 'b--', 'LineWidth', 1.8, 'DisplayName', 'Estimated S1');
plot(time_vec, x_true_ekf(:, 2), 'r-', 'LineWidth', 2.0, 'DisplayName', 'True S2');
plot(time_vec, x_hat_hist(:, 2), 'r--', 'LineWidth', 1.8, 'DisplayName', 'Estimated S2');
ylabel('Insulin Mass (mU)', 'FontWeight', 'bold');
title('Subcutaneous Insulin Depots (S1, S2)', 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, t_end]);

% Subplot 2: Plasma Insulin Concentration I
subplot(2, 2, 2);
plot(time_vec, x_true_ekf(:, 3), 'k-', 'LineWidth', 2.0, 'DisplayName', 'True Plasma Insulin (I)'); hold on;
plot(time_vec, x_hat_hist(:, 3), 'm--', 'LineWidth', 2.0, 'DisplayName', 'Estimated Insulin I (mU/L)');
ylabel('Plasma Insulin (mU/L)', 'FontWeight', 'bold');
title(sprintf('Plasma Insulin Concentration (RMSE: %.2f mU/L)', rmse_I), 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, t_end]);

% Subplot 3: Remote Insulin Actions (x1, x2, x3)
subplot(2, 2, 3);
plot(time_vec, x_true_ekf(:, 4)*1e3, 'b-', 'LineWidth', 1.8, 'DisplayName', 'True x1 \times 10^3'); hold on;
plot(time_vec, x_hat_hist(:, 4)*1e3, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Estimated x1');
plot(time_vec, x_true_ekf(:, 6), 'Color', [0.4 0.2 0.6], 'LineWidth', 1.8, 'DisplayName', 'True x3 (EGP)');
plot(time_vec, x_hat_hist(:, 6), '--', 'Color', [0.4 0.2 0.6], 'LineWidth', 1.5, 'DisplayName', 'Estimated x3');
xlabel('Time (min)', 'FontWeight', 'bold'); ylabel('Remote Action Level', 'FontWeight', 'bold');
title('Remote Insulin Actions (x1 Transport, x3 EGP)', 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, t_end]);

% Subplot 4: Glucose Compartments (Q1, Q2)
subplot(2, 2, 4);
plot(time_vec, x_true_ekf(:, 7), 'k-', 'LineWidth', 2.0, 'DisplayName', 'True Q1 (Accessible)'); hold on;
plot(time_vec, x_hat_hist(:, 7), 'g--', 'LineWidth', 1.8, 'DisplayName', 'Estimated Q1');
plot(time_vec, x_true_ekf(:, 8), 'Color', [0.8 0.4 0.1], 'LineWidth', 2.0, 'DisplayName', 'True Q2 (Peripheral)');
plot(time_vec, x_hat_hist(:, 8), '--', 'Color', [0.8 0.4 0.1], 'LineWidth', 1.8, 'DisplayName', 'Estimated Q2');
xlabel('Time (min)', 'FontWeight', 'bold'); ylabel('Glucose Mass (mmol)', 'FontWeight', 'bold');
title('Tissue Glucose Masses (Q1 Accessible, Q2 Non-Accessible)', 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, t_end]);

save_clean(fig2, fullfile(fig_dir, 'ekf_hidden_states_tracking.png'));

% --- FIGURE 3: Estimation Error & Covariance Bounds ---
fig3 = figure('Color', 'w', 'Position', [120, 100, 1000, 500], 'Visible', 'on');
sigma_G = sqrt(P_diag_hist(:, 7) ./ (p.V_G^2)) .* p.MMOLL_TO_MGDL;
plot(time_vec, G_err, 'b-', 'LineWidth', 1.8, 'DisplayName', sprintf('Estimation Error G - G_{hat} (RMSE: %.2f)', rmse_G));
hold on;
plot(time_vec,  2*sigma_G, 'r--', 'LineWidth', 1.2, 'DisplayName', '+2\sigma EKF Confidence Bound');
plot(time_vec, -2*sigma_G, 'r--', 'LineWidth', 1.2, 'DisplayName', '-2\sigma EKF Confidence Bound');
yline(0, 'k-', 'LineWidth', 1.0);
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Glucose Estimation Error (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title('Extended Kalman Filter Estimation Error with 95% (\pm 2\sigma) Covariance Bounds', ...
    'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast'); grid on; xlim([0, t_end]);
ylim([-max(15, ceil(max(abs(G_err)))), max(15, ceil(max(abs(G_err))))]);

save_clean(fig3, fullfile(fig_dir, 'ekf_estimation_error_and_cgm.png'));

% Save full simulation dataset
mat_path = fullfile(data_dir, 'ekf_simulation_results.mat');
results = struct();
results.meal_carbs_g = meal_carbs_g;
results.meal_dur     = meal_dur;
results.cgm_noise_sigma = cgm_noise_sigma;
results.time         = time_vec;
results.G_ekf_true   = G_ekf_true;
results.G_ekf_hat    = G_ekf_hat;
results.G_ideal      = G_ideal;
results.G_ol         = G_ol;
results.y_cgm_noisy  = y_cgm_noisy;
results.x_true       = x_true_ekf;
results.x_hat        = x_hat_hist;
results.P_diag       = P_diag_hist;
results.u_ekf        = u_deliv_ekf;
results.u_ideal      = u_deliv_ideal;
results.u_ol         = u_deliv_ol;
results.m_ekf        = m_ekf;
results.m_ideal      = m_ideal;
results.m_ol         = m_ol;
results.rmse_G       = rmse_G;
results.rmse_S1      = rmse_S1;
results.rmse_S2      = rmse_S2;
results.rmse_I       = rmse_I;

save(mat_path, 'results');
fprintf('\nEKF Simulation dataset successfully saved to: %s\n', mat_path);
fprintf('===============================================================================================\n');

metrics = struct('EKF', m_ekf, 'Ideal', m_ideal, 'OpenLoop', m_ol, 'RMSE_G', rmse_G);
end

%% Clean save helper
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
