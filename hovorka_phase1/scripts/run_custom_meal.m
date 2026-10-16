function [m_OL, m_NP, m_MP] = run_custom_meal(meal_carbs_g, meal_dur, meal_start_time)
% RUN_CUSTOM_MEAL Simulates the Hovorka T1D Artificial Pancreas under a custom meal.
% Mimics commercial insulin pumps delivering a prandial meal bolus at realistic
% high delivery speed (150 U/h motor delivery) while keeping the automated MPC
% basal controller safely constrained to standard basal bounds [0.0, 4.0] U/h.
%
% Usage:
%   run_custom_meal()              -> Interactively prompts for carbs & duration (defaults: 50g, 15 min)
%   run_custom_meal(50)            -> 50g CHO over 15 min at t=120min
%   run_custom_meal(100, 20)       -> 100g CHO over 20 min at t=120min
%   run_custom_meal(75, 15, 180)   -> 75g CHO over 15 min at t=180min

clc;
fprintf('========================================================================\n');
fprintf('  HOVORKA ARTIFICIAL PANCREAS: 150 U/h PRANDIAL BOLUS + MPC BASAL\n');
fprintf('========================================================================\n');

%% 1. Interactive / Default Inputs
if nargin < 1 || isempty(meal_carbs_g)
    try
        user_val = input('Enter meal carbohydrate amount in grams [Default: 50]: ', 's');
        if isempty(strtrim(user_val))
            meal_carbs_g = 50;
        else
            meal_carbs_g = str2double(user_val);
            if isnan(meal_carbs_g) || meal_carbs_g <= 0
                fprintf('  Invalid input. Defaulting to 50g CHO.\n');
                meal_carbs_g = 50;
            end
        end
    catch
        meal_carbs_g = 50;
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
    meal_start_time = 120; % min (standard 2h baseline)
end

%% 2. Setup Paths & Controller
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts')
    project_root = fileparts(this_dir);
else
    project_root = this_dir;
end
addpath(fullfile(project_root, 'scripts'));

% Design MPC controller with strict basal bounds [0.0, 4.0] U/h (unchanged)
[mpc_obj, plant_ss_d, p] = design_mpc_controller();

% Simulation Timeline
Ts = 5;               % Control interval: 5 min
t_end = 720;          % Total duration: 720 min (12 hours)
N_steps = t_end / Ts; % 144 steps
time_vec = (0:N_steps)' * Ts;

meal_rate_gmin = meal_carbs_g / meal_dur;

% Discrete meal disturbance array
d_meal_gmin = zeros(N_steps + 1, 1);
for k = 1:(N_steps + 1)
    t_curr = time_vec(k);
    if t_curr >= meal_start_time && t_curr < (meal_start_time + meal_dur)
        d_meal_gmin(k) = meal_rate_gmin;
    end
end

ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
Np = mpc_obj.PredictionHorizon;
target_mgdl = 100;

%% 3. Prandial Bolus Calculation & Delivery Speed (150 U/h)
% Patient Insulin-to-Carbohydrate Ratio (ICR): 24.0 g/U for Hovorka 70kg adult
ICR = 24.0; % g/U
bolus_fraction = 0.85; % 85% upfront prandial bolus (commercial HCL standard)
bolus_U = (meal_carbs_g / ICR) * bolus_fraction;

% Bolus motor delivery speed: 150 U/h
bolus_speed_Uh = 150.0; % U/h
% Bolus delivery duration in minutes:
t_bolus_min = (bolus_U / bolus_speed_Uh) * 60; % e.g. 0.71 min (42.5 sec) for 50g CHO

fprintf('\nSimulation Scenario:\n');
fprintf('  - Meal Carbohydrates:   %.1f grams\n', meal_carbs_g);
fprintf('  - Meal Duration:        %.1f min (Rate: %.2f g/min)\n', meal_dur, meal_rate_gmin);
fprintf('  - Prandial Bolus Dose:  %.2f Units (85%% of %.1fg @ ICR=%.0f)\n', bolus_U, meal_carbs_g, ICR);
fprintf('  - Bolus Delivery Speed: %.1f U/h (Duration: %.1f seconds)\n', bolus_speed_Uh, t_bolus_min * 60);
fprintf('  - MPC Basal Bounds:     [%.1f, %.1f] U/h (Kept safe & unchanged)\n', ...
    mpc_obj.ManipulatedVariables(1).Min, mpc_obj.ManipulatedVariables(1).Max);

%% 4. Case 1: Open-Loop Baseline (150 U/h Bolus + Fixed Basal)
fprintf('\n[1/3] Simulating Case 1: Open-Loop Baseline (Bolus + Basal)...\n');
x_OL = zeros(N_steps + 1, 11);
u_OL = zeros(N_steps, 1);
x_OL(1, :) = p.ss.x0';

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    u_basal_Uh = p.ss.u_basal_Uh;
    
    if t_k == meal_start_time
        % Stage 1: Bolus at 150 U/h for t_bolus_min
        t_split = t_k + t_bolus_min;
        u_stage1_mumin = (u_basal_Uh + bolus_speed_Uh) * p.UH_TO_MUMIN;
        [~, seg1] = ode15s(@(t, x) hovorka_ode(t, x, u_stage1_mumin, d_k_mmolmin, p), ...
            [t_k, t_split], x_OL(k, :)', ode_opts);
        
        % Stage 2: Return to basal for remainder of 5-min step
        u_stage2_mumin = u_basal_Uh * p.UH_TO_MUMIN;
        [~, seg2] = ode15s(@(t, x) hovorka_ode(t, x, u_stage2_mumin, d_k_mmolmin, p), ...
            [t_split, t_next], seg1(end, :)', ode_opts);
        
        x_OL(k+1, :) = seg2(end, :);
        u_OL(k) = u_basal_Uh + (bolus_U / (Ts / 60));
    else
        u_OL(k) = u_basal_Uh;
        u_k_mumin = u_basal_Uh * p.UH_TO_MUMIN;
        [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_OL(k, :)', ode_opts);
        x_OL(k+1, :) = x_seg(end, :);
    end
end

%% 5. Case 2: MPC Without Meal Preview (150 U/h Bolus + Reactive Basal)
fprintf('[2/3] Simulating Case 2: MPC Without Meal Preview...\n');
x_NP = zeros(N_steps + 1, 11);
u_NP = zeros(N_steps, 1);
x_NP(1, :) = p.ss.x0';

xmpc_NP = mpcstate(mpc_obj);

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    
    y_meas_k = x_NP(k, 11) * p.MMOLL_TO_MGDL;
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    xmpc_NP.Plant = x_NP(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_NP, y_meas_k, target_mgdl, d_curr);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        u_stage1_mumin = (u_opt + bolus_speed_Uh) * p.UH_TO_MUMIN;
        [~, seg1] = ode15s(@(t, x) hovorka_ode(t, x, u_stage1_mumin, d_k_mmolmin, p), ...
            [t_k, t_split], x_NP(k, :)', ode_opts);
        
        u_stage2_mumin = u_opt * p.UH_TO_MUMIN;
        [~, seg2] = ode15s(@(t, x) hovorka_ode(t, x, u_stage2_mumin, d_k_mmolmin, p), ...
            [t_split, t_next], seg1(end, :)', ode_opts);
        
        x_NP(k+1, :) = seg2(end, :);
        u_NP(k) = u_opt + (bolus_U / (Ts / 60));
    else
        u_NP(k) = u_opt;
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_NP(k, :)', ode_opts);
        x_NP(k+1, :) = x_seg(end, :);
    end
end

%% 6. Case 3: MPC With Meal Preview (150 U/h Bolus + Proactive Basal)
fprintf('[3/3] Simulating Case 3: MPC With Meal Preview...\n');
x_MP = zeros(N_steps + 1, 11);
u_MP = zeros(N_steps, 1);
x_MP(1, :) = p.ss.x0';

xmpc_MP = mpcstate(mpc_obj);

for k = 1:N_steps
    t_k = time_vec(k);
    t_next = time_vec(k+1);
    
    y_meas_k = x_MP(k, 11) * p.MMOLL_TO_MGDL;
    
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        idx_future = min(k + h, N_steps + 1);
        v_preview(h + 1) = d_meal_gmin(idx_future);
    end
    
    xmpc_MP.Plant = x_MP(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_MP, y_meas_k, target_mgdl, v_preview);
    d_curr = d_meal_gmin(k);
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    if t_k == meal_start_time
        t_split = t_k + t_bolus_min;
        u_stage1_mumin = (u_opt + bolus_speed_Uh) * p.UH_TO_MUMIN;
        [~, seg1] = ode15s(@(t, x) hovorka_ode(t, x, u_stage1_mumin, d_k_mmolmin, p), ...
            [t_k, t_split], x_MP(k, :)', ode_opts);
        
        u_stage2_mumin = u_opt * p.UH_TO_MUMIN;
        [~, seg2] = ode15s(@(t, x) hovorka_ode(t, x, u_stage2_mumin, d_k_mmolmin, p), ...
            [t_split, t_next], seg1(end, :)', ode_opts);
        
        x_MP(k+1, :) = seg2(end, :);
        u_MP(k) = u_opt + (bolus_U / (Ts / 60));
    else
        u_MP(k) = u_opt;
        u_k_mumin = u_opt * p.UH_TO_MUMIN;
        [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
            [t_k, t_next], x_MP(k, :)', ode_opts);
        x_MP(k+1, :) = x_seg(end, :);
    end
end

%% 7. Compute Clinical Metrics
G_OL_mgdl = (x_OL(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_NP_mgdl = (x_NP(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_MP_mgdl = (x_MP(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;

calc_metrics = @(G, u) struct( ...
    'PeakG', max(G), ...
    't_Peak', time_vec(find(G == max(G), 1)), ...
    'Excursion', max(G) - 100, ...
    'MinG', min(G), ...
    'TIR', 100 * mean(G >= 70 & G <= 180), ...
    'TAR', 100 * mean(G > 180), ...
    'TBR', 100 * mean(G < 70), ...
    'MaxU', max(u), ...
    'TotalU', sum(u) * (Ts / 60) ...
);

m_OL = calc_metrics(G_OL_mgdl, u_OL);
m_NP = calc_metrics(G_NP_mgdl, u_NP);
m_MP = calc_metrics(G_MP_mgdl, u_MP);

fprintf('\n========================================================================\n');
fprintf('  CLINICAL COMPARISON SUMMARY: %.1fg CHO MEAL (150 U/h Bolus Speed)\n', meal_carbs_g);
fprintf('========================================================================\n');
fprintf('%-28s %-16s %-18s %-18s\n', 'Metric', 'Open-Loop (Bolus)', 'MPC No-Preview', 'MPC Meal-Preview');
fprintf('------------------------------------------------------------------------\n');
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Peak Glucose (mg/dL)', m_OL.PeakG, m_NP.PeakG, m_MP.PeakG);
fprintf('%-28s %-16.1f %-18.1f %-18.1f\n', 'Time to Peak (min)', m_OL.t_Peak, m_NP.t_Peak, m_MP.t_Peak);
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Excursion above 100 (mg/dL)', m_OL.Excursion, m_NP.Excursion, m_MP.Excursion);
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Minimum Glucose (mg/dL)', m_OL.MinG, m_NP.MinG, m_MP.MinG);
fprintf('%-28s %-16.1f %-18.1f %-18.1f\n', 'Time in Range 70-180 (%)', m_OL.TIR, m_NP.TIR, m_MP.TIR);
fprintf('%-28s %-16.1f %-18.1f %-18.1f\n', 'Time above 180 (%)', m_OL.TAR, m_NP.TAR, m_MP.TAR);
fprintf('%-28s %-16.1f %-18.1f %-18.1f\n', 'Time below 70 (%)', m_OL.TBR, m_NP.TBR, m_MP.TBR);
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Prandial Bolus (U)', bolus_U, bolus_U, bolus_U);
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Total Insulin Delivered (U)', m_OL.TotalU, m_NP.TotalU, m_MP.TotalU);
fprintf('========================================================================\n');

%% 8. Directories Setup
fig_dir  = fullfile(project_root, 'results', 'figures', 'custom_input');
data_dir = fullfile(project_root, 'results', 'data', 'custom_input');
if ~exist(fig_dir, 'dir'),  mkdir(fig_dir);  end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

save_clean = @(fh, fp) clean_save_helper(fh, fp);

%% 9. Generate & Save Figures
% --- Figure 1: Combined 3-Panel Comparison ---
fig_comp = figure('Color', 'w', 'Position', [100, 80, 1050, 880], 'Visible', 'on');

% Subplot 1: Blood Glucose Response
subplot(3, 1, 1);
plot(time_vec, G_OL_mgdl, 'r--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL, TIR: %.1f%%)', m_OL.PeakG, m_OL.TIR));
hold on;
plot(time_vec, G_NP_mgdl, 'b-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('MPC No Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_NP.PeakG, m_NP.TIR));
plot(time_vec, G_MP_mgdl, 'g-', 'LineWidth', 2.2, ...
    'DisplayName', sprintf('MPC Meal Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_MP.PeakG, m_MP.TIR));
yline(100, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Target (100 mg/dL)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180 mg/dL)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70 mg/dL)');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Target Range [70-180]');
ylabel('Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Hovorka Artificial Pancreas: Custom Meal (%.1fg CHO, 150 U/h Bolus Speed)', ...
    meal_carbs_g), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, t_end]);
ylim_max = max(220, ceil(max([m_OL.PeakG, m_NP.PeakG, m_MP.PeakG, 180]) / 20) * 20 + 20);
ylim([50, ylim_max]);

% Subplot 2: Insulin Delivery Profile
subplot(3, 1, 2);
u_time = time_vec(1:end-1);
stairs(u_time, u_OL, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Open-Loop (Bolus + Basal)');
hold on;
stairs(u_time, u_NP, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No Preview (Max: %.1f U/h)', m_NP.MaxU));
stairs(u_time, u_MP, 'g-',  'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal Preview (Max: %.1f U/h)', m_MP.MaxU));
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'MPC Basal Ceiling (4.0 U/h)');
ylabel('Insulin Rate (U/h)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Insulin Delivery: %.2f U Bolus at 150 U/h (%.1f s) + Safe MPC Basal [0.0, 4.0] U/h', ...
    bolus_U, t_bolus_min * 60), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on; xlim([0, t_end]);
max_u_plot = max(6.0, ceil(max([m_OL.MaxU, m_NP.MaxU, m_MP.MaxU]) / 5) * 5 + 5);
ylim([0, max_u_plot]);

% Subplot 3: Meal Intake Profile
subplot(3, 1, 3);
stem(time_vec, d_meal_gmin, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'Marker', 'none', ...
    'DisplayName', sprintf('Custom Meal (%.1fg CHO over %.0f min)', meal_carbs_g, meal_dur));
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Meal Rate (g/min)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Prandial Glucose Disturbance (Total CHO: %.1fg at %.2f g/min)', ...
    meal_carbs_g, meal_rate_gmin), 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on; xlim([0, t_end]);

comp_fig_path = fullfile(fig_dir, 'custom_meal_comparison.png');
save_clean(fig_comp, comp_fig_path);

% --- Figure 2: Individual Case 1 Open-Loop ---
fig_c1 = figure('Color', 'w', 'Position', [120, 100, 900, 600], 'Visible', 'off');
subplot(2, 1, 1);
plot(time_vec, G_OL_mgdl, 'r-', 'LineWidth', 2.0); hold on;
yline(100, 'k:'); yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':'); yline(70, 'Color', [0.49 0.18 0.56], 'LineStyle', ':');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
title(sprintf('Case 1: Open-Loop Baseline (Peak: %.1f mg/dL, TIR: %.1f%%)', m_OL.PeakG, m_OL.TIR), 'FontWeight', 'bold');
grid on; xlim([0 t_end]); ylim([50 ylim_max]);
subplot(2, 1, 2);
stairs(u_time, u_OL, 'r-', 'LineWidth', 2.0); hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)');
xlabel('Time (min)', 'FontWeight', 'bold'); ylabel('Insulin (U/h)', 'FontWeight', 'bold');
title(sprintf('Open-Loop Delivery (Bolus: %.2f U at 150 U/h, Total: %.2f U)', bolus_U, m_OL.TotalU), 'FontWeight', 'bold');
grid on; xlim([0 t_end]); ylim([0 max_u_plot]);
save_clean(fig_c1, fullfile(fig_dir, 'custom_meal_case1_open_loop.png'));
close(fig_c1);

% --- Figure 3: Individual Case 2 MPC No-Preview ---
fig_c2 = figure('Color', 'w', 'Position', [140, 120, 900, 600], 'Visible', 'off');
subplot(2, 1, 1);
plot(time_vec, G_NP_mgdl, 'b-', 'LineWidth', 2.0); hold on;
yline(100, 'k:'); yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':'); yline(70, 'Color', [0.49 0.18 0.56], 'LineStyle', ':');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
title(sprintf('Case 2: MPC Without Meal Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_NP.PeakG, m_NP.TIR), 'FontWeight', 'bold');
grid on; xlim([0 t_end]); ylim([50 ylim_max]);
subplot(2, 1, 2);
stairs(u_time, u_NP, 'b-', 'LineWidth', 2.0); hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)');
xlabel('Time (min)', 'FontWeight', 'bold'); ylabel('Insulin (U/h)', 'FontWeight', 'bold');
title(sprintf('MPC No-Preview Delivery (Max: %.2f U/h, Total: %.2f U)', m_NP.MaxU, m_NP.TotalU), 'FontWeight', 'bold');
grid on; xlim([0 t_end]); ylim([0 max_u_plot]);
save_clean(fig_c2, fullfile(fig_dir, 'custom_meal_case2_mpc_no_preview.png'));
close(fig_c2);

% --- Figure 4: Individual Case 3 MPC Meal-Preview ---
fig_c3 = figure('Color', 'w', 'Position', [160, 140, 900, 600], 'Visible', 'off');
subplot(2, 1, 1);
plot(time_vec, G_MP_mgdl, 'g-', 'LineWidth', 2.0); hold on;
yline(100, 'k:'); yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':'); yline(70, 'Color', [0.49 0.18 0.56], 'LineStyle', ':');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
title(sprintf('Case 3: MPC With Meal Preview (Peak: %.1f mg/dL, TIR: %.1f%%)', m_MP.PeakG, m_MP.TIR), 'FontWeight', 'bold');
grid on; xlim([0 t_end]); ylim([50 ylim_max]);
subplot(2, 1, 2);
stairs(u_time, u_MP, 'g-', 'LineWidth', 2.0); hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)');
xlabel('Time (min)', 'FontWeight', 'bold'); ylabel('Insulin (U/h)', 'FontWeight', 'bold');
title(sprintf('MPC Meal-Preview Delivery (Max: %.2f U/h, Total: %.2f U)', m_MP.MaxU, m_MP.TotalU), 'FontWeight', 'bold');
grid on; xlim([0 t_end]); ylim([0 max_u_plot]);
save_clean(fig_c3, fullfile(fig_dir, 'custom_meal_case3_mpc_meal_preview.png'));
close(fig_c3);

% Save results dataset
mat_path = fullfile(data_dir, 'custom_meal_results.mat');
save(mat_path, 'meal_carbs_g', 'meal_dur', 'meal_start_time', 'bolus_U', ...
    'time_vec', 'G_OL_mgdl', 'G_NP_mgdl', 'G_MP_mgdl', ...
    'u_OL', 'u_NP', 'u_MP', 'd_meal_gmin', 'm_OL', 'm_NP', 'm_MP', 'p');
fprintf('Results MAT file saved to: %s\n', mat_path);
fprintf('========================================================================\n');
end

%% Clean save helper (overwrites cleanly)
function clean_save_helper(fig_handle, file_path)
    if exist(file_path, 'file')
        try, delete(file_path); catch, end
    end
    try
        exportgraphics(fig_handle, file_path, 'Resolution', 300);
    catch
        saveas(fig_handle, file_path);
    end
    fprintf('Saved: %s\n', file_path);
end
