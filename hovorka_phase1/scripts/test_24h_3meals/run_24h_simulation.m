% RUN_24H_SIMULATION 24-Hour Multi-Meal Artificial Pancreas Closed-Loop Study
% Evaluates Hovorka (2004) model over 1440 minutes with 3 meals of varying carb content:
%   - Meal 1 (Breakfast): 45g CHO at t = 180 min (07:00) over 15 min (3.0 g/min)
%   - Meal 2 (Lunch):     70g CHO at t = 480 min (12:00) over 20 min (3.5 g/min)
%   - Meal 3 (Dinner):    60g CHO at t = 840 min (18:00) over 20 min (3.0 g/min)
%   - Overnight fasting:  t = 860 to 1440 min
%
% Compares 3 branches:
%   1. Case 1: Open-Loop Baseline (Fixed Basal)
%   2. Case 2: MPC Without Meal Preview (Reactive)
%   3. Case 3: MPC With Meal Preview (Proactive)

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  24-HOUR MULTI-MEAL STUDY (HOVORKA AP MODEL)\n');
fprintf('  Breakfast (45g), Lunch (70g), Dinner (60g), Overnight Recovery\n');
fprintf('========================================================================\n');

% 1. Robust project paths
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir)); % up from scripts/test_24h_3meals to project_root
addpath(fullfile(project_root, 'scripts'));
addpath(fullfile(project_root, 'models'));

% Ensure output directories exist
fig_dir  = fullfile(project_root, 'results', 'figures', '24h_3meals');
data_dir = fullfile(project_root, 'results', 'data');
if ~exist(fig_dir, 'dir'),  mkdir(fig_dir);  end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

% 2. Load Hovorka Parameters & Controller
p = hovorka_params(70, 100);
[mpc_obj, plant_ss_d] = design_mpc_controller();

% 3. Simulation Timeline Setup
Ts = 5;                  % MPC sample time (minutes)
t_end = 1440;            % 24 hours (minutes)
time_vec = 0:Ts:t_end;   % 289 points
N_steps = length(time_vec) - 1; % 288 control intervals

% Construct 3-Meal Carbohydrate Disturbance Profile d_meal (g/min)
d_meal_gmin = zeros(size(time_vec));
for i = 1:length(time_vec)
    t = time_vec(i);
    % Meal 1 (Breakfast, 45g): t in [180, 195] min -> 3.0 g/min
    if t >= 180 && t < 195
        d_meal_gmin(i) = 3.0;
    % Meal 2 (Lunch, 70g): t in [480, 500] min -> 3.5 g/min
    elseif t >= 480 && t < 500
        d_meal_gmin(i) = 3.5;
    % Meal 3 (Dinner, 60g): t in [840, 860] min -> 3.0 g/min
    elseif t >= 840 && t < 860
        d_meal_gmin(i) = 3.0;
    end
end

target_mgdl = 100;
Np = mpc_obj.PredictionHorizon;
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

%% =====================================================================
% CASE 1: Open-Loop Baseline (Fixed Basal Infusion)
%% =====================================================================
fprintf('Simulating Case 1: Open-Loop Baseline (Fixed Basal = %.3f U/h)...\n', p.ss.u_basal_Uh);
x_OL = zeros(N_steps + 1, 11);
u_OL = p.ss.u_basal_Uh * ones(N_steps, 1);
x_OL(1, :) = p.ss.x0';

for k = 1:N_steps
    t_span_k = [time_vec(k), time_vec(k+1)];
    u_k_mumin = u_OL(k) * p.UH_TO_MUMIN;
    d_k_mmolmin = d_meal_gmin(k) * (1000 / p.MOL_WT_GLUCOSE);
    
    [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
        t_span_k, x_OL(k, :)', ode_opts);
    x_OL(k+1, :) = x_seg(end, :);
end

%% =====================================================================
% CASE 2: MPC Without Meal Preview (Reactive)
%% =====================================================================
fprintf('Simulating Case 2: MPC Without Meal Preview (Reactive)...\n');
x_NP = zeros(N_steps + 1, 11);
u_NP = zeros(N_steps, 1);
x_NP(1, :) = p.ss.x0';
xmpc_NP = mpcstate(mpc_obj);

for k = 1:N_steps
    t_span_k = [time_vec(k), time_vec(k+1)];
    y_meas_k = x_NP(k, 11) * p.MMOLL_TO_MGDL;
    d_curr = d_meal_gmin(k);
    
    xmpc_NP.Plant = x_NP(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_NP, y_meas_k, target_mgdl, d_curr);
    u_NP(k) = u_opt;
    
    u_k_mumin = u_opt * p.UH_TO_MUMIN;
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
        t_span_k, x_NP(k, :)', ode_opts);
    x_NP(k+1, :) = x_seg(end, :);
end

%% =====================================================================
% CASE 3: MPC With Meal Preview (Proactive Lookahead)
%% =====================================================================
fprintf('Simulating Case 3: MPC With Meal Preview (Proactive)...\n');
x_MP = zeros(N_steps + 1, 11);
u_MP = zeros(N_steps, 1);
x_MP(1, :) = p.ss.x0';
xmpc_MP = mpcstate(mpc_obj);

for k = 1:N_steps
    t_span_k = [time_vec(k), time_vec(k+1)];
    y_meas_k = x_MP(k, 11) * p.MMOLL_TO_MGDL;
    
    % Lookahead preview vector over prediction horizon (Np+1 samples)
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        idx_future = min(k + h, N_steps + 1);
        v_preview(h + 1) = d_meal_gmin(idx_future);
    end
    
    xmpc_MP.Plant = x_MP(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_MP, y_meas_k, target_mgdl, v_preview);
    u_MP(k) = u_opt;
    
    u_k_mumin = u_opt * p.UH_TO_MUMIN;
    d_k_mmolmin = d_meal_gmin(k) * (1000 / p.MOL_WT_GLUCOSE);
    
    [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
        t_span_k, x_MP(k, :)', ode_opts);
    x_MP(k+1, :) = x_seg(end, :);
end

%% =====================================================================
% 4. Compute Comprehensive 24-Hour Clinical Metrics
%% =====================================================================
G_OL_mgdl = (x_OL(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_NP_mgdl = (x_NP(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_MP_mgdl = (x_MP(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;

calc_metrics_24h = @(G, u) struct( ...
    'MeanG', mean(G), ...
    'PeakG', max(G), ...
    'MinG', min(G), ...
    'TIR', 100 * mean(G >= 70 & G <= 180), ...
    'TAR', 100 * mean(G > 180), ...
    'TBR', 100 * mean(G < 70), ...
    'MaxU', max(u), ...
    'TotalU', sum(u) * (Ts / 60));

m_OL = calc_metrics_24h(G_OL_mgdl, u_OL);
m_NP = calc_metrics_24h(G_NP_mgdl, u_NP);
m_MP = calc_metrics_24h(G_MP_mgdl, u_MP);

fprintf('\n====================================================================================\n');
fprintf('  24-HOUR CLINICAL METRICS SUMMARY: 3 MEALS (45g Breakfast, 70g Lunch, 60g Dinner)\n');
fprintf('====================================================================================\n');
fprintf('%-28s %-16s %-18s %-18s\n', 'Metric', 'Open-Loop', 'MPC No-Preview', 'MPC Meal-Preview');
fprintf('------------------------------------------------------------------------------------\n');
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Mean Glucose (mg/dL)', m_OL.MeanG, m_NP.MeanG, m_MP.MeanG);
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Peak Glucose (mg/dL)', m_OL.PeakG, m_NP.PeakG, m_MP.PeakG);
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Minimum Glucose (mg/dL)', m_OL.MinG, m_NP.MinG, m_MP.MinG);
fprintf('%-28s %-16.1f %-18.1f %-18.1f\n', 'Time in Range 70-180 (%)', m_OL.TIR, m_NP.TIR, m_MP.TIR);
fprintf('%-28s %-16.1f %-18.1f %-18.1f\n', 'Time above 180 (%)', m_OL.TAR, m_NP.TAR, m_MP.TAR);
fprintf('%-28s %-16.1f %-18.1f %-18.1f\n', 'Time below 70 (%)', m_OL.TBR, m_NP.TBR, m_MP.TBR);
fprintf('%-28s %-16.3f %-18.3f %-18.3f\n', 'Max Insulin Infusion (U/h)', m_OL.MaxU, m_NP.MaxU, m_MP.MaxU);
fprintf('%-28s %-16.2f %-18.2f %-18.2f\n', 'Total Insulin Delivered (U)', m_OL.TotalU, m_NP.TotalU, m_MP.TotalU);
fprintf('====================================================================================\n');

%% =====================================================================
% Helper: Clean Figure Saver
%% =====================================================================
function save_clean_figure(fig_handle, file_path)
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

%% =====================================================================
% 5. FIGURE 1: Comprehensive 3-Panel 24-Hour Comparison
%% =====================================================================
f1 = figure('Color', 'w', 'Position', [80, 50, 1100, 900]);

% Subplot 1: Glucose Profiles
subplot(3, 1, 1);
plot(time_vec, G_OL_mgdl, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Case 1: Open-Loop (Peak: %.1f)', m_OL.PeakG));
hold on;
plot(time_vec, G_NP_mgdl, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('Case 2: MPC No-Preview (Peak: %.1f)', m_NP.PeakG));
plot(time_vec, G_MP_mgdl, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('Case 3: MPC Meal-Preview (Peak: %.1f)', m_MP.PeakG));
yline(100, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Target (100 mg/dL)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70)');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Target Range [70-180]');
title('24-Hour Closed-Loop Glucose Regulation: Hovorka AP Model (3 Meals: 45g, 70g, 60g)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
xlim([0, t_end]); ylim([50, 320]); grid on;

% Subplot 2: Insulin Infusion Profiles
subplot(3, 1, 2);
u_time = time_vec(1:end-1);
stairs(u_time, u_OL, 'r--', 'LineWidth', 1.6, 'DisplayName', sprintf('Case 1: Fixed Basal (%.3f U/h)', p.ss.u_basal_Uh));
hold on;
stairs(u_time, u_NP, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('Case 2: MPC No-Preview (Max: %.2f U/h)', m_NP.MaxU));
stairs(u_time, u_MP, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('Case 3: MPC Meal-Preview (Max: %.2f U/h)', m_MP.MaxU));
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Max Limit (4.0 U/h)');
title('Insulin Infusion Delivery Profiles (U/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Insulin Infusion (U/h)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, t_end]); ylim([0, 4.5]); grid on;

% Subplot 3: Meal Intakes
subplot(3, 1, 3);
stem(time_vec, d_meal_gmin, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'Marker', 'none', 'DisplayName', 'Meal Intake Rate (g/min)');
hold on;
text(187, 3.2, 'Breakfast (45g)', 'FontWeight', 'bold', 'Color', [0.7 0.2 0.1]);
text(490, 3.7, 'Heavy Lunch (70g)', 'FontWeight', 'bold', 'Color', [0.7 0.2 0.1]);
text(850, 3.2, 'Dinner (60g)', 'FontWeight', 'bold', 'Color', [0.7 0.2 0.1]);
text(1100, 1.0, 'Overnight Fasting Period', 'FontWeight', 'bold', 'Color', [0.2 0.4 0.7]);
title('Prandial Carbohydrate Intakes Across 24 Hours', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time of Day (minutes from midnight / 0 to 1440 min)', 'FontWeight', 'bold');
ylabel('Meal Rate (g/min)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, t_end]); ylim([0, 4.5]); grid on;

save_clean_figure(f1, fullfile(fig_dir, 'all_cases_comparison_24h.png'));

%% =====================================================================
% 6. INDIVIDUAL SCOPE FIGURES (1 for each case)
%% =====================================================================
% Case 1 Open Loop 24h
f_c1 = figure('Color', 'w', 'Position', [100, 100, 850, 600]);
subplot(2, 1, 1);
plot(time_vec, G_OL_mgdl, 'r-', 'LineWidth', 2.0);
hold on;
yline(100, 'k--', 'Target (100 mg/dL)', 'LineWidth', 1.2);
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2);
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2);
patch([0 t_end t_end 0], [70 70 180 180], [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title(sprintf('Case 1: Open-Loop Baseline 24h (Peak: %.1f mg/dL, Mean: %.1f mg/dL)', m_OL.PeakG, m_OL.MeanG), 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Glucose (mg/dL)', 'FontWeight', 'bold');
xlim([0, t_end]); ylim([50, 320]); grid on;

subplot(2, 1, 2);
plot(u_time, u_OL, 'r-', 'LineWidth', 2.0);
hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2);
title('Insulin Infusion Rate (Fixed Basal = 0.392 U/h)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (U/h)', 'FontWeight', 'bold');
xlim([0, t_end]); ylim([0, 4.5]); grid on;
save_clean_figure(f_c1, fullfile(fig_dir, 'case1_open_loop_24h.png'));

% Case 2 MPC No-Preview 24h
f_c2 = figure('Color', 'w', 'Position', [150, 120, 850, 600]);
subplot(2, 1, 1);
plot(time_vec, G_NP_mgdl, 'b-', 'LineWidth', 2.0);
hold on;
yline(100, 'k--', 'Target (100 mg/dL)', 'LineWidth', 1.2);
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2);
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2);
patch([0 t_end t_end 0], [70 70 180 180], [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title(sprintf('Case 2: MPC Without Preview 24h (Peak: %.1f mg/dL, Min: %.1f mg/dL, TIR: %.1f%%)', ...
    m_NP.PeakG, m_NP.MinG, m_NP.TIR), 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Glucose (mg/dL)', 'FontWeight', 'bold');
xlim([0, t_end]); ylim([50, 320]); grid on;

subplot(2, 1, 2);
stairs(u_time, u_NP, 'b-', 'LineWidth', 2.0);
hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2);
title(sprintf('Insulin Infusion Rate (Max: %.2f U/h)', m_NP.MaxU), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (U/h)', 'FontWeight', 'bold');
xlim([0, t_end]); ylim([0, 4.5]); grid on;
save_clean_figure(f_c2, fullfile(fig_dir, 'case2_mpc_no_preview_24h.png'));

% Case 3 MPC Meal-Preview 24h
f_c3 = figure('Color', 'w', 'Position', [200, 140, 850, 600]);
subplot(2, 1, 1);
plot(time_vec, G_MP_mgdl, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.0);
hold on;
yline(100, 'k--', 'Target (100 mg/dL)', 'LineWidth', 1.2);
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2);
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2);
patch([0 t_end t_end 0], [70 70 180 180], [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title(sprintf('Case 3: MPC With Meal Preview 24h (Peak: %.1f mg/dL, Min: %.1f mg/dL, TIR: %.1f%%)', ...
    m_MP.PeakG, m_MP.MinG, m_MP.TIR), 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Glucose (mg/dL)', 'FontWeight', 'bold');
xlim([0, t_end]); ylim([50, 320]); grid on;

subplot(2, 1, 2);
stairs(u_time, u_MP, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.0);
hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2);
title(sprintf('Insulin Infusion Rate (Max: %.2f U/h)', m_MP.MaxU), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (U/h)', 'FontWeight', 'bold');
xlim([0, t_end]); ylim([0, 4.5]); grid on;
save_clean_figure(f_c3, fullfile(fig_dir, 'case3_mpc_meal_preview_24h.png'));

% Save workspace data
mat_path = fullfile(data_dir, 'simulation_24h_3meals.mat');
save(mat_path, 'time_vec', 'u_time', 'd_meal_gmin', ...
    'G_OL_mgdl', 'u_OL', 'm_OL', ...
    'G_NP_mgdl', 'u_NP', 'm_NP', ...
    'G_MP_mgdl', 'u_MP', 'm_MP', 'p');
fprintf('\n24-Hour simulation data successfully saved to: %s\n', mat_path);
fprintf('========================================================================\n');
