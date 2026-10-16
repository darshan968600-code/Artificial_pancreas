% RUN_AND_SAVE_ALL_CASES Simulates the 3-Way Hovorka Comparison Simulink Model
% and exports labeled figures for each case + the combined comparison into results/figures/.

clc; clear; close all;
fprintf('=========================================================\n');
fprintf('  RUNNING 3-WAY COMPARISON & EXPORTING ALL SCOPE PLOTS\n');
fprintf('=========================================================\n');

% Robust project root resolution
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts')
    project_root = fileparts(this_dir);
else
    project_root = this_dir;
end
addpath(fullfile(project_root, 'scripts'));
addpath(fullfile(project_root, 'models'));

% Ensure directories exist
fig_dir  = fullfile(project_root, 'results', 'figures');
data_dir = fullfile(project_root, 'results', 'data');
if ~exist(fig_dir, 'dir'),  mkdir(fig_dir);  end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

% Load parameters and controller into base workspace
p = hovorka_params(70, 100);
[mpc_obj, plant_ss_d] = design_mpc_controller();
assignin('base', 'p', p);
assignin('base', 'mpc_obj', mpc_obj);
assignin('base', 'plant_ss_d', plant_ss_d);

% Simulate the 3-way Simulink model
model_name = 'hovorka_3way_comparison';
model_path = fullfile(project_root, 'models', [model_name, '.slx']);
fprintf('Simulating %s (720 min / 12 hours)...\n', model_path);
load_system(model_path);
sim_out = sim(model_name);

% Extract signals with their respective time vectors
G_OL_val = squeeze(sim_out.G_OL_sim.Data);
t_G_OL   = squeeze(sim_out.G_OL_sim.Time);
I_OL_val = squeeze(sim_out.I_OL_sim.Data);
t_I_OL   = squeeze(sim_out.I_OL_sim.Time);

G_NP_val = squeeze(sim_out.G_NP_sim.Data);
t_G_NP   = squeeze(sim_out.G_NP_sim.Time);
I_NP_val = squeeze(sim_out.I_NP_sim.Data);
t_I_NP   = squeeze(sim_out.I_NP_sim.Time);

G_MP_val = squeeze(sim_out.G_MP_sim.Data);
t_G_MP   = squeeze(sim_out.G_MP_sim.Time);
I_MP_val = squeeze(sim_out.I_MP_sim.Data);
t_I_MP   = squeeze(sim_out.I_MP_sim.Time);

%% Helper function to save figure cleanly
function save_clean_plot(fig_handle, file_path)
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

%% =========================================================
% FIGURE 1: Scope 1 - Case 1 Open-Loop Baseline
%% =========================================================
f1 = figure('Color', 'w', 'Position', [100, 100, 850, 600]);
subplot(2, 1, 1);
plot(t_G_OL, G_OL_val, 'r-', 'LineWidth', 2.0);
hold on;
yline(100, 'k--', 'Target (100 mg/dL)', 'LineWidth', 1.2);
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70)');
patch([0 720 720 0], [70 70 180 180], [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title(sprintf('Scope 1: Case 1 - Open-Loop Baseline (Peak: %.1f mg/dL)', max(G_OL_val)), 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Glucose (mg/dL)', 'FontWeight', 'bold');
xlim([0, 720]); ylim([50, 300]); grid on;

subplot(2, 1, 2);
plot(t_I_OL, I_OL_val, 'r-', 'LineWidth', 2.0);
hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.2);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Max Limit (4.0 U/h)');
title('Insulin Infusion Rate (Fixed Basal = 0.392 U/h)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (U/h)', 'FontWeight', 'bold');
xlim([0, 720]); ylim([0, 4.5]); grid on;

save_clean_plot(f1, fullfile(fig_dir, 'case1_open_loop.png'));

%% =========================================================
% FIGURE 2: Scope 2 - Case 2 MPC Without Meal Preview
%% =========================================================
f2 = figure('Color', 'w', 'Position', [150, 120, 850, 600]);
subplot(2, 1, 1);
plot(t_G_NP, G_NP_val, 'b-', 'LineWidth', 2.0);
hold on;
yline(100, 'k--', 'Target (100 mg/dL)', 'LineWidth', 1.2);
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70)');
patch([0 720 720 0], [70 70 180 180], [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title(sprintf('Scope 2: Case 2 - MPC Without Preview (Peak: %.1f mg/dL, Min: %.1f mg/dL)', ...
    max(G_NP_val), min(G_NP_val)), 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Glucose (mg/dL)', 'FontWeight', 'bold');
xlim([0, 720]); ylim([50, 300]); grid on;

subplot(2, 1, 2);
stairs(t_I_NP, I_NP_val, 'b-', 'LineWidth', 2.0);
hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Max Limit (4.0 U/h)');
title(sprintf('Insulin Infusion Rate (Max: %.2f U/h)', max(I_NP_val)), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (U/h)', 'FontWeight', 'bold');
xlim([0, 720]); ylim([0, 4.5]); grid on;

save_clean_plot(f2, fullfile(fig_dir, 'case2_mpc_no_preview.png'));

%% =========================================================
% FIGURE 3: Scope 3 - Case 3 MPC With Meal Preview
%% =========================================================
f3 = figure('Color', 'w', 'Position', [200, 140, 850, 600]);
subplot(2, 1, 1);
plot(t_G_MP, G_MP_val, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.0);
hold on;
yline(100, 'k--', 'Target (100 mg/dL)', 'LineWidth', 1.2);
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70)');
patch([0 720 720 0], [70 70 180 180], [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
title(sprintf('Scope 3: Case 3 - MPC With Meal Preview (Peak: %.1f mg/dL, Min: %.1f mg/dL)', ...
    max(G_MP_val), min(G_MP_val)), 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Glucose (mg/dL)', 'FontWeight', 'bold');
xlim([0, 720]); ylim([50, 300]); grid on;

subplot(2, 1, 2);
stairs(t_I_MP, I_MP_val, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.0);
hold on;
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Max Limit (4.0 U/h)');
title(sprintf('Insulin Infusion Rate (Max: %.2f U/h)', max(I_MP_val)), 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (U/h)', 'FontWeight', 'bold');
xlim([0, 720]); ylim([0, 4.5]); grid on;

save_clean_plot(f3, fullfile(fig_dir, 'case3_mpc_meal_preview.png'));

%% =========================================================
% FIGURE 4: Scope 4 - Combined 3-Way Comparison
%% =========================================================
f4 = figure('Color', 'w', 'Position', [250, 80, 1000, 750]);
subplot(2, 1, 1);
plot(t_G_OL, G_OL_val, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Case 1: Open-Loop (Peak: %.1f)', max(G_OL_val)));
hold on;
plot(t_G_NP, G_NP_val, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('Case 2: MPC No-Preview (Peak: %.1f)', max(G_NP_val)));
plot(t_G_MP, G_MP_val, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('Case 3: MPC Meal-Preview (Peak: %.1f)', max(G_MP_val)));
yline(100, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Target (100)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70)');
patch([0 720 720 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Target Range [70-180]');
title('Scope 4: Combined Artificial Pancreas Comparison (Hovorka Model, 50g CHO Meal)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
xlim([0, 720]); ylim([50, 300]); grid on;

subplot(2, 1, 2);
plot(t_I_OL, I_OL_val, 'r--', 'LineWidth', 1.6, 'DisplayName', 'Case 1: Open-Loop (Basal)');
hold on;
stairs(t_I_NP, I_NP_val, 'b-',  'LineWidth', 2.0, 'DisplayName', 'Case 2: MPC No-Preview');
stairs(t_I_MP, I_MP_val, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', 'Case 3: MPC Meal-Preview');
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Max Limit (4.0 U/h)');
title('Insulin Infusion Delivery Profiles (U/h)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Simulation Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin Infusion (U/h)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, 720]); ylim([0, 4.5]); grid on;

save_clean_plot(f4, fullfile(fig_dir, 'all_cases_comparison.png'));

% Save workspace data
mat_path = fullfile(data_dir, 'simulink_3way_data.mat');
save(mat_path, 'sim_out', 'G_OL_val', 't_G_OL', 'I_OL_val', 't_I_OL', ...
    'G_NP_val', 't_G_NP', 'I_NP_val', 't_I_NP', ...
    'G_MP_val', 't_G_MP', 'I_MP_val', 't_I_MP', 'p');
fprintf('Simulation data successfully saved to: %s\n', mat_path);
fprintf('=========================================================\n');
