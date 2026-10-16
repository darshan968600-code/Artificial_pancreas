%% ANALYZE_BERGMAN_RESULTS Analyzes simulation results and exports publication figures
% Matches and recreates the outputs of AP_results.m with enhanced plotting.

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  BERGMAN ARTIFICIAL PANCREAS: RESULTS ANALYSIS & EXPORT\n');
fprintf('========================================================================\n');

% 1. Dynamic Paths
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(this_dir);
data_dir  = fullfile(project_root, 'results', 'data');
fig_dir   = fullfile(project_root, 'results', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

mat_file = fullfile(data_dir, 'AP_simulation_data.mat');
if ~exist(mat_file, 'file')
    error('Simulation data file not found: %s\nPlease run run_all_bergman_simulations.m first.', mat_file);
end

load(mat_file, 'results');
fprintf('Simulation data loaded successfully from: %s\n\n', mat_file);

%% Helper: Clean Plot Exporter
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

%% =====================================================================
% 1. NORMAL MEAL METRICS (300 min)
%% =====================================================================
G_ol_norm = results.normal.openLoop.glucose.Data;
t_ol_norm = results.normal.openLoop.glucose.Time / 60; % convert seconds to minutes
I_ol_norm = results.normal.openLoop.insulin.Data;

G_np_norm = results.normal.mpcNoPreview.glucose.Data;
t_np_norm = results.normal.mpcNoPreview.glucose.Time / 60;
I_np_norm = results.normal.mpcNoPreview.insulin.Data;

G_mp_norm = results.normal.mpcPreview.glucose.Data;
t_mp_norm = results.normal.mpcPreview.glucose.Time / 60;
I_mp_norm = results.normal.mpcPreview.insulin.Data;

% Peak glucose
G_peak_ol = max(G_ol_norm);
G_peak_np = max(G_np_norm);
G_peak_mp = max(G_mp_norm);

% Min glucose
G_min_ol = min(G_ol_norm);
G_min_np = min(G_np_norm);
G_min_mp = min(G_mp_norm);

% Time to peak (minutes)
[~, idx] = max(G_ol_norm); t_peak_ol = t_ol_norm(idx);
[~, idx] = max(G_np_norm); t_peak_np = t_np_norm(idx);
[~, idx] = max(G_mp_norm); t_peak_mp = t_mp_norm(idx);

% Max insulin
I_max_ol = max(I_ol_norm);
I_max_np = max(I_np_norm);
I_max_mp = max(I_mp_norm);

% Peak reduction %
red_np = (G_peak_ol - G_peak_np) / G_peak_ol * 100;
red_mp = (G_peak_ol - G_peak_mp) / G_peak_ol * 100;

%% =====================================================================
% 2. LARGE MEAL METRICS (600 min)
%% =====================================================================
G_ol_large = results.large.openLoop.glucose.Data;
t_ol_large = results.large.openLoop.glucose.Time / 60;
I_ol_large = results.large.openLoop.insulin.Data;

G_np_large = results.large.mpcNoPreview.glucose.Data;
t_np_large = results.large.mpcNoPreview.glucose.Time / 60;
I_np_large = results.large.mpcNoPreview.insulin.Data;

G_mp_large = results.large.mpcPreview.glucose.Data;
t_mp_large = results.large.mpcPreview.glucose.Time / 60;
I_mp_large = results.large.mpcPreview.insulin.Data;

G_peak_ol_L = max(G_ol_large);
G_peak_np_L = max(G_np_large);
G_peak_mp_L = max(G_mp_large);

G_min_ol_L = min(G_ol_large);
G_min_np_L = min(G_np_large);
G_min_mp_L = min(G_mp_large);

[~, idx] = max(G_ol_large); t_peak_ol_L = t_ol_large(idx);
[~, idx] = max(G_np_large); t_peak_np_L = t_np_large(idx);
[~, idx] = max(G_mp_large); t_peak_mp_L = t_mp_large(idx);

I_max_ol_L = max(I_ol_large);
I_max_np_L = max(I_np_large);
I_max_mp_L = max(I_mp_large);

red_np_L = (G_peak_ol_L - G_peak_np_L) / G_peak_ol_L * 100;
red_mp_L = (G_peak_ol_L - G_peak_mp_L) / G_peak_ol_L * 100;

%% =====================================================================
% 3. EXTENDED 600-MIN METRICS
%% =====================================================================
G_np_ext = results.extended.mpcNoPreview.glucose.Data;
t_np_ext = results.extended.mpcNoPreview.glucose.Time / 60;
I_np_ext = results.extended.mpcNoPreview.insulin.Data;

G_mp_ext = results.extended.mpcPreview.glucose.Data;
t_mp_ext = results.extended.mpcPreview.glucose.Time / 60;
I_mp_ext = results.extended.mpcPreview.insulin.Data;

G_min_np_ext = min(G_np_ext);
G_min_mp_ext = min(G_mp_ext);
[~, idx] = min(G_np_ext); t_min_np_ext = t_np_ext(idx);
[~, idx] = min(G_mp_ext); t_min_mp_ext = t_mp_ext(idx);
G_fin_np_ext = G_np_ext(end);
G_fin_mp_ext = G_mp_ext(end);

%% =====================================================================
% DISPLAY PERFORMANCE SUMMARY
%% =====================================================================
fprintf('===== NORMAL MEAL (D_peak = 0.6 mg/dL/min) =====\n');
fprintf('Open-loop peak glucose       = %8.3f mg/dL (at %.1f min)\n', G_peak_ol, t_peak_ol);
fprintf('No-preview MPC peak glucose  = %8.3f mg/dL (at %.1f min)\n', G_peak_np, t_peak_np);
fprintf('Preview MPC peak glucose     = %8.3f mg/dL (at %.1f min)\n', G_peak_mp, t_peak_mp);
fprintf('Open-loop maximum insulin    = %8.3f uU/mL\n', I_max_ol);
fprintf('No-preview MPC max insulin   = %8.3f uU/mL\n', I_max_np);
fprintf('Preview MPC max insulin      = %8.3f uU/mL\n', I_max_mp);
fprintf('No-preview peak reduction    = %8.2f %%\n', red_np);
fprintf('Preview peak reduction       = %8.2f %%\n\n', red_mp);

fprintf('===== LARGE MEAL (D_peak = 1.2 mg/dL/min) =====\n');
fprintf('Open-loop peak glucose       = %8.3f mg/dL (at %.1f min)\n', G_peak_ol_L, t_peak_ol_L);
fprintf('No-preview MPC peak glucose  = %8.3f mg/dL (at %.1f min)\n', G_peak_np_L, t_peak_np_L);
fprintf('Preview MPC peak glucose     = %8.3f mg/dL (at %.1f min)\n', G_peak_mp_L, t_peak_mp_L);
fprintf('Open-loop maximum insulin    = %8.3f uU/mL\n', I_max_ol_L);
fprintf('No-preview MPC max insulin   = %8.3f uU/mL\n', I_max_np_L);
fprintf('Preview MPC max insulin      = %8.3f uU/mL\n', I_max_mp_L);
fprintf('No-preview peak reduction    = %8.2f %%\n', red_np_L);
fprintf('Preview peak reduction       = %8.2f %%\n\n', red_mp_L);

fprintf('===== EXTENDED 600-MINUTE RESULTS =====\n');
fprintf('No-preview minimum glucose   = %8.3f mg/dL (at %.1f min), final = %8.3f mg/dL\n', G_min_np_ext, t_min_np_ext, G_fin_np_ext);
fprintf('Preview minimum glucose      = %8.3f mg/dL (at %.1f min), final = %8.3f mg/dL\n\n', G_min_mp_ext, t_min_mp_ext, G_fin_mp_ext);

% Summary Table
experiments = { ...
    'Normal meal - Open-loop'; ...
    'Normal meal - MPC no preview'; ...
    'Normal meal - MPC with preview'; ...
    'Large meal - Open-loop'; ...
    'Large meal - MPC no preview'; ...
    'Large meal - MPC with preview'};

peakG = [G_peak_ol; G_peak_np; G_peak_mp; G_peak_ol_L; G_peak_np_L; G_peak_mp_L];
peakT = [t_peak_ol; t_peak_np; t_peak_mp; t_peak_ol_L; t_peak_np_L; t_peak_mp_L];
minG  = [G_min_ol; G_min_np; G_min_mp; G_min_ol_L; G_min_np_L; G_min_mp_L];
maxI  = [I_max_ol; I_max_np; I_max_mp; I_max_ol_L; I_max_np_L; I_max_mp_L];
redP  = [0.0; red_np; red_mp; 0.0; red_np_L; red_mp_L];

summary_table = table(experiments, peakG, peakT, minG, maxI, redP, ...
    'VariableNames', {'Experiment', 'PeakGlucose_mg_dL', 'PeakTime_min', 'MinimumGlucose_mg_dL', 'MaximumInsulin_uU_mL', 'PeakReduction_percent'});
disp(summary_table);

csv_file = fullfile(data_dir, 'bergman_final_summary.csv');
writetable(summary_table, csv_file);
fprintf('Summary table saved to: %s\n\n', csv_file);

%% =====================================================================
% FIGURE 1: Normal Meal Comparison (300 min)
%% =====================================================================
f1 = figure('Color', 'w', 'Position', [100, 100, 900, 650]);
subplot(2, 1, 1);
plot(t_ol_norm, G_ol_norm, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL)', G_peak_ol));
hold on;
plot(t_np_norm, G_np_norm, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No-Preview (Peak: %.1f mg/dL)', G_peak_np));
plot(t_mp_norm, G_mp_norm, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal-Preview (Peak: %.1f mg/dL)', G_peak_mp));
yline(90, 'k:', 'Basal (90 mg/dL)', 'LineWidth', 1.2);
title('Bergman Minimal Model: Normal Meal Response (D_{peak} = 0.6 mg/dL/min at t = 60 min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, 300]); grid on;

subplot(2, 1, 2);
plot(t_ol_norm, I_ol_norm, 'r--', 'LineWidth', 1.6, 'DisplayName', 'Open-Loop (Basal = 7.0 uU/mL)');
hold on;
plot(t_np_norm, I_np_norm, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No-Preview (Max: %.2f uU/mL)', I_max_np));
plot(t_mp_norm, I_mp_norm, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal-Preview (Max: %.2f uU/mL)', I_max_mp));
title('Plasma Insulin Concentration Profiles (uU/mL)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Simulation Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (uU/mL)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, 300]); grid on;

save_clean_plot(f1, fullfile(fig_dir, 'normal_meal_comparison.png'));

%% =====================================================================
% FIGURE 2: Large Meal Comparison (600 min)
%% =====================================================================
f2 = figure('Color', 'w', 'Position', [150, 120, 900, 650]);
subplot(2, 1, 1);
plot(t_ol_large, G_ol_large, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL)', G_peak_ol_L));
hold on;
plot(t_np_large, G_np_large, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No-Preview (Peak: %.1f mg/dL)', G_peak_np_L));
plot(t_mp_large, G_mp_large, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal-Preview (Peak: %.1f mg/dL)', G_peak_mp_L));
yline(90, 'k:', 'Basal (90 mg/dL)', 'LineWidth', 1.2);
title('Bergman Minimal Model: Large Meal Response (D_{peak} = 1.2 mg/dL/min at t = 60 min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, 600]); grid on;

subplot(2, 1, 2);
plot(t_ol_large, I_ol_large, 'r--', 'LineWidth', 1.6, 'DisplayName', 'Open-Loop (Basal = 7.0 uU/mL)');
hold on;
plot(t_np_large, I_np_large, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No-Preview (Max: %.2f uU/mL)', I_max_np_L));
plot(t_mp_large, I_mp_large, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal-Preview (Max: %.2f uU/mL)', I_max_mp_L));
title('Plasma Insulin Concentration Profiles (uU/mL)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Simulation Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (uU/mL)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, 600]); grid on;

save_clean_plot(f2, fullfile(fig_dir, 'large_meal_comparison.png'));

%% =====================================================================
% FIGURE 3: Extended Recovery (600 min)
%% =====================================================================
f3 = figure('Color', 'w', 'Position', [200, 140, 900, 650]);
subplot(2, 1, 1);
plot(t_np_ext, G_np_ext, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No-Preview (Min: %.1f mg/dL)', G_min_np_ext));
hold on;
plot(t_mp_ext, G_mp_ext, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal-Preview (Min: %.1f mg/dL)', G_min_mp_ext));
yline(90, 'k:', 'Basal (90 mg/dL)', 'LineWidth', 1.2);
title('Bergman Model: Extended 600-Minute Recovery Under Closed-Loop MPC', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Blood Glucose (mg/dL)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, 600]); grid on;

subplot(2, 1, 2);
plot(t_np_ext, I_np_ext, 'b-',  'LineWidth', 2.0, 'DisplayName', 'MPC No-Preview');
hold on;
plot(t_mp_ext, I_mp_ext, 'Color', [0.1 0.65 0.2], 'LineWidth', 2.2, 'DisplayName', 'MPC Meal-Preview');
title('Plasma Insulin Concentration Profiles (uU/mL)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Simulation Time (minutes)', 'FontWeight', 'bold');
ylabel('Insulin (uU/mL)', 'FontWeight', 'bold');
legend('Location', 'northeast');
xlim([0, 600]); grid on;

save_clean_plot(f3, fullfile(fig_dir, 'extended_recovery_comparison.png'));

fprintf('All analysis plots saved successfully to: %s\n', fig_dir);
fprintf('========================================================================\n');
