% VALIDATE_OPEN_LOOP Simulates open-loop meal and bolus challenges
% on the Hovorka (2004) T1D model to validate physiological realism.

clc; clear; close all;
fprintf('=========================================================\n');
fprintf('  STEP 2: HOVORKA OPEN-LOOP MEAL & BOLUS VALIDATION\n');
fprintf('=========================================================\n');

% Robust project root and path resolution
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts')
    project_root = fileparts(this_dir);
else
    project_root = this_dir;
end
addpath(fullfile(project_root, 'scripts'));

% Load nominal patient parameters (70 kg adult, 100 mg/dL target)
p = hovorka_params(70, 100);

% Simulation time span: 0 to 720 min (12 hours)
t_start = 0;
t_end   = 720; % min
tspan   = [t_start, t_end];

%% Meal Profile: 50g Carbohydrate meal at t = 120 min over 15 min duration
meal_time = 120;       % min
meal_dur  = 15;        % min
meal_carbs_g = 50;     % grams
meal_carbs_mmol = meal_carbs_g * (1000 / p.MOL_WT_GLUCOSE); % ~277.5 mmol
meal_rate_mmol_min = meal_carbs_mmol / meal_dur;             % mmol/min

D_meal_func = @(t) (t >= meal_time & t <= (meal_time + meal_dur)) * meal_rate_mmol_min;

%% Case A: Basal Only (Unannounced meal, no bolus)
u_basal_mumin = p.ss.u_basal_mumin;
u_caseA_func  = @(t) u_basal_mumin;

fprintf('Running Case A: 50g Meal with Basal Insulin Only (%.2f U/h)...\n', p.ss.u_basal_Uh);
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1.0);
[t_A, x_A] = ode15s(@(t, x) hovorka_ode(t, x, u_caseA_func(t), D_meal_func(t), p), tspan, p.ss.x0, opts);

%% Case B: Basal + Standard Manual Bolus (Standard clinical practice)
% Standard adult Insulin-to-Carb Ratio (ICR): ~12 g/U
% Bolus for 50g meal: 50 / 12 = ~4.17 U = 4170 mU
% Delivered as a subcutaneous square infusion over 10 min starting at t = 120 min
bolus_U = meal_carbs_g / 12.0; % ~4.17 U
bolus_dur = 10; % min
bolus_rate_mumin = (bolus_U * p.U_TO_MU) / bolus_dur; % mU/min

u_caseB_func = @(t) u_basal_mumin + ...
    (t >= meal_time & t <= (meal_time + bolus_dur)) * bolus_rate_mumin;

fprintf('Running Case B: 50g Meal + Manual Bolus (%.2f U over %d min)...\n', bolus_U, bolus_dur);
[t_B, x_B] = ode15s(@(t, x) hovorka_ode(t, x, u_caseB_func(t), D_meal_func(t), p), tspan, p.ss.x0, opts);

%% Compute Signals for Analysis
% Case A
G_A_mgdl      = (x_A(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_meas_A_mgdl = x_A(:, 11) .* p.MMOLL_TO_MGDL;
I_A           = x_A(:, 3);
UG_A          = x_A(:, 10) ./ p.tau_D;
u_A_Uh        = arrayfun(@(t) u_caseA_func(t) * p.MUMIN_TO_UH, t_A);

% Case B
G_B_mgdl      = (x_B(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_meas_B_mgdl = x_B(:, 11) .* p.MMOLL_TO_MGDL;
I_B           = x_B(:, 3);
UG_B          = x_B(:, 10) ./ p.tau_D;
u_B_Uh        = arrayfun(@(t) u_caseB_func(t) * p.MUMIN_TO_UH, t_B);

%% Calculate Key Metrics
[G_peak_A, idx_peak_A] = max(G_A_mgdl);
t_peak_A = t_A(idx_peak_A);

[G_peak_B, idx_peak_B] = max(G_B_mgdl);
t_peak_B = t_B(idx_peak_B);

fprintf('---------------------------------------------------------\n');
fprintf('  RESULTS SUMMARY: 50g MEAL OPEN-LOOP CHALLENGE\n');
fprintf('---------------------------------------------------------\n');
fprintf('Case A (Basal Only):\n');
fprintf('  Peak Glucose:        %.2f mg/dL\n', G_peak_A);
fprintf('  Excursion above 100: +%.2f mg/dL\n', G_peak_A - 100);
fprintf('  Time to Peak:        %.1f min (%.1f min post-meal)\n', t_peak_A, t_peak_A - meal_time);
fprintf('  Final Glucose (12h): %.2f mg/dL\n', G_A_mgdl(end));

fprintf('\nCase B (Basal + Manual Bolus):\n');
fprintf('  Peak Glucose:        %.2f mg/dL\n', G_peak_B);
fprintf('  Excursion above 100: +%.2f mg/dL\n', G_peak_B - 100);
fprintf('  Time to Peak:        %.1f min (%.1f min post-meal)\n', t_peak_B, t_peak_B - meal_time);
fprintf('  Final Glucose (12h): %.2f mg/dL\n', G_B_mgdl(end));
fprintf('---------------------------------------------------------\n');

%% Generate Publication-Quality Plots
figure('Color', 'w', 'Position', [100, 100, 950, 800]);

% Subplot 1: Plasma Glucose Concentration
subplot(4, 1, 1);
plot(t_A, G_A_mgdl, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Case A: Basal Only');
hold on;
plot(t_B, G_B_mgdl, 'b-', 'LineWidth', 2.0, 'DisplayName', sprintf('Case B: Basal + %.1f U Bolus', bolus_U));
yline(100, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Target (100 mg/dL)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia Threshold (180 mg/dL)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia Threshold (70 mg/dL)');
patch([0 t_end t_end 0], [70 70 180 180], [0.9 0.95 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Target Range (70-180 mg/dL)');
ylabel('Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title('Hovorka Model Open-Loop Validation: 50g Carbohydrate Meal Challenge', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, t_end]); ylim([60, max(G_peak_A + 20, 250)]);

% Subplot 2: Insulin Infusion Rate
subplot(4, 1, 2);
plot(t_A, u_A_Uh, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Case A (Basal 0.39 U/h)');
hold on;
plot(t_B, u_B_Uh, 'b-', 'LineWidth', 1.8, 'DisplayName', sprintf('Case B (Basal + Bolus)'));
ylabel('Insulin Infusion (U/h)', 'FontSize', 11, 'FontWeight', 'bold');
grid on; xlim([0, t_end]);
legend('Location', 'northeast');

% Subplot 3: Plasma Insulin Concentration
subplot(4, 1, 3);
plot(t_A, I_A, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Case A (Plasma Insulin)');
hold on;
plot(t_B, I_B, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Case B (Plasma Insulin)');
ylabel('Plasma Insulin (mU/L)', 'FontSize', 11, 'FontWeight', 'bold');
grid on; xlim([0, t_end]);
legend('Location', 'northeast');

% Subplot 4: Gut Glucose Appearance Rate (Meal)
subplot(4, 1, 4);
plot(t_A, UG_A, 'm-', 'LineWidth', 1.8, 'DisplayName', 'Meal Glucose Appearance U_G');
xlabel('Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('U_G (mmol/min)', 'FontSize', 11, 'FontWeight', 'bold');
grid on; xlim([0, t_end]);
legend('Location', 'northeast');

% Ensure target directories exist
fig_dir  = fullfile(project_root, 'results', 'figures');
data_dir = fullfile(project_root, 'results', 'data');
if ~exist(fig_dir, 'dir'),  mkdir(fig_dir);  end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

fig_path = fullfile(fig_dir, 'open_loop_validation.png');
mat_path = fullfile(data_dir, 'open_loop_validation.mat');

% Cleanly overwrite previous image if present
if exist(fig_path, 'file')
    try, delete(fig_path); catch, end
end

try
    exportgraphics(gcf, fig_path, 'Resolution', 300);
catch
    saveas(gcf, fig_path);
end
save(mat_path, 't_A', 'x_A', 't_B', 'x_B', 'p', 'G_peak_A', 't_peak_A', 'G_peak_B', 't_peak_B');
fprintf('Figure successfully saved to: %s\n', fig_path);
fprintf('Data successfully saved to:   %s\n', mat_path);
fprintf('=========================================================\n');
