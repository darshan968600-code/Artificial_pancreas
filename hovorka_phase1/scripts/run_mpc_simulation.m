% RUN_MPC_SIMULATION Runs closed-loop comparison between:
%   1. Open-Loop Baseline (Basal only)
%   2. MPC Without Meal Preview (Reactive Closed-Loop)
%   3. MPC With Meal Preview (Proactive Feedforward Closed-Loop)
% on the full non-linear Hovorka (2004) T1D plant.

clc; clear; close all;
fprintf('=========================================================\n');
fprintf('  STEP 3 & 4: HOVORKA MPC CLOSED-LOOP COMPARISON STUDY\n');
fprintf('=========================================================\n');

% Robust project root and path resolution
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts')
    project_root = fileparts(this_dir);
else
    project_root = this_dir;
end
addpath(fullfile(project_root, 'scripts'));

% 1. Design and obtain MPC controller and parameters
[mpc_obj, plant_ss_d, p] = design_mpc_controller();

% 2. Simulation Timeline
Ts = 5;               % Control interval: 5 min
t_end = 720;          % Total duration: 720 min (12 hours)
N_steps = t_end / Ts; % 144 steps
time_vec = (0:N_steps)' * Ts;

% 3. Meal Disturbance Definition: 50g CHO at t = 120 min over 15 min
meal_start_time = 120; % min
meal_dur        = 15;  % min
meal_carbs_g    = 50;  % grams
meal_rate_gmin  = meal_carbs_g / meal_dur; % 3.333 g/min

% Discrete meal profile for every step
d_meal_gmin = zeros(N_steps + 1, 1);
for k = 1:(N_steps + 1)
    t_curr = time_vec(k);
    if t_curr >= meal_start_time && t_curr < (meal_start_time + meal_dur)
        d_meal_gmin(k) = meal_rate_gmin;
    end
end

% ODE integration options
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
Np = mpc_obj.PredictionHorizon;
target_mgdl = 100;

%% =========================================================
% CASE 1: Open-Loop Baseline (Fixed Basal Infusion)
%% =========================================================
fprintf('\nSimulating Case 1: Open-Loop Baseline...\n');
x_OL = zeros(N_steps + 1, 11);
u_OL = zeros(N_steps, 1);
x_OL(1, :) = p.ss.x0';

for k = 1:N_steps
    t_span_k = [time_vec(k), time_vec(k+1)];
    u_k_Uh = p.ss.u_basal_Uh;
    u_OL(k) = u_k_Uh;
    
    d_k_gmin = d_meal_gmin(k);
    u_k_mumin = u_k_Uh * p.UH_TO_MUMIN;
    d_k_mmolmin = d_k_gmin * (1000 / p.MOL_WT_GLUCOSE);
    
    [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
        t_span_k, x_OL(k, :)', ode_opts);
    x_OL(k+1, :) = x_seg(end, :);
end

%% =========================================================
% CASE 2: MPC Without Meal Preview (Reactive)
%% =========================================================
fprintf('Simulating Case 2: MPC Without Meal Preview (Reactive)...\n');
x_NP = zeros(N_steps + 1, 11);
u_NP = zeros(N_steps, 1);
x_NP(1, :) = p.ss.x0';

% Initialize MPC state
xmpc_NP = mpcstate(mpc_obj);

for k = 1:N_steps
    t_span_k = [time_vec(k), time_vec(k+1)];
    
    % Measured plant output (CGM glucose in mg/dL)
    y_meas_k = x_NP(k, 11) * p.MMOLL_TO_MGDL;
    
    % Current measured disturbance (no future preview passed)
    d_curr = d_meal_gmin(k);
    v_signal = d_curr; % Scalar: MPC assumes constant or zero future change
    
    % Compute optimal insulin infusion command from MPC with IOB state feedback
    xmpc_NP.Plant = x_NP(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_NP, y_meas_k, target_mgdl, v_signal);
    u_NP(k) = u_opt;
    
    % Apply to non-linear Hovorka ODE plant over [t_k, t_k+1]
    u_k_mumin = u_opt * p.UH_TO_MUMIN;
    d_k_mmolmin = d_curr * (1000 / p.MOL_WT_GLUCOSE);
    
    [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
        t_span_k, x_NP(k, :)', ode_opts);
    x_NP(k+1, :) = x_seg(end, :);
end

%% =========================================================
% CASE 3: MPC With Meal Preview (Proactive Feedforward)
%% =========================================================
fprintf('Simulating Case 3: MPC With Meal Preview (Proactive)...\n');
x_MP = zeros(N_steps + 1, 11);
u_MP = zeros(N_steps, 1);
x_MP(1, :) = p.ss.x0';

% Initialize MPC state
xmpc_MP = mpcstate(mpc_obj);

for k = 1:N_steps
    t_span_k = [time_vec(k), time_vec(k+1)];
    
    % Measured plant output (CGM glucose in mg/dL)
    y_meas_k = x_MP(k, 11) * p.MMOLL_TO_MGDL;
    
    % Build preview sequence over prediction horizon (Np+1 samples)
    v_preview = zeros(Np + 1, 1);
    for h = 0:Np
        idx_future = min(k + h, N_steps + 1);
        v_preview(h + 1) = d_meal_gmin(idx_future);
    end
    
    % Compute optimal insulin infusion command with full preview & state feedback
    xmpc_MP.Plant = x_MP(k, :)';
    u_opt = mpcmove(mpc_obj, xmpc_MP, y_meas_k, target_mgdl, v_preview);
    u_MP(k) = u_opt;
    
    % Apply to non-linear Hovorka ODE plant over [t_k, t_k+1]
    u_k_mumin = u_opt * p.UH_TO_MUMIN;
    d_k_mmolmin = d_meal_gmin(k) * (1000 / p.MOL_WT_GLUCOSE);
    
    [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
        t_span_k, x_MP(k, :)', ode_opts);
    x_MP(k+1, :) = x_seg(end, :);
end

%% =========================================================
% 4. Compute Comprehensive Clinical Metrics
%% =========================================================
% Plasma glucose in mg/dL for all cases
G_OL_mgdl = (x_OL(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_NP_mgdl = (x_NP(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
G_MP_mgdl = (x_MP(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;

% Metrics function
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
fprintf('  CLINICAL COMPARISON SUMMARY TABLE: 50g CHO MEAL AT t = 120 min\n');
fprintf('========================================================================\n');
fprintf('%-28s %-14s %-16s %-16s\n', 'Metric', 'Open-Loop', 'MPC No-Preview', 'MPC Meal-Preview');
fprintf('------------------------------------------------------------------------\n');
fprintf('%-28s %-14.2f %-16.2f %-16.2f\n', 'Peak Glucose (mg/dL)', m_OL.PeakG, m_NP.PeakG, m_MP.PeakG);
fprintf('%-28s %-14.1f %-16.1f %-16.1f\n', 'Time to Peak (min)', m_OL.t_Peak, m_NP.t_Peak, m_MP.t_Peak);
fprintf('%-28s %-14.2f %-16.2f %-16.2f\n', 'Excursion above 100 (mg/dL)', m_OL.Excursion, m_NP.Excursion, m_MP.Excursion);
fprintf('%-28s %-14.2f %-16.2f %-16.2f\n', 'Minimum Glucose (mg/dL)', m_OL.MinG, m_NP.MinG, m_MP.MinG);
fprintf('%-28s %-14.1f %-16.1f %-16.1f\n', 'Time in Range 70-180 (%)', m_OL.TIR, m_NP.TIR, m_MP.TIR);
fprintf('%-28s %-14.1f %-16.1f %-16.1f\n', 'Time above 180 (%)', m_OL.TAR, m_NP.TAR, m_MP.TAR);
fprintf('%-28s %-14.1f %-16.1f %-16.1f\n', 'Time below 70 (%)', m_OL.TBR, m_NP.TBR, m_MP.TBR);
fprintf('%-28s %-14.3f %-16.3f %-16.3f\n', 'Max Insulin Infusion (U/h)', m_OL.MaxU, m_NP.MaxU, m_MP.MaxU);
fprintf('%-28s %-14.2f %-16.2f %-16.2f\n', 'Total Insulin Delivered (U)', m_OL.TotalU, m_NP.TotalU, m_MP.TotalU);
fprintf('========================================================================\n');

%% =========================================================
% 5. Generate High-Resolution Comparison Figure
%% =========================================================
figure('Color', 'w', 'Position', [100, 80, 1000, 850]);

% Subplot 1: Blood Glucose Response
subplot(3, 1, 1);
plot(time_vec, G_OL_mgdl, 'r--', 'LineWidth', 2.0, 'DisplayName', sprintf('Open-Loop (Peak: %.1f mg/dL)', m_OL.PeakG));
hold on;
plot(time_vec, G_NP_mgdl, 'b-', 'LineWidth', 2.0, 'DisplayName', sprintf('MPC No Preview (Peak: %.1f mg/dL)', m_NP.PeakG));
plot(time_vec, G_MP_mgdl, 'g-', 'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal Preview (Peak: %.1f mg/dL)', m_MP.PeakG));
yline(100, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Target (100 mg/dL)');
yline(180, 'Color', [0.85 0.33 0.1], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hyperglycemia (180 mg/dL)');
yline(70,  'Color', [0.49 0.18 0.56], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Hypoglycemia (70 mg/dL)');
patch([0 t_end t_end 0], [70 70 180 180], [0.85 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Target Range [70-180]');
ylabel('Blood Glucose (mg/dL)', 'FontSize', 11, 'FontWeight', 'bold');
title('Hovorka Artificial Pancreas: Closed-Loop Model Predictive Control', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northeast', 'NumColumns', 2);
grid on; xlim([0, t_end]); ylim([50, 300]);

% Subplot 2: Insulin Infusion Command
subplot(3, 1, 2);
u_time = time_vec(1:end-1);
stairs(u_time, u_OL, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Open-Loop (0.39 U/h)');
hold on;
stairs(u_time, u_NP, 'b-',  'LineWidth', 2.0, 'DisplayName', sprintf('MPC No Preview (Max: %.2f U/h)', m_NP.MaxU));
stairs(u_time, u_MP, 'g-',  'LineWidth', 2.2, 'DisplayName', sprintf('MPC Meal Preview (Max: %.2f U/h)', m_MP.MaxU));
yline(p.ss.u_basal_Uh, 'k--', 'Basal (0.39 U/h)', 'LineWidth', 1.0);
yline(4.0, 'Color', [0.6 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Max Limit (4.0 U/h)');
ylabel('Insulin Infusion (U/h)', 'FontSize', 11, 'FontWeight', 'bold');
title('Insulin Pump Delivery Profile', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on; xlim([0, t_end]); ylim([0, 4.5]);

% Subplot 3: Meal Disturbance
subplot(3, 1, 3);
stem(time_vec, d_meal_gmin, 'Color', [0.85 0.33 0.1], 'LineWidth', 2.0, 'Marker', 'none', 'DisplayName', 'Meal Intake (50g CHO)');
xlabel('Simulation Time (minutes)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Meal Rate (g/min)', 'FontSize', 11, 'FontWeight', 'bold');
title('Meal Disturbance (50g CHO at t = 120 min)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast');
grid on; xlim([0, t_end]);

% Ensure target directories exist in project_root
fig_dir  = fullfile(project_root, 'results', 'figures');
data_dir = fullfile(project_root, 'results', 'data');
if ~exist(fig_dir, 'dir'),  mkdir(fig_dir);  end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

fig_path = fullfile(fig_dir, 'mpc_comparison.png');
mat_path = fullfile(data_dir, 'mpc_comparison_results.mat');

% Cleanly overwrite previous image if present
if exist(fig_path, 'file')
    try, delete(fig_path); catch, end
end

try
    exportgraphics(gcf, fig_path, 'Resolution', 300);
catch
    saveas(gcf, fig_path);
end
save(mat_path, 'time_vec', 'G_OL_mgdl', 'G_NP_mgdl', 'G_MP_mgdl', ...
    'u_OL', 'u_NP', 'u_MP', 'd_meal_gmin', 'm_OL', 'm_NP', 'm_MP', 'p');

fprintf('Comparison figure successfully saved to: %s\n', fig_path);
fprintf('Comparison data successfully saved to:   %s\n', mat_path);
fprintf('=========================================================\n');
