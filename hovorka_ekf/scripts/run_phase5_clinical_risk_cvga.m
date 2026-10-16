function [cohort_data, metrics_table] = run_phase5_clinical_risk_cvga()
% RUN_PHASE5_CLINICAL_RISK_CVGA - Phase 5 of Master Engineering Roadmap
% Standard Clinical Risk Metrics & Control-Variability Grid Analysis (CVGA)
%
% Integrates simulation results across all phases (Nominal, Sensor Noise,
% Exercise +40% SI, Sick Day -35% SI, Volume Mismatch, and 24h Multi-Meal)
% and computes consensus clinical safety metrics:
%   1. Control-Variability Grid Analysis (CVGA) 2D Risk Grid (Magni / Kovatchev)
%   2. Kovatchev Glycemic Risk Indices (LBGI, HBGI, BGRI)
%   3. Consensus International Targets (TIR, TAR L1/L2, TBR L1/L2, CV%, eA1c, GMI, J-Index)
%   4. Pairwise Clinical Risk Scorecard: Open-Loop vs EKF Output-Feedback MPC
%
% Usage:
%   run_phase5_clinical_risk_cvga()

fprintf('========================================================================\n');
fprintf('  HOVORKA ARTIFICIAL PANCREAS: PHASE 5 - CLINICAL RISK & CVGA VALIDATION\n');
fprintf('  Aggregated Cross-Scenario Cohort: Standard Clinical Benchmarking Suite\n');
fprintf('========================================================================\n');

% Path configuration
script_dir = fileparts(mfilename('fullpath'));
root_dir   = fullfile(script_dir, '..');
fig_dir    = fullfile(root_dir, 'results', 'figures', 'phase5');
data_dir   = fullfile(root_dir, 'results', 'data', 'phase5');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

%% 1. Load Simulation Data Across All Prior Phases
p1_file = fullfile(root_dir, 'results', 'data', 'phase1', 'ekf_simulation_results.mat');
p2_file = fullfile(root_dir, 'results', 'data', 'phase2', 'phase2_sensor_realism_results.mat');
p3_file = fullfile(root_dir, 'results', 'data', 'phase3', 'phase3_robustness_results.mat');
p4_file = fullfile(root_dir, 'results', 'data', 'phase4', 'phase4_24h_circadian_results.mat');

fprintf('[1/4] Ingesting clinical simulation databases from Phases 1-4...\n');
d1 = load(p1_file);
d2 = load(p2_file);
d3 = load(p3_file);
d4 = load(p4_file);

%% 2. Assemble the Cross-Scenario Cohort
% We define 6 representative daily-life clinical challenges:
%   Challenge 1: Nominal 75g Single Meal (Phase 1)
%   Challenge 2: Dexcom G7 Noise + Dropout + Lag (Phase 2)
%   Challenge 3: Post-Exercise Hypoglycemia Stress (+40% SI) (Phase 3)
%   Challenge 4: Sick Day Hyperglycemia Stress (-35% SI) (Phase 3)
%   Challenge 5: Distribution Volume Mismatch (Weight Error) (Phase 3)
%   Challenge 6: 24-Hour Circadian Multi-Meal Protocol (Phase 4)

cohort = struct();

% --- Challenge 1: Nominal 75g Meal ---
cohort(1).name     = 'Nominal Single Meal (75g)';
cohort(1).short    = 'Nominal';
cohort(1).G_ol     = d1.results.G_ol;
cohort(1).G_ekf    = d1.results.G_ekf_true;
cohort(1).G_ideal  = d1.results.G_ideal;

% --- Challenge 2: Realistic Sensor Noise (Phase 2) ---
cohort(2).name     = 'Sensor Realism (ARMA + Dropout)';
cohort(2).short    = 'Sensor Noise';
cohort(2).G_ol     = d2.results.G_ol;
cohort(2).G_ekf    = d2.results.G_blood_ekf;
cohort(2).G_ideal  = d2.results.G_ideal;

% --- Helper for Phase 3 Open-Loop Baseline Simulation ---
% Phase 3 tested EKF-MPC under parameter perturbations; let's simulate the
% exact matching Open-Loop responses for Exercise, Sick Day, and Vol Mismatch.
p_nom = hovorka_params(70, 100);
time_vec_p3 = d3.time_vec;
N_steps_p3 = length(time_vec_p3) - 1;
Ts_p3 = 5;
d_meal_p3 = zeros(N_steps_p3 + 1, 1);
for k = 1:length(time_vec_p3)
    if time_vec_p3(k) >= 30 && time_vec_p3(k) < 45
        d_meal_p3(k) = 75.0 / 15.0;
    end
end
ICR = 24.0;
u_bolus_p3 = 75.0 / ICR; % 3.125 U upfront bolus for 100% open-loop

% Function handle for simulating an open-loop perturbed run
sim_ol_perturbed = @(p_pert) simulate_openloop_case(p_pert, time_vec_p3, d_meal_p3, u_bolus_p3, Ts_p3);

% Challenge 3: Post-Exercise (+40% SI)
p_ex = p_nom;
p_ex.S_IT = p_nom.S_IT * 1.40;
p_ex.S_ID = p_nom.S_ID * 1.40;
p_ex.S_IE = p_nom.S_IE * 1.40;
p_ex.kb1  = p_ex.ka1 * p_ex.S_IT;
p_ex.kb2  = p_ex.ka2 * p_ex.S_ID;
p_ex.kb3  = p_ex.ka3 * p_ex.S_IE;
G_ol_ex = sim_ol_perturbed(p_ex);

cohort(3).name     = 'Post-Exercise (+40% S_I)';
cohort(3).short    = 'Exercise';
cohort(3).G_ol     = G_ol_ex;
cohort(3).G_ekf    = d3.all_results{2}.G_blood;
cohort(3).G_ideal  = [];

% Challenge 4: Sick Day (-35% SI)
p_sick = p_nom;
p_sick.S_IT = p_nom.S_IT * 0.65;
p_sick.S_ID = p_nom.S_ID * 0.65;
p_sick.S_IE = p_nom.S_IE * 0.65;
p_sick.kb1  = p_sick.ka1 * p_sick.S_IT;
p_sick.kb2  = p_sick.ka2 * p_sick.S_ID;
p_sick.kb3  = p_sick.ka3 * p_sick.S_IE;
G_ol_sick = sim_ol_perturbed(p_sick);

cohort(4).name     = 'Sick Day (-35% S_I)';
cohort(4).short    = 'Sick Day';
cohort(4).G_ol     = G_ol_sick;
cohort(4).G_ekf    = d3.all_results{3}.G_blood;
cohort(4).G_ideal  = [];

% Challenge 5: Distribution Volume Error (-20% VG, +20% VI)
p_vol = p_nom;
p_vol.V_G = p_nom.V_G * 0.80;
p_vol.V_I = p_nom.V_I * 1.20;
G_ol_vol = sim_ol_perturbed(p_vol);

cohort(5).name     = 'Volume Error (-20% V_G, +20% V_I)';
cohort(5).short    = 'Vol Mismatch';
cohort(5).G_ol     = G_ol_vol;
cohort(5).G_ekf    = d3.all_results{4}.G_blood;
cohort(5).G_ideal  = [];

% Challenge 6: 24-Hour Circadian Multi-Meal (Phase 4)
cohort(6).name     = '24h Circadian Multi-Meal Protocol';
cohort(6).short    = '24h Daily Life';
cohort(6).G_ol     = d4.results.G_ol;
cohort(6).G_ekf    = d4.results.G_blood_ekf;
cohort(6).G_ideal  = d4.results.G_ideal;

N_cohort = length(cohort);

%% 3. Compute Standard Clinical Risk Metrics for All Cohorts
fprintf('[2/4] Calculating Kovatchev Glycemic Risk Indices & Consensus Metrics...\n');

clear metrics_ol metrics_ekf;

for c = 1:N_cohort
    metrics_ol(c)  = compute_full_clinical_metrics(cohort(c).G_ol);
    metrics_ekf(c) = compute_full_clinical_metrics(cohort(c).G_ekf);
end

% Print Comprehensive Clinical Risk Summary Table
fprintf('\n========================================================================================================\n');
fprintf('  PHASE 5 CLINICAL SCORECARD: KOVATCHEV RISK INDICES & CONSENSUS BENCHMARKS\n');
fprintf('========================================================================================================\n');
fprintf('%-24s | %-12s | %-7s %-7s %-7s | %-6s %-6s %-6s | %-6s %-6s\n', ...
    'Clinical Challenge', 'Mode', 'LBGI', 'HBGI', 'BGRI', 'TIR%', 'TAR%', 'TBR%', 'eA1c%', '%CV');
fprintf('--------------------------------------------------------------------------------------------------------\n');

for c = 1:N_cohort
    m_ol  = metrics_ol(c);
    m_ekf = metrics_ekf(c);
    fprintf('%-24s | Open-Loop    | %7.2f %7.2f %7.2f | %6.1f %6.1f %6.1f | %6.2f %6.1f\n', ...
        cohort(c).name, m_ol.LBGI, m_ol.HBGI, m_ol.BGRI, m_ol.TIR, m_ol.TAR_L1+m_ol.TAR_L2, m_ol.TBR_L1+m_ol.TBR_L2, m_ol.eA1c, m_ol.CV);
    fprintf('%-24s | EKF-MPC (CL) | %7.2f %7.2f %7.2f | %6.1f %6.1f %6.1f | %6.2f %6.1f\n', ...
        '', m_ekf.LBGI, m_ekf.HBGI, m_ekf.BGRI, m_ekf.TIR, m_ekf.TAR_L1+m_ekf.TAR_L2, m_ekf.TBR_L1+m_ekf.TBR_L2, m_ekf.eA1c, m_ekf.CV);
    fprintf('--------------------------------------------------------------------------------------------------------\n');
end

% Cohort Mean Averages
avg_ol_lbgi = mean([metrics_ol.LBGI]); avg_ol_hbgi = mean([metrics_ol.HBGI]); avg_ol_bgri = mean([metrics_ol.BGRI]);
avg_ol_tir  = mean([metrics_ol.TIR]);  avg_ol_tar  = mean([metrics_ol.TAR_L1] + [metrics_ol.TAR_L2]);
avg_ol_tbr  = mean([metrics_ol.TBR_L1] + [metrics_ol.TBR_L2]); avg_ol_cv = mean([metrics_ol.CV]);

avg_ekf_lbgi = mean([metrics_ekf.LBGI]); avg_ekf_hbgi = mean([metrics_ekf.HBGI]); avg_ekf_bgri = mean([metrics_ekf.BGRI]);
avg_ekf_tir  = mean([metrics_ekf.TIR]);  avg_ekf_tar  = mean([metrics_ekf.TAR_L1] + [metrics_ekf.TAR_L2]);
avg_ekf_tbr  = mean([metrics_ekf.TBR_L1] + [metrics_ekf.TBR_L2]); avg_ekf_cv = mean([metrics_ekf.CV]);

fprintf('COHORT AVERAGE:          | Open-Loop    | %7.2f %7.2f %7.2f | %6.1f %6.1f %6.1f | %6.2f %6.1f\n', ...
    avg_ol_lbgi, avg_ol_hbgi, avg_ol_bgri, avg_ol_tir, avg_ol_tar, avg_ol_tbr, mean([metrics_ol.eA1c]), avg_ol_cv);
fprintf('COHORT AVERAGE:          | EKF-MPC (CL) | %7.2f %7.2f %7.2f | %6.1f %6.1f %6.1f | %6.2f %6.1f\n', ...
    avg_ekf_lbgi, avg_ekf_hbgi, avg_ekf_bgri, avg_ekf_tir, avg_ekf_tar, avg_ekf_tbr, mean([metrics_ekf.eA1c]), avg_ekf_cv);
fprintf('========================================================================================================\n');
fprintf('  CLINICAL RISK REDUCTION SUMMARY:\n');
fprintf('  * Overall Blood Glucose Risk Index (BGRI): Dropped from %.2f (High) to %.2f (Low) (-%.1f%% reduction)\n', ...
    avg_ol_bgri, avg_ekf_bgri, 100*(avg_ol_bgri - avg_ekf_bgri)/avg_ol_bgri);
fprintf('  * High Blood Glucose Index (HBGI):         Dropped from %.2f to %.2f (-%.1f%% reduction)\n', ...
    avg_ol_hbgi, avg_ekf_hbgi, 100*(avg_ol_hbgi - avg_ekf_hbgi)/avg_ol_hbgi);
fprintf('  * Cohort Mean Time-in-Range (TIR):          Boosted from %.1f%% to %.1f%% (+%.1f%% absolute gain)\n', ...
    avg_ol_tir, avg_ekf_tir, avg_ekf_tir - avg_ol_tir);
fprintf('  * Cohort Glycemic Variability (CV%%):        Decreased from %.1f%% to %.1f%%\n', ...
    avg_ol_cv, avg_ekf_cv);
fprintf('========================================================================================================\n\n');

%% 4. Generate Publication-Quality Figures

fprintf('[3/4] Generating Publication Figures (Interactive Display Mode)...\n');

% -------------------------------------------------------------------------
% FIGURE 1: Control-Variability Grid Analysis (CVGA) 2D Grid Plot
% -------------------------------------------------------------------------
fig1 = figure(1501);
set(fig1, 'Name', 'Phase 5: Control-Variability Grid Analysis (CVGA)', ...
          'Color', 'w', 'Position', [80, 80, 1100, 850], 'Visible', 'on');
clf(fig1);
ax1 = axes('Parent', fig1);
hold(ax1, 'on');

% Define CVGA coordinate space
x_lims = [40, 130];
y_lims = [90, 320];

% Draw Standard CVGA Zones with Semi-Transparent Alpha Patches
% Zone A: Accurate Control / Optimum (Green) [Gmin: 90-110, Gmax: 110-180] extended to Gmin >= 90
patch([90, 130, 130, 90], [90, 90, 180, 180], [0.55 0.90 0.55], 'FaceAlpha', 0.45, 'EdgeColor', [0.2 0.7 0.2], 'LineWidth', 1.5);
text(107, 135, sprintf('ZONE A\n(Optimal Control)'), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.0 0.45 0.0]);

% Lower Zone B: Benign Hypoglycemia [Gmin: 70-90, Gmax: 90-180] (Light Green/Yellow)
patch([70, 90, 90, 70], [90, 90, 180, 180], [0.85 0.95 0.55], 'FaceAlpha', 0.35, 'EdgeColor', [0.6 0.8 0.2], 'LineWidth', 1.2);
text(80, 135, sprintf('ZONE B\n(Benign Hypo)'), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.35 0.5 0.1]);

% Upper Zone B: Benign Hyperglycemia [Gmin: 90-130, Gmax: 180-260]
patch([90, 130, 130, 90], [180, 180, 260, 260], [0.95 0.95 0.55], 'FaceAlpha', 0.35, 'EdgeColor', [0.8 0.8 0.2], 'LineWidth', 1.2);
text(110, 220, sprintf('ZONE B\n(Benign Hyper)'), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.1]);

% Bivariate Zone B: [Gmin: 70-90, Gmax: 180-260]
patch([70, 90, 90, 70], [180, 180, 260, 260], [1.0 0.92 0.60], 'FaceAlpha', 0.35, 'EdgeColor', [0.9 0.7 0.2], 'LineWidth', 1.2);
text(80, 220, sprintf('ZONE B\n(Bivariate Dev)'), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.6 0.4 0.1]);

% Lower Zone C: Over-Correction / Hypoglycemia Risk [Gmin: 54-70, Gmax: 90-180] (Orange)
patch([54, 70, 70, 54], [90, 90, 180, 180], [1.0 0.75 0.45], 'FaceAlpha', 0.35, 'EdgeColor', [0.9 0.5 0.1], 'LineWidth', 1.2);
text(62, 135, sprintf('ZONE C\n(Over-Correction)'), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.7 0.3 0.0]);

% Upper Zone C: Under-Correction / Hyperglycemia Risk [Gmin: 90-130, Gmax: 260-320]
patch([90, 130, 130, 90], [260, 260, 320, 320], [1.0 0.75 0.45], 'FaceAlpha', 0.35, 'EdgeColor', [0.9 0.5 0.1], 'LineWidth', 1.2);
text(110, 290, sprintf('ZONE C\n(Under-Correction)'), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.7 0.3 0.0]);

% Zone D: Failure to Control [Gmin: 54-70, Gmax: 180-260] (Light Red)
patch([54, 70, 70, 54], [180, 180, 260, 260], [1.0 0.55 0.55], 'FaceAlpha', 0.35, 'EdgeColor', [0.8 0.2 0.2], 'LineWidth', 1.2);
text(62, 220, sprintf('ZONE D\n(Control Failure)'), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.65 0.1 0.1]);

% Zone E: Extreme Hypo Hazard [Gmin: 40-54, Gmax: 90-320] (Dark Red)
patch([40, 54, 54, 40], [90, 90, 320, 320], [0.95 0.35 0.35], 'FaceAlpha', 0.40, 'EdgeColor', [0.6 0.1 0.1], 'LineWidth', 1.5);
text(47, 205, sprintf('ZONE E\n(Severe Hypo\nHazard < 54)'), 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.5 0.0 0.0]);

% Zone E: Extreme Hyper Hazard [Gmax: 300-320, Gmin: 54-90]
patch([54, 90, 90, 54], [260, 260, 320, 320], [0.95 0.35 0.35], 'FaceAlpha', 0.35, 'EdgeColor', [0.6 0.1 0.1], 'LineWidth', 1.2);
text(72, 290, sprintf('ZONE E\n(Extreme Hyper)'), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', [0.5 0.0 0.0]);

% Plot Migration Arrows from Open-Loop to EKF-MPC
colors = lines(N_cohort);
h_ol_plot = [];
h_ekf_plot = [];

for c = 1:N_cohort
    x_ol  = metrics_ol(c).Nadir;
    y_ol  = metrics_ol(c).Peak;
    x_ekf = metrics_ekf(c).Nadir;
    y_ekf = metrics_ekf(c).Peak;
    
    % Draw transition vector
    quiver(ax1, x_ol, y_ol, x_ekf - x_ol, y_ekf - y_ol, 0, ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1.8, 'MaxHeadSize', 0.6, 'LineStyle', '--');
    
    % Open-Loop marker (Square, hollow)
    h_ol = plot(ax1, x_ol, y_ol, 's', 'MarkerSize', 11, 'MarkerFaceColor', [0.85 0.2 0.2], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % EKF-MPC marker (Circle, filled)
    h_ekf = plot(ax1, x_ekf, y_ekf, 'o', 'MarkerSize', 11, 'MarkerFaceColor', [0.1 0.55 0.95], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Annotate scenario labels
    text(x_ekf + 1.2, y_ekf - 2.0, sprintf('(%d) %s', c, cohort(c).short), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.1 0.2 0.5]);
    
    if c == 1
        h_ol_plot = h_ol;
        h_ekf_plot = h_ekf;
    end
end

grid(ax1, 'on');
box(ax1, 'on');
set(ax1, 'XLim', x_lims, 'YLim', y_lims, 'FontSize', 11);
xlabel(ax1, 'Minimum Blood Glucose G_{min} (mg/dL)  [Hypoglycemia Risk \rightarrow]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax1, '[\leftarrow Hyperglycemia Risk]  Maximum Blood Glucose G_{max} (mg/dL)', 'FontSize', 12, 'FontWeight', 'bold');
title(ax1, 'Control-Variability Grid Analysis (CVGA): Open-Loop vs EKF Output-Feedback MPC', 'FontSize', 14, 'FontWeight', 'bold');

legend(ax1, [h_ol_plot, h_ekf_plot], ...
    {'Open-Loop Baselines (MDI Therapy)', 'Proposed EKF-MPC (Closed-Loop AID)'}, ...
    'Location', 'northwest', 'FontSize', 11);

% Save Figure 1
fig1_path = fullfile(fig_dir, 'phase5_cvga_risk_grid.png');
drawnow;
if exist(fig1_path, 'file'), try delete(fig1_path); catch; end; end
exportgraphics(fig1, fig1_path, 'Resolution', 300);
fprintf('  Saved Figure 1: %s\n', fig1_path);


% -------------------------------------------------------------------------
% FIGURE 2: Kovatchev Glycemic Risk Indices (LBGI & HBGI) Bar Comparison
% -------------------------------------------------------------------------
fig2 = figure(1502);
set(fig2, 'Name', 'Phase 5: Kovatchev Glycemic Risk Indices (LBGI & HBGI)', ...
          'Color', 'w', 'Position', [120, 120, 1100, 700], 'Visible', 'on');
clf(fig2);

% Subplot 1: Low Blood Glucose Index (Hypo Risk)
subplot(1, 2, 1);
lbgi_data = [[metrics_ol.LBGI]', [metrics_ekf.LBGI]'];
b1 = bar(lbgi_data, 'grouped');
b1(1).FaceColor = [0.85 0.3 0.3];
b1(2).FaceColor = [0.2 0.6 0.9];
hold on;
yline(1.1, 'k--', 'Minimal Risk Threshold (< 1.1)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(2.5, 'm:', 'Low Risk Boundary (1.1 - 2.5)', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
yline(5.0, 'r--', 'High Risk Boundary (> 5.0)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
set(gca, 'XTick', 1:N_cohort, 'XTickLabel', {cohort.short}, 'XTickLabelRotation', 30, 'FontSize', 10);
ylabel('Low Blood Glucose Index (LBGI)', 'FontSize', 12, 'FontWeight', 'bold');
title('Hypoglycemia Risk (LBGI)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Open-Loop Baseline', 'Proposed EKF-MPC'}, 'Location', 'northwest', 'FontSize', 10);
grid on; box on;
ylim([0, max(max(lbgi_data)) * 1.25 + 0.5]);

% Subplot 2: High Blood Glucose Index (Hyper Risk)
subplot(1, 2, 2);
hbgi_data = [[metrics_ol.HBGI]', [metrics_ekf.HBGI]'];
b2 = bar(hbgi_data, 'grouped');
b2(1).FaceColor = [0.9 0.4 0.2];
b2(2).FaceColor = [0.15 0.7 0.3];
hold on;
yline(5.0, 'k--', 'Minimal Risk Threshold (< 5.0)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(10.0, 'm:', 'Low Risk Boundary (5.0 - 10.0)', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
yline(15.0, 'r--', 'High Risk Boundary (> 15.0)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
set(gca, 'XTick', 1:N_cohort, 'XTickLabel', {cohort.short}, 'XTickLabelRotation', 30, 'FontSize', 10);
ylabel('High Blood Glucose Index (HBGI)', 'FontSize', 12, 'FontWeight', 'bold');
title('Hyperglycemia Risk (HBGI)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Open-Loop Baseline', 'Proposed EKF-MPC'}, 'Location', 'northeast', 'FontSize', 10);
grid on; box on;
ylim([0, max(max(hbgi_data)) * 1.25 + 1.0]);

sgtitle('Kovatchev Glycemic Risk Indices Across Stress Cohort', 'FontSize', 15, 'FontWeight', 'bold');

% Save Figure 2
fig2_path = fullfile(fig_dir, 'phase5_kovatchev_risk_indices.png');
drawnow;
if exist(fig2_path, 'file'), try delete(fig2_path); catch; end; end
exportgraphics(fig2, fig2_path, 'Resolution', 300);
fprintf('  Saved Figure 2: %s\n', fig2_path);


% -------------------------------------------------------------------------
% FIGURE 3: International Consensus Clinical Metrics Scorecard
% -------------------------------------------------------------------------
fig3 = figure(1503);
set(fig3, 'Name', 'Phase 5: International Clinical Consensus Metrics Scorecard', ...
          'Color', 'w', 'Position', [160, 160, 1150, 750], 'Visible', 'on');
clf(fig3);

% Panel 1: Time in Range (TIR 70-180 mg/dL)
subplot(2, 2, 1);
tir_data = [[metrics_ol.TIR]', [metrics_ekf.TIR]'];
b_tir = bar(tir_data, 'grouped');
b_tir(1).FaceColor = [0.8 0.4 0.4];
b_tir(2).FaceColor = [0.2 0.75 0.3];
hold on;
yline(70, 'r--', 'International Target (TIR > 70%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
set(gca, 'XTick', 1:N_cohort, 'XTickLabel', {cohort.short}, 'FontSize', 9);
ylabel('Time in Range (%)', 'FontSize', 11, 'FontWeight', 'bold');
title('Time-in-Range (70 - 180 mg/dL)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0, 115]); grid on; box on;
legend({'Open-Loop', 'EKF-MPC'}, 'Location', 'southwest', 'FontSize', 9);

% Panel 2: Glycemic Variability (%CV)
subplot(2, 2, 2);
cv_data = [[metrics_ol.CV]', [metrics_ekf.CV]'];
b_cv = bar(cv_data, 'grouped');
b_cv(1).FaceColor = [0.85 0.5 0.2];
b_cv(2).FaceColor = [0.2 0.55 0.85];
hold on;
yline(36, 'k--', 'Clinical Stability Target (\leq 36%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
set(gca, 'XTick', 1:N_cohort, 'XTickLabel', {cohort.short}, 'FontSize', 9);
ylabel('Coefficient of Variation (%CV)', 'FontSize', 11, 'FontWeight', 'bold');
title('Glycemic Variability (%CV = \sigma / \mu)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0, max(max(cv_data)) * 1.3]); grid on; box on;
legend({'Open-Loop', 'EKF-MPC'}, 'Location', 'northwest', 'FontSize', 9);

% Panel 3: Composite Blood Glucose Risk Index (BGRI = LBGI + HBGI)
subplot(2, 2, 3);
bgri_data = [[metrics_ol.BGRI]', [metrics_ekf.BGRI]'];
b_bgri = bar(bgri_data, 'grouped');
b_bgri(1).FaceColor = [0.75 0.25 0.4];
b_bgri(2).FaceColor = [0.25 0.7 0.6];
hold on;
yline(6.1, 'g--', 'Minimal Risk (\leq 6.1)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(12.5, 'm:', 'Moderate Risk (> 12.5)', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
set(gca, 'XTick', 1:N_cohort, 'XTickLabel', {cohort.short}, 'FontSize', 9);
ylabel('Composite Risk Index (BGRI)', 'FontSize', 11, 'FontWeight', 'bold');
title('Total Glycemic Risk (BGRI = LBGI + HBGI)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0, max(max(bgri_data)) * 1.3]); grid on; box on;
legend({'Open-Loop', 'EKF-MPC'}, 'Location', 'northwest', 'FontSize', 9);

% Panel 4: Estimated HbA1c (eA1c %)
subplot(2, 2, 4);
ea1c_data = [[metrics_ol.eA1c]', [metrics_ekf.eA1c]'];
b_ea1c = bar(ea1c_data, 'grouped');
b_ea1c(1).FaceColor = [0.9 0.4 0.2];
b_ea1c(2).FaceColor = [0.3 0.65 0.3];
hold on;
yline(7.0, 'r--', 'ADA Target (< 7.0%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(5.7, 'k:', 'Normal Glucose (< 5.7%)', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
set(gca, 'XTick', 1:N_cohort, 'XTickLabel', {cohort.short}, 'FontSize', 9);
ylabel('Estimated HbA1c (%)', 'FontSize', 11, 'FontWeight', 'bold');
title('Estimated Long-Term Glycemia (eA1c)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([4.5, 7.5]); grid on; box on;
legend({'Open-Loop', 'EKF-MPC'}, 'Location', 'southwest', 'FontSize', 9);

sgtitle('International Clinical Consensus Scorecard: Cross-Cohort Evaluation', 'FontSize', 15, 'FontWeight', 'bold');

% Save Figure 3
fig3_path = fullfile(fig_dir, 'phase5_consensus_metrics_scorecard.png');
drawnow;
if exist(fig3_path, 'file'), try delete(fig3_path); catch; end; end
exportgraphics(fig3, fig3_path, 'Resolution', 300);
fprintf('  Saved Figure 3: %s\n', fig3_path);

%% 5. Save Aggregated Phase 5 Results Dataset
fprintf('[4/4] Archiving Phase 5 clinical risk dataset...\n');
metrics_table = struct('cohort', cohort, 'metrics_ol', metrics_ol, 'metrics_ekf', metrics_ekf);
save_file = fullfile(data_dir, 'phase5_clinical_risk_results.mat');
save(save_file, 'cohort', 'metrics_ol', 'metrics_ekf', 'metrics_table');
fprintf('Phase 5 clinical risk dataset saved to: %s\n', save_file);
fprintf('========================================================================\n\n');

cohort_data = cohort;
end

%% ========================================================================
% HELPER FUNCTIONS
% =========================================================================

function m = compute_full_clinical_metrics(G_mgdl)
% COMPUTE_FULL_CLINICAL_METRICS - Complete Kovatchev and International Consensus Metrics
    G = max(G_mgdl, 1.0); % avoid negative/zero in logarithms
    N = length(G);
    
    % Kovatchev Symmetric Transform
    f = 1.509 * ((log(G)).^1.084 - 5.381);
    r = 10 * (f.^2);
    rl = zeros(size(r));
    rh = zeros(size(r));
    rl(f < 0) = r(f < 0);
    rh(f > 0) = r(f > 0);
    
    m.LBGI = mean(rl);
    m.HBGI = mean(rh);
    m.BGRI = m.LBGI + m.HBGI;
    
    % Basic Statistics
    m.Mean = mean(G);
    m.Std  = std(G);
    m.CV   = (m.Std / m.Mean) * 100;
    m.Peak = max(G);
    m.Nadir = min(G);
    
    % International Consensus Time-in-Ranges
    m.TIR    = (sum(G >= 70 & G <= 180) / N) * 100;
    m.TAR_L1 = (sum(G > 180 & G <= 250) / N) * 100;
    m.TAR_L2 = (sum(G > 250) / N) * 100;
    m.TBR_L1 = (sum(G >= 54 & G < 70) / N) * 100;
    m.TBR_L2 = (sum(G < 54) / N) * 100;
    
    % Long-Term Glycemia Estimates
    m.eA1c    = (m.Mean + 46.7) / 28.7;
    m.GMI     = 3.31 + 0.02392 * m.Mean;
    m.J_Index = 0.00106 * (m.Mean + m.Std)^2;
end

function G_blood = simulate_openloop_case(p_patient, time_vec, d_meal_vec, u_bolus, Ts)
% SIMULATE_OPENLOOP_CASE - Simulates open-loop standard MDI therapy under parameter mismatch
    N_steps = length(time_vec) - 1;
    x = p_patient.ss.x0;
    u_basal = p_patient.ss.u_basal_Uh;
    G_blood = zeros(N_steps + 1, 1);
    G_blood(1) = (x(7) / p_patient.V_G) * p_patient.MMOLL_TO_MGDL;
    
    bolus_speed_Uh = 150.0; % High-speed commercial pump bolus rate (150 U/h)
    t_bolus_min = (u_bolus / bolus_speed_Uh) * 60;
    ode_opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
    
    for k = 1:N_steps
        t_k = time_vec(k);
        t_next = time_vec(k+1);
        d_k_gmin = d_meal_vec(k);
        d_k_mmolmin = d_k_gmin * (1000 / p_patient.MOL_WT_GLUCOSE);
        
        if t_k == 30 % Meal start time
            t_split = t_k + t_bolus_min;
            u_bolus_rate_mumin = (u_basal + bolus_speed_Uh) * p_patient.UH_TO_MUMIN;
            [~, s1] = ode15s(@(t, x_cur) hovorka_ode(t, x_cur, u_bolus_rate_mumin, d_k_mmolmin, p_patient), ...
                [t_k, t_split], x, ode_opts);
            u_basal_mumin = u_basal * p_patient.UH_TO_MUMIN;
            [~, s2] = ode15s(@(t, x_cur) hovorka_ode(t, x_cur, u_basal_mumin, d_k_mmolmin, p_patient), ...
                [t_split, t_next], s1(end, :)', ode_opts);
            x = s2(end, :)';
        else
            u_k_mumin = u_basal * p_patient.UH_TO_MUMIN;
            [~, s] = ode15s(@(t, x_cur) hovorka_ode(t, x_cur, u_k_mumin, d_k_mmolmin, p_patient), ...
                [t_k, t_next], x, ode_opts);
            x = s(end, :)';
        end
        G_blood(k+1) = (x(7) / p_patient.V_G) * p_patient.MMOLL_TO_MGDL;
    end
end
