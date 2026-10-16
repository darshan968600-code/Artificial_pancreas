%% RUN_ALL_BERGMAN_SIMULATIONS Master simulation suite for Bergman Model
% Recreates the complete set of experiments from the original Artificial_pancreas project:
%   1. Normal Meal (300 min): Open-Loop, MPC No-Preview, MPC Meal-Preview
%   2. Large Meal (600 min): Open-Loop, MPC No-Preview, MPC Meal-Preview
%   3. Extended Normal Meal (600 min): MPC No-Preview, MPC Meal-Preview
%
% Saves results to results/data/AP_simulation_data.mat with identical variable structure.

clc; clear; close all;
fprintf('========================================================================\n');
fprintf('  RUNNING COMPLETE BERGMAN ARTIFICIAL PANCREAS SIMULATION SUITE\n');
fprintf('========================================================================\n');

% 1. Dynamic path setup
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(this_dir); % C:\Pancreas_project\bergman_model
model_dir = fullfile(project_root, 'models');
data_dir  = fullfile(project_root, 'results', 'data');
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

addpath(fullfile(project_root, 'scripts'));
addpath(model_dir);

% 2. Load Controllers into Base Workspace
load(fullfile(model_dir, 'mpc_baseline.mat'), 'mpcobj');
load(fullfile(model_dir, 'mpc_meal_preview.mat'), 'mpcobj_md');
assignin('base', 'mpcobj', mpcobj);
assignin('base', 'mpcobj_md', mpcobj_md);
evalin('base', 'clear out');

% Initialize results structure
results = struct();
results.normal   = struct();
results.large    = struct();
results.extended = struct();

stopTimeNormal   = '18000'; % 300 minutes (in seconds)
stopTimeLarge    = '36000'; % 600 minutes (in seconds)
stopTimeExtended = '36000'; % 600 minutes (in seconds)

%% =====================================================================
% PART 1: NORMAL MEAL SIMULATIONS (300 min / 18,000 sec)
%% =====================================================================
fprintf('\n[1/3] Running Normal-Meal Simulations (D_peak = 0.6 mg/dL/min)...\n');

% 1. Open-Loop Normal Meal
model = 'Artificial_Pancreas_OpenLoop_Comparison';
fprintf('  Simulating %s...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeNormal);
out = sim(model);
results.normal.openLoop.glucose = out.logsout.getElement(1).Values;
results.normal.openLoop.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

% 2. MPC No-Preview Normal Meal
model = 'Artificial_Pancreas_MPC_NoPreview';
fprintf('  Simulating %s...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeNormal);
out = sim(model);
results.normal.mpcNoPreview.glucose = out.logsout.getElement(1).Values;
results.normal.mpcNoPreview.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

% 3. MPC With Meal Preview Normal Meal
model = 'Artificial_Pancreas_MPC_MealPreview';
fprintf('  Simulating %s...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeNormal);
out = sim(model);
results.normal.mpcPreview.glucose = out.logsout.getElement(1).Values;
results.normal.mpcPreview.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

%% =====================================================================
% PART 2: LARGE MEAL SIMULATIONS (600 min / 36,000 sec)
%% =====================================================================
fprintf('\n[2/3] Running Large-Meal Simulations (D_peak = 1.2 mg/dL/min)...\n');

% 1. Open-Loop Large Meal
model = 'Artificial_Pancreas_OpenLoop_LargeMeal';
fprintf('  Simulating %s...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeLarge);
out = sim(model);
results.large.openLoop.glucose = out.logsout.getElement(1).Values;
results.large.openLoop.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

% 2. MPC No-Preview Large Meal
model = 'Artificial_Pancreas_MPC_NoPreview_LargeMeal';
fprintf('  Simulating %s...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeLarge);
out = sim(model);
results.large.mpcNoPreview.glucose = out.logsout.getElement(1).Values;
results.large.mpcNoPreview.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

% 3. MPC With Meal Preview Large Meal
model = 'Artificial_Pancreas_MPC_MealPreview_LargeMeal';
fprintf('  Simulating %s...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeLarge);
out = sim(model);
results.large.mpcPreview.glucose = out.logsout.getElement(1).Values;
results.large.mpcPreview.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

%% =====================================================================
% PART 3: EXTENDED 600-MINUTE RECOVERY SIMULATIONS
%% =====================================================================
fprintf('\n[3/3] Running Extended 600-minute Simulations...\n');

% 1. Extended No-Preview
model = 'Artificial_Pancreas_MPC_NoPreview';
fprintf('  Simulating %s (Extended 600 min)...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeExtended);
out = sim(model);
results.extended.mpcNoPreview.glucose = out.logsout.getElement(1).Values;
results.extended.mpcNoPreview.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

% 2. Extended With Preview
model = 'Artificial_Pancreas_MPC_MealPreview';
fprintf('  Simulating %s (Extended 600 min)...\n', model);
load_system(fullfile(model_dir, [model, '.slx']));
set_param(model, 'StopTime', stopTimeExtended);
out = sim(model);
results.extended.mpcPreview.glucose = out.logsout.getElement(1).Values;
results.extended.mpcPreview.insulin = out.logsout.getElement(2).Values;
close_system(model, 0);

% Save results structure
mat_file = fullfile(data_dir, 'AP_simulation_data.mat');
save(mat_file, 'results');
fprintf('\nAll Bergman simulations completed successfully!\n');
fprintf('Data saved to:\n  %s\n', mat_file);
fprintf('========================================================================\n');
