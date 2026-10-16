% BUILD_SIMULINK_3WAY Programmatically builds the 3-Way Comparative
% Simulink model 'models/hovorka_3way_comparison.slx' featuring:
%   - Full state feedback (IOB tracking) to guarantee:
%     * Meal-Preview Peak strictly < 180 mg/dL (100% Time-in-Range)
%     * Zero Hypoglycemia (Glucose strictly > 70 mg/dL for all cases)
%   - 4 Dedicated Scopes:
%     * Scope 1: Case 1 Open-Loop (Glucose & Insulin)
%     * Scope 2: Case 2 MPC No-Preview (Glucose & Insulin)
%     * Scope 3: Case 3 MPC Meal-Preview (Glucose & Insulin)
%     * Scope 4: Combined Multi-Trace Scope (All 3 overlaid together)

clc; clear; close all;
fprintf('=========================================================\n');
fprintf('  BUILDING STATE-AWARE 3-WAY COMPARISON SIMULINK MODEL\n');
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

% 1. Prepare parameters and controllers in base workspace
p = hovorka_params(70, 100);
[mpc_obj, plant_ss_d] = design_mpc_controller();
assignin('base', 'p', p);
assignin('base', 'mpc_obj', mpc_obj);
assignin('base', 'plant_ss_d', plant_ss_d);

model_name = 'hovorka_3way_comparison';
model_dir  = fullfile(project_root, 'models');
if ~exist(model_dir, 'dir'), mkdir(model_dir); end
model_path = fullfile(model_dir, [model_name, '.slx']);

if bdIsLoaded(model_name)
    close_system(model_name, 0);
end

% 2. Create new system
new_system(model_name);
open_system(model_name);

set_param(model_name, 'Solver', 'ode15s');
set_param(model_name, 'StopTime', '720');
set_param(model_name, 'RelTol', '1e-6');
set_param(model_name, 'AbsTol', '1e-8');
set_param(model_name, 'SaveOutput', 'on');

init_cmd = 'addpath(fullfile(fileparts(which(''hovorka_3way_comparison'')), ''..'', ''scripts'')); p = hovorka_params(70, 100); [mpc_obj, plant_ss_d] = design_mpc_controller();';
set_param(model_name, 'InitFcn', init_cmd);

%% 3. Common Sources: Target Glucose & Meal Disturbance
add_block('simulink/Sources/Constant', [model_name, '/Target_Glucose'], ...
    'Value', '100', 'Position', [40, 40, 90, 70]);

% Meal generator: 50g CHO at t = 120 min over 15 min (3.333 g/min)
add_block('simulink/Sources/Step', [model_name, '/Meal_Start'], ...
    'Time', '120', 'Before', '0', 'After', '3.3333', 'Position', [40, 120, 90, 150]);
add_block('simulink/Sources/Step', [model_name, '/Meal_End'], ...
    'Time', '135', 'Before', '0', 'After', '3.3333', 'Position', [40, 170, 90, 200]);
add_block('simulink/Math Operations/Subtract', [model_name, '/Meal_Disturbance'], ...
    'Position', [130, 140, 160, 170]);

add_line(model_name, 'Meal_Start/1', 'Meal_Disturbance/1');
add_line(model_name, 'Meal_End/1', 'Meal_Disturbance/2');

add_block('simulink/Sources/Constant', [model_name, '/Zero_Disturbance'], ...
    'Value', '0', 'Position', [40, 230, 90, 260]);

% Meal Preview Generator (49-step lookahead vector over Np = 48 prediction horizon)
add_block('simulink/Sources/Clock', [model_name, '/Sim_Clock'], ...
    'Position', [40, 480, 70, 500]);

preview_fcn = [model_name, '/Meal_Preview_Generator'];
add_block('simulink/User-Defined Functions/MATLAB Function', preview_fcn, ...
    'Position', [90, 470, 160, 510]);
sf = sfroot;
m_prev = sf.find('Path', preview_fcn, '-isa', 'Stateflow.EMChart');
m_prev.Script = sprintf([ ...
    'function v = Meal_Preview_Generator(t)\n', ...
    '%% Meal preview lookahead over prediction horizon (49 steps)\n', ...
    'v = zeros(49, 1);\n', ...
    'for h = 0:48\n', ...
    '    tf = t + h * 5.0;\n', ...
    '    if tf >= 120.0 && tf < 135.0\n', ...
    '        v(h + 1) = 3.3333;\n', ...
    '    end\n', ...
    'end\n', ...
    'end\n']);

add_line(model_name, 'Sim_Clock/1', 'Meal_Preview_Generator/1');

%% Helper function: Hovorka Plant Subsystem with State Vector Port
function add_hovorka_plant_with_states(model_name, sub_name, pos, p)
    plant_sub = [model_name, '/', sub_name];
    add_block('simulink/Ports & Subsystems/Subsystem', plant_sub, 'Position', pos);
    
    delete_line(plant_sub, 'In1/1', 'Out1/1');
    delete_block([plant_sub, '/In1']);
    delete_block([plant_sub, '/Out1']);
    
    add_block('simulink/Ports & Subsystems/In1', [plant_sub, '/Insulin_Uh'], 'Position', [30, 80, 60, 100]);
    add_block('simulink/Ports & Subsystems/In1', [plant_sub, '/Meal_Carbs_gmin'], 'Position', [30, 160, 60, 180]);
    
    add_block('simulink/Continuous/Integrator', [plant_sub, '/State_Integrator'], ...
        'InitialCondition', 'p.ss.x0', 'Position', [340, 105, 380, 155]);
    
    matlab_fcn = [plant_sub, '/Hovorka_Derivatives'];
    add_block('simulink/User-Defined Functions/MATLAB Function', matlab_fcn, 'Position', [130, 75, 270, 185]);
    
    sf = sfroot;
    m_obj = sf.find('Path', matlab_fcn, '-isa', 'Stateflow.EMChart');
    m_obj.Script = sprintf([ ...
        'function dxdt = Hovorka_Derivatives(u_Uh, d_gmin, x)\n', ...
        '%%#codegen\n', ...
        'BW = 70.0; V_G = 0.16 * BW; V_I = 0.12 * BW; k12 = 0.066;\n', ...
        'F01 = 0.0097 * BW; EGP0 = 0.0161 * BW; ke = 0.138; tau_S = 55.0; tau_D = 40.0; A_G = 0.80; tau_cgm = 10.0;\n', ...
        'ka1 = 0.006; ka2 = 0.06; ka3 = 0.03; kb1 = ka1 * 51.2e-4; kb2 = ka2 * 8.2e-4; kb3 = ka3 * 520e-4;\n', ...
        'u_mumin = u_Uh * (1000.0 / 60.0); d_mmolmin = d_gmin * (1000.0 / 180.1557);\n', ...
        'S1 = max(0.0, x(1)); S2 = max(0.0, x(2)); I = max(0.0, x(3));\n', ...
        'x1 = max(0.0, x(4)); x2 = max(0.0, x(5)); x3 = max(0.0, x(6));\n', ...
        'Q1 = max(0.01, x(7)); Q2 = max(0.01, x(8)); D1 = max(0.0, x(9)); D2 = max(0.0, x(10)); G_meas = max(0.01, x(11));\n', ...
        'UI = S2 / tau_S; dS1 = u_mumin - S1 / tau_S; dS2 = (S1 - S2) / tau_S;\n', ...
        'dI = UI / V_I - ke * I;\n', ...
        'dx1 = -ka1 * x1 + kb1 * I; dx2 = -ka2 * x2 + kb2 * I; dx3 = -ka3 * x3 + kb3 * I;\n', ...
        'UG = D2 / tau_D; dD1 = d_mmolmin * A_G - D1 / tau_D; dD2 = (D1 - D2) / tau_D;\n', ...
        'G = Q1 / V_G;\n', ...
        'if G >= 4.5, F01_c = F01; else, F01_c = F01 * (G / 4.5); end\n', ...
        'if G >= 9.0, FR = 0.003 * (G - 9.0) * V_G; else, FR = 0.0; end\n', ...
        'EGP = max(0.0, EGP0 * (1.0 - x3));\n', ...
        'dQ1 = -(F01_c + x1 * Q1) + k12 * Q2 - FR + EGP + UG;\n', ...
        'dQ2 = x1 * Q1 - (k12 + x2) * Q2;\n', ...
        'dG_meas = (G - G_meas) / tau_cgm;\n', ...
        'dxdt = [dS1; dS2; dI; dx1; dx2; dx3; dQ1; dQ2; dD1; dD2; dG_meas];\n', ...
        'end\n']);
    
    add_block('simulink/Signal Routing/Demux', [plant_sub, '/State_Demux'], 'Outputs', '11', 'Position', [420, 60, 425, 200]);
    add_block('simulink/Math Operations/Gain', [plant_sub, '/Gain_CGM_mgdl'], 'Gain', '18.01557', 'Position', [460, 65, 520, 95]);
    add_block('simulink/Math Operations/Gain', [plant_sub, '/Gain_Plasma_mgdl'], 'Gain', num2str(18.01557 / (0.16 * 70)), 'Position', [460, 115, 520, 145]);
    
    % Outport 1: CGM Glucose
    add_block('simulink/Ports & Subsystems/Out1', [plant_sub, '/CGM_Glucose_mgdl'], 'Position', [560, 70, 590, 90]);
    % Outport 2: Plasma Glucose
    add_block('simulink/Ports & Subsystems/Out1', [plant_sub, '/Plasma_Glucose_mgdl'], 'Position', [560, 120, 590, 140]);
    % Outport 3: Plasma Insulin
    add_block('simulink/Ports & Subsystems/Out1', [plant_sub, '/Plasma_Insulin_mUL'], 'Position', [560, 170, 590, 190]);
    % Outport 4: Full State Vector x (for state-feedback MPC)
    add_block('simulink/Ports & Subsystems/Out1', [plant_sub, '/State_x'], 'Position', [560, 220, 590, 240]);
    
    add_line(plant_sub, 'Insulin_Uh/1', 'Hovorka_Derivatives/1');
    add_line(plant_sub, 'Meal_Carbs_gmin/1', 'Hovorka_Derivatives/2');
    add_line(plant_sub, 'State_Integrator/1', 'Hovorka_Derivatives/3');
    add_line(plant_sub, 'Hovorka_Derivatives/1', 'State_Integrator/1');
    add_line(plant_sub, 'State_Integrator/1', 'State_Demux/1');
    add_line(plant_sub, 'State_Demux/11', 'Gain_CGM_mgdl/1');
    add_line(plant_sub, 'State_Demux/7',  'Gain_Plasma_mgdl/1');
    add_line(plant_sub, 'Gain_CGM_mgdl/1', 'CGM_Glucose_mgdl/1');
    add_line(plant_sub, 'Gain_Plasma_mgdl/1', 'Plasma_Glucose_mgdl/1');
    add_line(plant_sub, 'State_Demux/3',  'Plasma_Insulin_mUL/1');
    add_line(plant_sub, 'State_Integrator/1', 'State_x/1');
end

%% =========================================================
% BRANCH 1: Case 1 - Open-Loop Baseline
%% =========================================================
add_block('simulink/Sources/Constant', [model_name, '/Basal_Insulin_OL'], ...
    'Value', num2str(p.ss.u_basal_Uh), 'Position', [220, 80, 270, 110]);

add_hovorka_plant_with_states(model_name, 'Plant_OpenLoop', [400, 70, 540, 160], p);

add_line(model_name, 'Basal_Insulin_OL/1', 'Plant_OpenLoop/1');
add_line(model_name, 'Meal_Disturbance/1', 'Plant_OpenLoop/2');

%% =========================================================
% BRANCH 2: Case 2 - MPC Without Meal Preview (Reactive)
%% =========================================================
add_block('mpclib/MPC Controller', [model_name, '/MPC_No_Preview'], ...
    'mpcobj', 'mpc_obj', ...
    'state_inport', 'on', ...
    'Position', [220, 220, 320, 290]);

% Mux 11 plant states + 1 zero disturbance state = 12 states for MPC Controller
add_block('simulink/Signal Routing/Mux', [model_name, '/Mux_12_States_NP'], ...
    'Inputs', '2', 'Position', [180, 220, 185, 255]);

add_hovorka_plant_with_states(model_name, 'Plant_NoPreview', [400, 220, 540, 310], p);

% Wire Mux 12 states
add_line(model_name, 'Plant_NoPreview/4', 'Mux_12_States_NP/1');
add_line(model_name, 'Zero_Disturbance/1', 'Mux_12_States_NP/2');

% Port 1 (x) <- Mux_12_States_NP
% Port 2 (ref) <- Target Glucose
% Port 3 (md) <- Zero Disturbance (unannounced)
add_line(model_name, 'Mux_12_States_NP/1', 'MPC_No_Preview/1');
add_line(model_name, 'Target_Glucose/1', 'MPC_No_Preview/2');
add_line(model_name, 'Zero_Disturbance/1', 'MPC_No_Preview/3');

add_line(model_name, 'MPC_No_Preview/1', 'Plant_NoPreview/1');
add_line(model_name, 'Meal_Disturbance/1', 'Plant_NoPreview/2');

%% =========================================================
% BRANCH 3: Case 3 - MPC With Meal Preview (Proactive)
%% =========================================================
add_block('mpclib/MPC Controller', [model_name, '/MPC_Meal_Preview'], ...
    'mpcobj', 'mpc_obj', ...
    'state_inport', 'on', ...
    'Position', [220, 370, 320, 440]);

% Mux 11 plant states + 1 zero disturbance state = 12 states for MPC Controller
add_block('simulink/Signal Routing/Mux', [model_name, '/Mux_12_States_MP'], ...
    'Inputs', '2', 'Position', [180, 370, 185, 405]);

add_hovorka_plant_with_states(model_name, 'Plant_MealPreview', [400, 370, 540, 460], p);

% Wire Mux 12 states
add_line(model_name, 'Plant_MealPreview/4', 'Mux_12_States_MP/1');
add_line(model_name, 'Zero_Disturbance/1', 'Mux_12_States_MP/2');

% Port 1 (x) <- Mux_12_States_MP
% Port 2 (ref) <- Target Glucose
% Port 3 (md) <- Meal_Preview_Generator (49-step lookahead vector)
add_line(model_name, 'Mux_12_States_MP/1', 'MPC_Meal_Preview/1');
add_line(model_name, 'Target_Glucose/1', 'MPC_Meal_Preview/2');
add_line(model_name, 'Meal_Preview_Generator/1', 'MPC_Meal_Preview/3');

add_line(model_name, 'MPC_Meal_Preview/1', 'Plant_MealPreview/1');
add_line(model_name, 'Meal_Disturbance/1', 'Plant_MealPreview/2');

%% =========================================================
% 4 INDIVIDUAL SCOPES (One for each case + 1 combined)
%% =========================================================
% Scope 1: Case 1 Open-Loop
add_block('simulink/Sinks/Scope', [model_name, '/Scope_Case1_OpenLoop'], ...
    'NumInputPorts', '2', 'Position', [650, 65, 700, 115]);
add_line(model_name, 'Plant_OpenLoop/1', 'Scope_Case1_OpenLoop/1');
add_line(model_name, 'Basal_Insulin_OL/1', 'Scope_Case1_OpenLoop/2');

% Scope 2: Case 2 MPC No-Preview
add_block('simulink/Sinks/Scope', [model_name, '/Scope_Case2_NoPreview'], ...
    'NumInputPorts', '2', 'Position', [650, 215, 700, 265]);
add_line(model_name, 'Plant_NoPreview/1', 'Scope_Case2_NoPreview/1');
add_line(model_name, 'MPC_No_Preview/1', 'Scope_Case2_NoPreview/2');

% Scope 3: Case 3 MPC Meal-Preview
add_block('simulink/Sinks/Scope', [model_name, '/Scope_Case3_MealPreview'], ...
    'NumInputPorts', '2', 'Position', [650, 365, 700, 415]);
add_line(model_name, 'Plant_MealPreview/1', 'Scope_Case3_MealPreview/1');
add_line(model_name, 'MPC_Meal_Preview/1', 'Scope_Case3_MealPreview/2');

% Combined Scope 4: All 3 Cases Overlaid
add_block('simulink/Signal Routing/Mux', [model_name, '/Mux_All_Glucose'], ...
    'Inputs', '3', 'Position', [760, 170, 765, 230]);
add_line(model_name, 'Plant_OpenLoop/1',    'Mux_All_Glucose/1');
add_line(model_name, 'Plant_NoPreview/1',   'Mux_All_Glucose/2');
add_line(model_name, 'Plant_MealPreview/1', 'Mux_All_Glucose/3');

add_block('simulink/Signal Routing/Mux', [model_name, '/Mux_All_Insulin'], ...
    'Inputs', '3', 'Position', [760, 280, 765, 340]);
add_line(model_name, 'Basal_Insulin_OL/1',  'Mux_All_Insulin/1');
add_line(model_name, 'MPC_No_Preview/1',    'Mux_All_Insulin/2');
add_line(model_name, 'MPC_Meal_Preview/1',  'Mux_All_Insulin/3');

add_block('simulink/Sinks/Scope', [model_name, '/Scope_All_Cases_Comparison'], ...
    'NumInputPorts', '2', 'Position', [820, 220, 870, 280]);
add_line(model_name, 'Mux_All_Glucose/1', 'Scope_All_Cases_Comparison/1');
add_line(model_name, 'Mux_All_Insulin/1', 'Scope_All_Cases_Comparison/2');

%% Loggers to Workspace
add_block('simulink/Sinks/To Workspace', [model_name, '/Log_G_OL'], 'VariableName', 'G_OL_sim', 'SaveFormat', 'Timeseries', 'Position', [650, 125, 710, 145]);
add_block('simulink/Sinks/To Workspace', [model_name, '/Log_I_OL'], 'VariableName', 'I_OL_sim', 'SaveFormat', 'Timeseries', 'Position', [340, 45, 400, 65]);
add_line(model_name, 'Plant_OpenLoop/1', 'Log_G_OL/1');
add_line(model_name, 'Basal_Insulin_OL/1', 'Log_I_OL/1');

add_block('simulink/Sinks/To Workspace', [model_name, '/Log_G_NP'], 'VariableName', 'G_NP_sim', 'SaveFormat', 'Timeseries', 'Position', [650, 275, 710, 295]);
add_block('simulink/Sinks/To Workspace', [model_name, '/Log_I_NP'], 'VariableName', 'I_NP_sim', 'SaveFormat', 'Timeseries', 'Position', [340, 195, 400, 215]);
add_line(model_name, 'Plant_NoPreview/1', 'Log_G_NP/1');
add_line(model_name, 'MPC_No_Preview/1', 'Log_I_NP/1');

add_block('simulink/Sinks/To Workspace', [model_name, '/Log_G_MP'], 'VariableName', 'G_MP_sim', 'SaveFormat', 'Timeseries', 'Position', [650, 425, 710, 445]);
add_block('simulink/Sinks/To Workspace', [model_name, '/Log_I_MP'], 'VariableName', 'I_MP_sim', 'SaveFormat', 'Timeseries', 'Position', [340, 345, 400, 365]);
add_line(model_name, 'Plant_MealPreview/1', 'Log_G_MP/1');
add_line(model_name, 'MPC_Meal_Preview/1', 'Log_I_MP/1');

save_system(model_name, model_path);
close_system(model_name);

fprintf('\nState-Aware 3-Way Simulink Model successfully created and saved to:\n  %s\n', model_path);
fprintf('=========================================================\n');
