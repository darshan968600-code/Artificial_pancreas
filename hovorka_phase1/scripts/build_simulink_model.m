% BUILD_SIMULINK_MODEL Programmatically constructs the Simulink model
% 'models/hovorka_artificial_pancreas.slx' for the Hovorka Artificial Pancreas.
% Safe State-Feedback Version: Full IOB tracking, zero algebraic loops.

clc; clear; close all;
fprintf('=========================================================\n');
fprintf('  BUILDING SAFE STATE-FEEDBACK SIMULINK MODEL (.slx)\n');
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

% 1. Prepare parameters and MPC controller in base workspace
p = hovorka_params(70, 100);
[mpc_obj, plant_ss_d] = design_mpc_controller();
assignin('base', 'p', p);
assignin('base', 'mpc_obj', mpc_obj);
assignin('base', 'plant_ss_d', plant_ss_d);

% Model name and file path
model_name = 'hovorka_artificial_pancreas';
model_dir  = fullfile(project_root, 'models');
if ~exist(model_dir, 'dir'), mkdir(model_dir); end
model_path = fullfile(model_dir, [model_name, '.slx']);

% Close if already open
if bdIsLoaded(model_name)
    close_system(model_name, 0);
end

% 2. Create new Simulink system
new_system(model_name);
open_system(model_name);

set_param(model_name, 'Solver', 'ode15s');
set_param(model_name, 'StopTime', '720');
set_param(model_name, 'RelTol', '1e-6');
set_param(model_name, 'AbsTol', '1e-8');
set_param(model_name, 'SaveOutput', 'on');
set_param(model_name, 'SignalLogging', 'on');
set_param(model_name, 'SignalLoggingName', 'logsout');

init_cmd = 'addpath(fullfile(fileparts(which(''hovorka_artificial_pancreas'')), ''..'', ''scripts'')); p = hovorka_params(70, 100); [mpc_obj, plant_ss_d] = design_mpc_controller();';
set_param(model_name, 'InitFcn', init_cmd);

%% 3. Add Reference Block (Target = 100 mg/dL)
add_block('simulink/Sources/Constant', [model_name, '/Target_Glucose'], ...
    'Value', '100', 'Position', [50, 100, 100, 130]);

%% 4. Add Meal Disturbance Generator (Pulse: 50g at t = 120 min over 15 min)
add_block('simulink/Sources/Step', [model_name, '/Meal_Start'], ...
    'Time', '120', 'Before', '0', 'After', '3.3333', 'Position', [50, 250, 100, 280]);
add_block('simulink/Sources/Step', [model_name, '/Meal_End'], ...
    'Time', '135', 'Before', '0', 'After', '3.3333', 'Position', [50, 310, 100, 340]);
add_block('simulink/Math Operations/Subtract', [model_name, '/Meal_Pulse'], ...
    'Position', [140, 270, 170, 300]);

add_block('simulink/Sources/Constant', [model_name, '/Zero_Dist_State'], ...
    'Value', '0', 'Position', [140, 330, 170, 350]);

% Meal Preview Generator (49-step lookahead vector over Np = 48 prediction horizon)
add_block('simulink/Sources/Clock', [model_name, '/Sim_Clock'], ...
    'Position', [50, 360, 80, 380]);

preview_fcn = [model_name, '/Meal_Preview_Generator'];
add_block('simulink/User-Defined Functions/MATLAB Function', preview_fcn, ...
    'Position', [100, 350, 170, 390]);
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

%% 5. Add MPC Controller Block with State Feedback
add_block('mpclib/MPC Controller', [model_name, '/MPC_Controller'], ...
    'mpcobj', 'mpc_obj', ...
    'state_inport', 'on', ...
    'Position', [250, 140, 370, 220]);

% Mux 11 plant states + 1 disturbance state = 12 states for MPC
add_block('simulink/Signal Routing/Mux', [model_name, '/Mux_12_States'], ...
    'Inputs', '2', 'Position', [210, 135, 215, 175]);

%% 6. Add Plant Subsystem (Hovorka Non-Linear Differential Equations)
plant_sub = [model_name, '/Hovorka_Plant'];
add_block('simulink/Ports & Subsystems/Subsystem', plant_sub, ...
    'Position', [480, 150, 660, 260]);

delete_line(plant_sub, 'In1/1', 'Out1/1');
delete_block([plant_sub, '/In1']);
delete_block([plant_sub, '/Out1']);

add_block('simulink/Ports & Subsystems/In1', [plant_sub, '/Insulin_Uh'], 'Position', [30, 80, 60, 100]);
add_block('simulink/Ports & Subsystems/In1', [plant_sub, '/Meal_Carbs_gmin'], 'Position', [30, 160, 60, 180]);

add_block('simulink/Continuous/Integrator', [plant_sub, '/State_Integrator'], ...
    'InitialCondition', 'p.ss.x0', 'Position', [340, 105, 380, 155]);

matlab_fcn = [plant_sub, '/Hovorka_Derivatives'];
add_block('simulink/User-Defined Functions/MATLAB Function', matlab_fcn, ...
    'Position', [130, 75, 270, 185]);

sf = sfroot;
m_obj = sf.find('Path', matlab_fcn, '-isa', 'Stateflow.EMChart');
fcn_code = sprintf([ ...
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
m_obj.Script = fcn_code;

add_block('simulink/Signal Routing/Demux', [plant_sub, '/State_Demux'], ...
    'Outputs', '11', 'Position', [420, 60, 425, 200]);

add_block('simulink/Math Operations/Gain', [plant_sub, '/Gain_CGM_mgdl'], ...
    'Gain', '18.01557', 'Position', [460, 65, 520, 95]);
add_block('simulink/Math Operations/Gain', [plant_sub, '/Gain_Plasma_mgdl'], ...
    'Gain', num2str(18.01557 / (0.16 * 70)), 'Position', [460, 115, 520, 145]);

add_block('simulink/Ports & Subsystems/Out1', [plant_sub, '/CGM_Glucose_mgdl'], 'Position', [560, 70, 590, 90]);
add_block('simulink/Ports & Subsystems/Out1', [plant_sub, '/Plasma_Glucose_mgdl'], 'Position', [560, 120, 590, 140]);
add_block('simulink/Ports & Subsystems/Out1', [plant_sub, '/Plasma_Insulin_mUL'], 'Position', [560, 170, 590, 190]);
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

%% 7. Add Scopes and To Workspace Logging
add_block('simulink/Sinks/Scope', [model_name, '/Glucose_Scope'], ...
    'Position', [730, 140, 770, 180]);
add_block('simulink/Sinks/Scope', [model_name, '/Insulin_Scope'], ...
    'Position', [420, 60, 460, 100]);

add_block('simulink/Sinks/To Workspace', [model_name, '/Log_Glucose'], ...
    'VariableName', 'G_sim', 'SaveFormat', 'Timeseries', 'Position', [730, 200, 790, 230]);
add_block('simulink/Sinks/To Workspace', [model_name, '/Log_Insulin'], ...
    'VariableName', 'I_sim', 'SaveFormat', 'Timeseries', 'Position', [420, 20, 480, 50]);

%% 8. Wire Top-Level Diagram
add_line(model_name, 'Meal_Start/1', 'Meal_Pulse/1');
add_line(model_name, 'Meal_End/1', 'Meal_Pulse/2');
add_line(model_name, 'Sim_Clock/1', 'Meal_Preview_Generator/1');

% 12-state vector: 11 plant states + 1 zero disturbance state
add_line(model_name, 'Hovorka_Plant/4', 'Mux_12_States/1');
add_line(model_name, 'Zero_Dist_State/1', 'Mux_12_States/2');

% MPC Controller connections
add_line(model_name, 'Mux_12_States/1', 'MPC_Controller/1');
add_line(model_name, 'Target_Glucose/1', 'MPC_Controller/2');
add_line(model_name, 'Meal_Preview_Generator/1', 'MPC_Controller/3');

% Plant connections
add_line(model_name, 'MPC_Controller/1', 'Hovorka_Plant/1');
add_line(model_name, 'Meal_Pulse/1', 'Hovorka_Plant/2');

% Scopes & Loggers
add_line(model_name, 'Hovorka_Plant/1', 'Glucose_Scope/1');
add_line(model_name, 'Hovorka_Plant/1', 'Log_Glucose/1');
add_line(model_name, 'MPC_Controller/1', 'Insulin_Scope/1');
add_line(model_name, 'MPC_Controller/1', 'Log_Insulin/1');

%% 9. Save and Close System
save_system(model_name, model_path);
close_system(model_name);

fprintf('\nSafe State-Feedback Simulink Model successfully created and saved to:\n  %s\n', model_path);
fprintf('=========================================================\n');
