function fix_model_callbacks()
%% FIX_MODEL_CALLBACKS Automatically configures PreLoadFcn and InitFcn on all
% Bergman Simulink models so that they load their required MPC controllers automatically
% without requiring the user to manually load variables into the MATLAB workspace.

clc; close all;
fprintf('========================================================================\n');
fprintf('  CONFIGURING AUTOMATIC MPC LOADING CALLBACKS ON BERGMAN MODELS\n');
fprintf('========================================================================\n');

this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(this_dir); % C:\Pancreas_project\bergman_model
model_dir = fullfile(project_root, 'models');

models = { ...
    'Artificial_Pancreas_OpenLoop_Comparison', ...
    'Artificial_Pancreas_MPC_NoPreview', ...
    'Artificial_Pancreas_MPC_MealPreview', ...
    'Artificial_Pancreas_OpenLoop_LargeMeal', ...
    'Artificial_Pancreas_MPC_NoPreview_LargeMeal', ...
    'Artificial_Pancreas_MPC_MealPreview_LargeMeal'};

% Robust callback code to load both controllers whenever model opens or simulates
init_code = sprintf([ ...
    '%% Automatic Controller Loader Callback\n', ...
    'mdl_name = bdroot;\n', ...
    'mdl_file = which(mdl_name);\n', ...
    'if ~isempty(mdl_file)\n', ...
    '    m_dir = fileparts(mdl_file);\n', ...
    'else\n', ...
    '    m_dir = pwd;\n', ...
    'end\n', ...
    'b_file = fullfile(m_dir, ''mpc_baseline.mat'');\n', ...
    'if exist(b_file, ''file'')\n', ...
    '    s_base = load(b_file, ''mpcobj'');\n', ...
    '    assignin(''base'', ''mpcobj'', s_base.mpcobj);\n', ...
    'end\n', ...
    'p_file = fullfile(m_dir, ''mpc_meal_preview.mat'');\n', ...
    'if exist(p_file, ''file'')\n', ...
    '    s_prev = load(p_file, ''mpcobj_md'');\n', ...
    '    assignin(''base'', ''mpcobj_md'', s_prev.mpcobj_md);\n', ...
    'end\n']);

for k = 1:length(models)
    m = models{k};
    m_path = fullfile(model_dir, [m, '.slx']);
    fprintf('Configuring callbacks for %s...\n', m);
    
    if bdIsLoaded(m)
        close_system(m, 0);
    end
    
    load_system(m_path);
    set_param(m, 'PreLoadFcn', init_code);
    set_param(m, 'InitFcn',    init_code);
    save_system(m, m_path);
    close_system(m);
    fprintf('  [OK] Callbacks configured and saved.\n');
end

fprintf('\n========================================================================\n');
fprintf('VERIFYING STANDALONE SIMULATION OF EACH MODEL (WITH EMPTY WORKSPACE):\n');
fprintf('========================================================================\n');

for k = 1:length(models)
    m = models{k};
    m_path = fullfile(model_dir, [m, '.slx']);
    
    % Clear base workspace completely
    evalin('base', 'clear all');
    
    try
        out = sim(m_path);
        G_val = out.logsout.getElement(1).Values.Data;
        I_val = out.logsout.getElement(2).Values.Data;
        fprintf('[PASS] %-45s Peak G = %6.2f mg/dL, Max I = %5.2f uU/mL\n', ...
            m, max(G_val), max(I_val));
    catch ME
        fprintf('[FAIL] %-45s ERROR: %s\n', m, ME.message);
    end
end

fprintf('========================================================================\n');
end
