% TEST_TUNING Tests various MPC tuning parameters to find the optimal
% clinical setup that prevents post-meal hypoglycemia.

clc; clear; close all;
this_dir = fileparts(mfilename('fullpath'));
if endsWith(this_dir, 'scripts'), project_root = fileparts(this_dir); else, project_root = this_dir; end
addpath(fullfile(project_root, 'scripts'));

p = hovorka_params(70, 100);

% Simulation parameters
Ts = 5;
t_end = 720;
N_steps = t_end / Ts;
time_vec = (0:N_steps)' * Ts;
meal_start_time = 120;
meal_dur = 15;
meal_carbs_g = 50;
meal_rate_gmin = meal_carbs_g / meal_dur;

d_meal_gmin = zeros(N_steps + 1, 1);
for k = 1:(N_steps + 1)
    if time_vec(k) >= meal_start_time && time_vec(k) < (meal_start_time + meal_dur)
        d_meal_gmin(k) = meal_rate_gmin;
    end
end

ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Test grid of parameters
u_max_list = [1.2, 1.5, 2.0];
w_u_list   = [5.0, 10.0, 20.0];

fprintf('Evaluating tuning combinations on Hovorka plant...\n');
fprintf('%-10s %-10s %-15s %-15s %-15s\n', 'u_max', 'W_u', 'Peak Glucose', 'Min Glucose', 'Total Insulin');
fprintf('-----------------------------------------------------------------------\n');

for u_max = u_max_list
    for w_u = w_u_list
        [mpc_obj, ~] = design_mpc_controller();
        mpc_obj.ManipulatedVariables(1).Max = u_max;
        mpc_obj.Weights.ManipulatedVariables = w_u;
        mpc_obj.Weights.ManipulatedVariablesRate = 20.0;
        
        % Run simulation with meal preview
        x_sim = zeros(N_steps + 1, 11);
        u_sim = zeros(N_steps, 1);
        x_sim(1, :) = p.ss.x0';
        xmpc = mpcstate(mpc_obj);
        
        for k = 1:N_steps
            t_span_k = [time_vec(k), time_vec(k+1)];
            y_meas_k = x_sim(k, 11) * p.MMOLL_TO_MGDL;
            
            v_preview = zeros(mpc_obj.PredictionHorizon + 1, 1);
            for h = 0:mpc_obj.PredictionHorizon
                idx_f = min(k + h, N_steps + 1);
                v_preview(h + 1) = d_meal_gmin(idx_f);
            end
            
            u_opt = mpcmove(mpc_obj, xmpc, y_meas_k, 100, v_preview);
            u_sim(k) = u_opt;
            
            u_k_mumin = u_opt * p.UH_TO_MUMIN;
            d_k_mmolmin = d_meal_gmin(k) * (1000 / p.MOL_WT_GLUCOSE);
            [~, x_seg] = ode15s(@(t, x) hovorka_ode(t, x, u_k_mumin, d_k_mmolmin, p), ...
                t_span_k, x_sim(k, :)', ode_opts);
            x_sim(k+1, :) = x_seg(end, :);
        end
        
        G_sim_mgdl = (x_sim(:, 7) ./ p.V_G) .* p.MMOLL_TO_MGDL;
        peak_g = max(G_sim_mgdl);
        min_g  = min(G_sim_mgdl);
        tot_u  = sum(u_sim) * (Ts / 60);
        
        fprintf('%-10.1f %-10.1f %-15.2f %-15.2f %-15.2f\n', ...
            u_max, w_u, peak_g, min_g, tot_u);
    end
end
