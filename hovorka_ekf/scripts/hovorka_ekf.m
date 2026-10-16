function [ekf, innov, K] = hovorka_ekf(action, ekf, varargin)
% HOVORKA_EKF Continuous-Discrete Extended Kalman Filter for the 11-State Hovorka T1D Model.
% Estimates all 11 unmeasured physiological states (including subcutaneous insulin depots
% S1, S2, plasma insulin I, and remote actions x1, x2, x3) from a single scalar noisy
% Continuous Glucose Monitor (CGM) reading y(k) [mg/dL].
%
% Usage:
%   ekf = hovorka_ekf('init', [], p, Ts, Q, R, x0, P0);
%   ekf = hovorka_ekf('predict', ekf, u_mumin, d_mmolmin, p);
%   [ekf, innov, K] = hovorka_ekf('update', ekf, y_meas_mgdl, p);

innov = [];
K = [];

switch lower(action)
    case 'init'
        p  = varargin{1};
        Ts = varargin{2};
        
        if length(varargin) >= 3 && ~isempty(varargin{3})
            Q = varargin{3};
        else
            % Default process noise covariance Q (tuned for metabolic variance)
            % States: [S1, S2, I, x1, x2, x3, Q1, Q2, D1, D2, G_meas]
            Q = diag([ ...
                1.0,    ... % S1: Subcutaneous depot 1 (mU^2)
                1.0,    ... % S2: Subcutaneous depot 2 (mU^2)
                0.25,   ... % I:  Plasma insulin (mU/L)^2
                1e-6,   ... % x1: Remote action 1
                1e-6,   ... % x2: Remote action 2
                1e-4,   ... % x3: Remote action 3 (EGP)
                4.0,    ... % Q1: Glucose compartment 1 (mmol^2)
                4.0,    ... % Q2: Glucose compartment 2 (mmol^2)
                1.0,    ... % D1: Gut compartment 1 (mmol^2)
                1.0,    ... % D2: Gut compartment 2 (mmol^2)
                0.05    ... % G_meas: Subcutaneous sensor glucose (mmol/L)^2
            ]);
        end
        
        if length(varargin) >= 4 && ~isempty(varargin{4})
            R = varargin{4};
        else
            % Default CGM measurement noise covariance R: sigma = 5.0 mg/dL -> R = 25 (mg/dL)^2
            R = 25.0; 
        end
        
        if length(varargin) >= 5 && ~isempty(varargin{5})
            x0 = varargin{5};
        else
            x0 = p.ss.x0; % Nominal steady state
        end
        
        if length(varargin) >= 6 && ~isempty(varargin{6})
            P0 = varargin{6};
        else
            % Initial estimation covariance
            P0 = diag([10, 10, 5, 1e-4, 1e-4, 1e-3, 20, 20, 10, 10, 1.0]);
        end
        
        ekf = struct();
        ekf.x_hat = x0(:);
        ekf.P     = P0;
        ekf.Q     = Q;
        ekf.R     = R;
        ekf.Ts    = Ts;
        ekf.nx    = 11;
        ekf.H     = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, p.MMOLL_TO_MGDL]; % Output matrix

    case 'predict'
        % ekf = hovorka_ekf('predict', ekf, u_mumin, d_mmolmin, p)
        u_mumin   = varargin{1};
        d_mmolmin = varargin{2};
        p         = varargin{3};
        
        Ts = ekf.Ts;
        x_curr = ekf.x_hat;
        
        % 1. Nonlinear State Propagation over [0, Ts] using Runge-Kutta 4th Order
        f = @(x) hovorka_ode(0, x, u_mumin, d_mmolmin, p);
        
        dt_sub = Ts / 5; % 1-minute internal sub-steps for high numerical precision
        x_prop = x_curr;
        for s = 1:5
            k1 = f(x_prop);
            k2 = f(x_prop + 0.5 * dt_sub * k1);
            k3 = f(x_prop + 0.5 * dt_sub * k2);
            k4 = f(x_prop + dt_sub * k3);
            x_prop = x_prop + (dt_sub / 6) * (k1 + 2*k2 + 2*k3 + k4);
        end
        ekf.x_hat = max(0, x_prop); % Physical positivity constraint
        
        % 2. Analytical Jacobian Evaluation at current state estimate
        F = hovorka_jacobian(ekf.x_hat, p);
        
        % 3. Discrete State Transition Matrix Phi via 1st/2nd-order matrix exponential
        % Phi = expm(F * Ts) approx I + F*Ts + 0.5*(F*Ts)^2
        I_nx = eye(ekf.nx);
        Phi = expm(F * Ts);
        
        % 4. Covariance Propagation
        % Qd = integral of Phi(tau)*Q*Phi(tau)' dtau approx Q * Ts
        Qd = ekf.Q * Ts;
        ekf.P = Phi * ekf.P * Phi' + Qd;
        
        % Ensure numerical symmetry of covariance matrix
        ekf.P = 0.5 * (ekf.P + ekf.P');

    case 'update'
        % [ekf, innov, K] = hovorka_ekf('update', ekf, y_meas_mgdl, p)
        y_meas_mgdl = varargin{1};
        p           = varargin{2};
        
        if isnan(y_meas_mgdl) || isempty(y_meas_mgdl)
            % SENSOR DROPOUT / PACKET LOSS FAIL-SAFE:
            % When no valid CGM reading arrives, skip measurement update.
            % The state estimate and covariance remain at their propagated prior values.
            % Covariance P remains expanded, reflecting growing uncertainty.
            innov = NaN;
            K = zeros(ekf.nx, 1);
            return;
        end
        
        H = ekf.H;
        x_prior = ekf.x_hat;
        P_prior = ekf.P;
        
        % 1. Innovation / Measurement Residual
        y_hat = H * x_prior;
        innov = y_meas_mgdl - y_hat;
        
        % 2. Innovation Covariance
        S = H * P_prior * H' + ekf.R;
        
        % 3. Kalman Gain
        K = (P_prior * H') / S;
        
        % 4. State Update
        x_post = x_prior + K * innov;
        ekf.x_hat = max(0, x_post); % Physiological bounds
        
        % 5. Joseph Form Covariance Update (Guarantees Positive Definiteness)
        I_KH = eye(ekf.nx) - K * H;
        ekf.P = I_KH * P_prior * I_KH' + K * ekf.R * K';
        ekf.P = 0.5 * (ekf.P + ekf.P');
        
    otherwise
        error('Unknown EKF action: %s. Use ''init'', ''predict'', or ''update''.', action);
end
end

%% Analytical Jacobian Function for Hovorka Model
function F = hovorka_jacobian(x, p)
% Evaluates continuous-time Jacobian F = df/dx analytically
F = zeros(11, 11);

S1 = x(1); S2 = x(2); I = x(3);
x1 = x(4); x2 = x(5); x3 = x(6);
Q1 = x(7); Q2 = x(8); D1 = x(9); D2 = x(10); G_meas = x(11);

% Subsystem 1: Subcutaneous Insulin (dS1, dS2)
F(1, 1) = -1 / p.tau_S;
F(2, 1) =  1 / p.tau_S;
F(2, 2) = -1 / p.tau_S;

% Subsystem 2: Plasma Insulin (dI)
F(3, 2) =  1 / (p.tau_S * p.V_I);
F(3, 3) = -p.ke;

% Subsystem 3: Insulin Action (dx1, dx2, dx3)
F(4, 3) =  p.kb1;  F(4, 4) = -p.ka1;
F(5, 3) =  p.kb2;  F(5, 5) = -p.ka2;
F(6, 3) =  p.kb3;  F(6, 6) = -p.ka3;

% Subsystem 4: Gut Absorption (dD1, dD2)
F(9, 9)   = -1 / p.tau_D;
F(10, 9)  =  1 / p.tau_D;
F(10, 10) = -1 / p.tau_D;

% Subsystem 5: Glucose Compartments (dQ1, dQ2)
G = Q1 / p.V_G;

% dF01_c / dQ1
if G >= 4.5
    dF01_dQ1 = 0;
else
    dF01_dQ1 = p.F01 / (4.5 * p.V_G);
end

% dFR / dQ1
if G >= 9.0
    dFR_dQ1 = 0.003;
else
    dFR_dQ1 = 0;
end

% Row 7: dQ1
F(7, 4)  = -Q1;                     % d(dQ1)/dx1
F(7, 6)  = -p.EGP0;                 % d(dQ1)/dx3
F(7, 7)  = -(dF01_dQ1 + x1 + dFR_dQ1); % d(dQ1)/dQ1
F(7, 8)  =  p.k12;                  % d(dQ1)/dQ2
F(7, 10) =  1 / p.tau_D;            % d(dQ1)/dD2 (UG)

% Row 8: dQ2
F(8, 4)  =  Q1;                     % d(dQ2)/dx1
F(8, 5)  = -Q2;                     % d(dQ2)/dx2
F(8, 7)  =  x1;                     % d(dQ2)/dQ1
F(8, 8)  = -(p.k12 + x2);           % d(dQ2)/dQ2

% Subsystem 6: CGM Delay (dG_meas)
F(11, 7)  =  1 / (p.V_G * p.tau_cgm);
F(11, 11) = -1 / p.tau_cgm;

end
