# Artificial Pancreas System: Hovorka & Bergman Models with MPC and EKF

This repository contains a comprehensive MATLAB & Simulink research codebase for an **Artificial Pancreas System** for Type 1 Diabetes (T1D). It models glucose-insulin regulatory dynamics, implements advanced control strategies (Model Predictive Control with and without meal preview), and provides state estimation using an Extended Kalman Filter (EKF) under realistic noisy CGM conditions.

---

## 🌟 Key Features

1. **Physiological Models**:
   - **Hovorka Model**: Nonlinear 8-state differential equation model capturing subcutaneous insulin absorption, plasma insulin kinetics, and three subsystem insulin actions (glucose distribution, utilization, and endogenous glucose production suppression).
   - **Bergman Minimal Model**: Classical 3-state glucose-insulin minimal model for comparative benchmarking.
2. **Control Strategies**:
   - **Linearized Model Predictive Control (MPC)**: Formulated state-space discrete-time MPC with input saturation constraints to prevent hypoglycemia.
   - **Meal Preview Feedforward MPC**: Anticipatory bolus/basal adjustment taking into account known meal carbohydrate disturbances.
   - **Open-Loop Basal Control**: Baseline insulin infusion comparisons.
3. **State Estimation**:
   - **Extended Kalman Filter (EKF)**: Continuous-discrete nonlinear state estimator for full unmeasured physiological state recovery under realistic Continuous Glucose Monitor (CGM) sensor noise and meal disturbances.
4. **Benchmarking & Validation**:
   - 24-hour simulation profiles featuring multi-meal challenges (breakfast, lunch, dinner).
   - Direct Hovorka vs. Bergman comparative trajectory analyses.
   - Simulink block diagrams and standalone MATLAB simulation scripts.

---

## 📂 Repository Structure

```text
Pancreas_project/
├── hovorka_phase1/                 # Hovorka model, open-loop, MPC controllers & 24h simulations
│   ├── models/                    # Simulink .slx models
│   ├── scripts/                   # ODE functions, MPC design, and test scripts
│   ├── results/                   # Figures and simulation .mat files
│   └── docs/                      # Model specifications and notes
├── bergman_model/                  # Bergman minimal model implementations and simulations
│   ├── models/                    # Simulink .slx models
│   ├── scripts/                   # MPC controllers and simulation runners
│   └── results/                   # Figures and simulation data
├── hovorka_bergman_comparision/    # Benchmark comparison between Hovorka and Bergman models
│   ├── scripts/                   # Comparative simulation scripts
│   └── results/                   # Side-by-side comparison plots
├── hovorka_ekf/                    # Extended Kalman Filter (EKF) state estimation
│   ├── models/                    # Simulink model with EKF
│   ├── scripts/                   # EKF continuous-discrete filter implementation
│   └── results/                   # State estimation error and trajectory figures
└── README.md                       # Project overview and instructions
```

---

## 🚀 Getting Started

### Prerequisites
- **MATLAB** (R2021a or newer recommended)
- **Simulink**
- **Control System Toolbox** (for state-space and MPC formulation)
- **Optimization Toolbox** (optional, for custom QP formulations)

### Running Simulations

1. **Clone the repository**:
   ```bash
   git clone https://github.com/<your-username>/<your-repo-name>.git
   cd <your-repo-name>
   ```

2. **Run Hovorka 24-Hour 3-Meal Simulation**:
   In MATLAB, navigate to `hovorka_phase1/scripts/test_24h_3meals` and run:
   ```matlab
   run_24h_simulation
   ```

3. **Run Hovorka vs. Bergman Benchmark**:
   In MATLAB, navigate to `hovorka_bergman_comparision/scripts` and run:
   ```matlab
   compare_hovorka_bergman
   ```

4. **Run Hovorka EKF Estimation**:
   In MATLAB, navigate to `hovorka_ekf/scripts` and run:
   ```matlab
   run_hovorka_ekf_simulation
   ```

---

## 📊 Summary of Results

- **Meal Preview MPC** significantly reduces postprandial glucose peaks compared to standard MPC and open-loop basal insulin, preventing glycemic excursions above 180 mg/dL.
- **EKF** provides fast convergence and state estimation across all 8 internal states of the Hovorka model despite noisy CGM readings and meal ingestion.

---

## 📄 License
This project is for research and educational purposes.
