# PRDC Pricing via Multi-Dimensional PDE

## Overview
This repository implements a full quantitative pricing pipeline for **Power Reverse Dual Currency (PRDC)** derivatives in the FX market. It consists of two independent modules:
- **Local Volatility Calibration** (Python) — calibrates the FX local volatility surface from market data.
- **Multi-Dimensional PDE Pricer** (C++) — prices FX derivatives via an ADI finite difference scheme on a 3D grid.

---

## Repository Structure

```text
prdc-pricing-multidim-pde/
├── local-vol-calibration/        # Python — Local volatility calibration
│   ├── lib/
│   │   ├── bsanalytic.py
│   │   ├── calibrator.py
│   │   ├── fxivolinterpolator.py
│   │   ├── surfaces.py
│   │   └── interpolator.py
│   ├── sim_lib/
│   │   └── StochasticSim_Multiprocessing.py
│   ├── marketdata_JSON_asof_04_30_2020/
│   └── FX_LocalVol_Calibration.ipynb
│
└── pde-pricer/                   # C++ — 3D ADI PDE pricer
    ├── grid3d.hpp
    ├── model_params.hpp
    ├── fd_coeffs.hpp
    ├── apply_operators.hpp
    ├── thomas_solver.hpp
    ├── implicit_solvers.hpp
    ├── boundary.hpp
    ├── initial_condition.hpp
    ├── douglas_step.hpp
    ├── interpolation.hpp
    ├── benchmark.hpp
    └── main.cpp
```

---

## Module 1 — Local Volatility Calibration (Python)

### Model & Calibration
The FX spot process follows a **Local Volatility model with 2 Stochastic Rates (LV2SR)** modeled as G1++ (Hull-White) processes. Calibration is performed slice by slice in the T-Forward measure using Monte Carlo simulations and the extended Dupire formula.

### Stack
- `numpy`: Numerical simulation
- `scipy`: Interpolation & optimization
- `pandas`: Market data handling

### How to run
```bash
cd local-vol-calibration
pip install -r requirements.txt
jupyter notebook FX_LocalVol_Calibration.ipynb
```

---

## Module 2 — Multi-Dimensional PDE Pricer (C++)

### Model & Numerical Scheme
The pricer solves the 3-factor Black-Scholes PDE on a 3D grid (Spot, Domestic Rate, Foreign Rate) using the **Douglas ADI finite difference scheme** with standard boundary conditions for FX Call/Put options. 

### Validation
Results are successfully benchmarked against the Garman-Kohlhagen analytical formula, confirming convergence accuracy.



### Build & Run
```bash
cd pde-pricer
mkdir build && cd build
cmake ..
make
./prdc_pricer
```

---
