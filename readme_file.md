# VBD MODE Interview Exercise

Implementing a 4th-order Runge-Kutta numerical solver for a system of ordinary differential equations 

## Exercise Description

This python project implements algorithms 1 and 2 from the VBD MODE interview exercise. It solves a system of three coupled ODEs representing population dynamics over 365-day period
### algorithm 1: RK4 Function
Classic 4th-order RK method with four stages (k₁, k₂, k₃, k₄)

### algorithm 2: System Function
Computes derivative vector dy = [dE, dJ, dA] for the ODE system

## Requirements

- Python 3.7+
- NumPy
- Matplotlib

## Usage

My code aims to 
1. Solve the ODE system using the RK4 method
2. Display progress during integration
3. Print final equilibrium values
4. Generate, display, and save a plot (`results.png`)


## System Description

The system consists of three state variables (E, J, A) governed by:

- **dE/dt** = βA - δ_E·E - μ_E·E
- **dJ/dt** = δ_E·E - δ_J·J - α·J² - μ_J·J
- **dA/dt** = ω·δ_J·J - μ_A·A

### Parameters
- β = 24
- δ_E = 0.6, μ_E = 0.15
- δ_J = 0.08, μ_J = 0.05
- α = 0.003
- ω = 0.5
- μ_A = 0.1

### Initial Conditions
- y₀ = [10, 0, 0]
- Time span: [0, 365] days
- Time step: h = 0.01

## Output description

The script produces:
1. **Console output**: Progress updates and final eqbrm values
2. **Plot window**: Three subplots for E(t), J(t), and A(t)
3. **Image file**: `results.png` 

### Typical Results
The system reaches equilibrium around t ≈ 150-200 days with approximate values:
- E ≈ 32,000
- J ≈ 2,500
- A ≈ 1,000


