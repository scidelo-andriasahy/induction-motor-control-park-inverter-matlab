# Three-Phase Induction Motor Modeling & Control — Park Reference Frames, AC Dimmers & V/f Inverter (MATLAB)

[![Portfolio](https://img.shields.io/badge/Portfolio-scidelo--andriasahy.github.io-0284c7?style=for-the-badge)](https://scidelo-andriasahy.github.io/?project=9)
[![MATLAB](https://img.shields.io/badge/Simulation-MATLAB_RK4_%26_FFT-e16737?style=for-the-badge)](https://www.mathworks.com/products/matlab.html)

> **Author:** Tsikomia Scidelo ANDRIASAHY — GSEA Engineering (ENSA Tangier)

![Induction Motor Park Frame Comparison](figures/fig_tp1_comparaison_reperes.png)

## Overview
Complete numerical modeling (RK4 integration) and power-electronic control study of a three-phase squirrel-cage induction machine (MAS):
- **Park Reference Frame Analysis:** Stator, rotor, and synchronous (`dq0`) dynamic simulations under no-load, rated load, and pull-out torque stall conditions.
- **AC Voltage Controllers (Gradateurs):** Phase-angle vs. burst-firing (wave-train) speed regulation and harmonic THD spectral analysis (`calc_fft_thd.m`).
- **PWM V/f Inverter Drive:** Constant flux `V/f` scalar control for variable-speed operation.
