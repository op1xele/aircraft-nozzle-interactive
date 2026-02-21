# Aircraft Nozzle Interactive
**Advanced Rocket Nozzle Design Suite**

A comprehensive desktop application for rocket nozzle design, analysis, and optimization using the Method of Characteristics (MOC).

---

## Overview

Aircraft Nozzle Interactive is a Python-based engineering tool that implements the Method of Characteristics for supersonic nozzle flow analysis. It provides a full design-to-analysis workflow including geometry generation, performance calculation, thermal analysis, structural analysis, weight estimation, and optimization — all within an interactive GUI.

---

## Supported Nozzle Types

| # | Type | Description |
|---|------|-------------|
| 1 | 2D Ideal | Two-dimensional ideal nozzle (MOC) |
| 2 | Axisymmetric | Axisymmetric nozzle (MOC) |
| 3 | Plug / Aerospike | C.C. Lee plug nozzle with internal/external expansion control |
| 4 | Cone | Conical nozzle |
| 5 | Wedge | Wedge nozzle |
| 6 | Bell | Bell-shaped nozzle (Rao approximation) |
| 7 | Truncated Ideal | Truncated ideal contour |
| 8 | Minimum Length | Minimum length nozzle |
| 9 | Throat Expansion | Throat expansion nozzle |

---

## Features

### Core Engine
- Method of Characteristics (MOC) solver for supersonic flow
- Isentropic flow relations with real gas effects (NASA 9-coefficient polynomials)
- Prandtl-Meyer expansion and shock wave calculations
- Plume analysis with shock reflection modeling
- Coalescence detection and reporting

### Engineering Modules
- **Viscous Correction** — Boundary layer displacement thickness (Bartz + Thwaites integral method, Mangler axisymmetric correction)
- **Thermal Analysis** — Wall temperature distribution via Bartz correlation, heat flux, heat transfer coefficient, regenerative cooling channel sizing
- **Structural Analysis** — Hoop stress, axial stress, thermal stress, Von Mises combined stress, safety factors, fatigue life (Goodman criterion)
- **Weight Estimation** — Shell mass integration by section, material comparison, thrust-to-weight ratio
- **Optimization** — Nelder-Mead simplex and genetic algorithm for target Isp/efficiency optimization

### Real Gas Effects
- Temperature-dependent gamma using NASA 9-coefficient thermodynamic polynomials
- Propellants: H2/O2, RP-1/LOX, CH4/LOX, N2O4/MMH, air
- Frozen vs. equilibrium flow comparison

### Output and Analysis Tools
- 138 tools accessible from the Tools menu
- Convergence monitor, grid independence study, error estimation
- Parametric study, sensitivity analysis, DOE setup
- Nozzle comparison, performance map generator
- Report generator (HTML, TXT, CSV, PDF)
- CFD export (OpenFOAM, Tecplot, VTK, CGNS, PLOT3D)
- DXF/STL/G-code geometry export

### Material Database
6 materials with full thermophysical properties:
- Inconel 718, Copper C18150, Titanium Ti-6Al-4V
- Carbon-Carbon, Stainless 316L, Rhenium

### Unit System
- Full SI / Imperial toggle (psia/kPa, R/K, in/cm, lbf/N, ft/s / m/s)

---

## Requirements

```
Python 3.10+
tkinter
numpy
matplotlib
scipy
reportlab
```

Install dependencies:

```
pip install numpy matplotlib scipy reportlab
```

---

## Running

```
python main.py
```

Or use the compiled executable:

```
Aircraft Nozzle Calculator.exe
```

---

## Project Structure

```
Aircraft Nozzle Interactive/
├── main.py
├── config.py
├── core/
│   ├── moc_engine.py
│   ├── isentropic.py
│   ├── performance.py
│   ├── geometry.py
│   ├── shock_calculator.py
│   ├── external_flow.py
│   ├── viscous.py
│   ├── thermal.py
│   ├── structural.py
│   ├── weight.py
│   ├── optimization.py
│   └── pdf_report.py
├── gui/
│   ├── main_window.py
│   ├── canvas.py
│   ├── input_panels.py
│   └── output_panels.py
└── utils/
    ├── helpers.py
    └── validators.py
```

---

## Method of Characteristics

The MOC solver implements the characteristic line method for supersonic flows:
- Right-running (C+) and left-running (C-) characteristics
- Internal flow field calculation
- Wall contour design for shock-free flow
- Reflection and cancellation of characteristics at the wall

Key parameters:
- num_rays: number of initial characteristic rays (10-60)
- mach_exit: design exit Mach number
- throat_height: throat half-height (in)
- gamma: specific heat ratio

---

## Validation

The MOC results have been compared against:
- Classical isentropic flow tables
- Prandtl-Meyer expansion theory
- Published nozzle contour data

---

## Contact and Collaboration

This project is open for collaboration with aerospace research institutions and industry partners.

For technical inquiries, collaboration proposals, or licensing discussions, please reach out via GitHub Issues or direct contact.

---

## License

MIT License — free for academic and commercial use with attribution.
