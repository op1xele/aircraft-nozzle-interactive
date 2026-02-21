# Aircraft Nozzle Interactive — Technical Documentation

**Version:** 2.0  
**Method:** Method of Characteristics (MOC)  
**Language:** Python 3.10+

---

## 1. Theoretical Background

### 1.1 Method of Characteristics (MOC)

The Method of Characteristics is an exact mathematical technique for solving hyperbolic partial differential equations governing supersonic flow. For two-dimensional irrotational isentropic flow, the governing equations reduce to characteristic compatibility equations.

**Characteristic equations (2D):**

Along C+ characteristics (right-running):
```
dθ + dν = 0   →   θ + ν = const (K+)
```

Along C− characteristics (left-running):
```
dθ − dν = 0   →   θ − ν = const (K−)
```

Where:
- `θ` = flow angle (rad)
- `ν` = Prandtl-Meyer function (rad)
- `K+`, `K−` = Riemann invariants

**Prandtl-Meyer function:**
```
ν(M) = √((γ+1)/(γ-1)) · arctan(√((γ-1)/(γ+1)·(M²-1))) − arctan(√(M²-1))
```

**Characteristic slope angles:**
```
C+: dy/dx = tan(θ − μ)
C−: dy/dx = tan(θ + μ)
```

Where `μ = arcsin(1/M)` is the Mach angle.

### 1.2 Isentropic Flow Relations

For a calorically perfect gas:

```
T/T₀ = (1 + (γ-1)/2 · M²)⁻¹
p/p₀ = (T/T₀)^(γ/(γ-1))
ρ/ρ₀ = (T/T₀)^(1/(γ-1))
A/A* = (1/M) · [(2/(γ+1)) · (1 + (γ-1)/2 · M²)]^((γ+1)/(2(γ-1)))
```

### 1.3 Real Gas Effects

Temperature-dependent gamma is computed using NASA 9-coefficient thermodynamic polynomials (McBride & Gordon 1994):

```
Cp/R = a₁/T² + a₂/T + a₃ + a₄T + a₅T² + a₆T³ + a₇T⁴
```

Valid ranges: 200–1000 K (low temperature), 1000–6000 K (high temperature)

The iterative solution at each Mach number:
```
T_k = T₀ / (1 + 0.5·(γ(T_k)-1)·M²)
```
Converges in ~30 iterations to within 0.1 K tolerance.

---

## 2. Numerical Implementation

### 2.1 MOC Solver (core/moc_engine.py)

**Initial data line:**  
The calculation begins at the throat where the flow is sonic (M=1). Initial rays are distributed from the wall to the centerline with uniform angular spacing.

**Interior point calculation:**  
For each interior point P defined by the intersection of a C+ characteristic from point A and a C− characteristic from point B:
```
K+ (from A) = θ_A + ν_A
K− (from B) = θ_B − ν_B

θ_P = (K+ + K−) / 2
ν_P = (K+ − K−) / 2
M_P = ν⁻¹(ν_P)
```

**Wall point calculation:**  
At the wall, the flow tangency condition `θ_wall = wall_slope` is applied. Combined with the incoming characteristic invariant, this gives the wall Mach number.

**Centerline point calculation:**  
Symmetry condition `θ = 0` applied on the centerline.

### 2.2 Viscous Correction (core/viscous.py)

Boundary layer displacement thickness `δ*` is computed using:

**Laminar (Blasius):**
```
δ* = 1.72 · x / √Re_x
```

**Turbulent (1/7 power law with Van Driest compressibility):**
```
δ* = 0.046 · x / Re_x^0.2 · (T_wall/T_e)^0.61
```

**Transition criterion (Michel):** Re_x > 5×10⁵

**Axisymmetric Mangler correction:**
```
δ*_axisym = δ*_2D · √(r_throat / r_local)
```

The physical wall is displaced outward by `δ*` to produce the correct inviscid flow area.

### 2.3 Thermal Analysis (core/thermal.py)

**Bartz correlation** for heat transfer coefficient:
```
h = (0.026/D_t^0.2) · (μ^0.2·Cp/Pr^0.6) · (p₀/c*)^0.8 · (D_t/R_c)^0.1 · (A_t/A)^0.9 · σ
```

**Wall temperature (1D radial conduction):**
```
q = (T_aw − T_coolant) / (1/h + t/k)
T_hot = T_aw − q/h
T_cold = T_hot − q·t/k
```

**Adiabatic wall temperature:**
```
T_aw = T₀ · (T/T₀ + r·(γ-1)/2·M²·T/T₀)
```
Where `r = Pr^(1/3)` is the recovery factor.

### 2.4 Structural Analysis (core/structural.py)

**Thin-wall pressure vessel (hoop stress):**
```
σ_hoop = (p_i − p_e) · r / t
σ_axial = σ_hoop / 2
```

**Thermal stress (biaxial constraint):**
```
σ_T = E · α · ΔT / (1 − ν)
```

**Von Mises criterion:**
```
σ_VM = √(σ₁² − σ₁·σ₂ + σ₂²)
```

**Fatigue life (Goodman criterion):**
```
σ_alt/Se + σ_mean/Su = 1/N_f
```

### 2.5 Performance Calculations (core/performance.py)

**Thrust:**
```
F = ṁ·V_e + (p_e − p_a)·A_e
```

**Specific impulse:**
```
Isp = F / (ṁ·g₀)
```

**Characteristic velocity:**
```
c* = p₀·A_t / ṁ
```

**Thrust coefficient:**
```
C_F = F / (p₀·A_t)
```

---

## 3. Nozzle Design Methods

### 3.1 Minimum Length Nozzle (prob=8)

Designed to achieve uniform parallel flow at the exit in the shortest possible length. The expansion section is designed using MOC with the Prandtl-Meyer angle equal to the maximum wall angle at the throat.

### 3.2 Bell Nozzle (prob=6)

Rao's approximate method: parabolic approximation of the optimal contour. Parameters:
- Initial wall angle at throat: `θ_i`  
- Exit wall angle: `θ_e`
- Length fraction relative to 15° cone

### 3.3 Plug/Aerospike Nozzle (prob=3, prob=5)

C.C. Lee method with `perint` parameter controlling internal/external expansion split:
```
φ_internal = ν_exit · perint
φ_external = ν_exit · (1 − perint)
```

### 3.4 2D Ideal Nozzle (prob=1)

Full MOC solution with:
- Uniform flow at throat (sonic line approximation)
- Characteristic net generation from throat to exit
- Wall shaped to produce uniform parallel supersonic flow at exit

---

## 4. Optimization Algorithms (core/optimization.py)

### 4.1 Nelder-Mead Simplex

Direct search method, no gradients required. Parameters: α=1.0 (reflection), γ=2.0 (expansion), ρ=0.5 (contraction), σ=0.5 (shrink).

Convergence criterion: objective function < tolerance (default 1×10⁻⁴)

### 4.2 Genetic Algorithm

Population-based global search:
- Population size: 20
- Selection: top 25% elitism
- Crossover: arithmetic (random α)
- Mutation rate: 20%, Gaussian noise σ=0.1
- Generations: configurable (default 30)

Design variables: [Mach_exit, num_rays]  
Objective: weighted sum of (Isp_error² + efficiency_error²)

---

## 5. Material Properties

| Material | ρ (kg/m³) | k (W/mK) | E (GPa) | σ_y (MPa) | T_max (K) |
|----------|-----------|----------|---------|-----------|-----------|
| Inconel 718 | 8200 | 11.4 | 200 | 1100 | 1200 |
| Copper C18150 | 8900 | 350 | 115 | 310 | 700 |
| Ti-6Al-4V | 4430 | 6.7 | 114 | 880 | 600 |
| Carbon-Carbon | 1800 | 40 | 70 | 200 | 2500 |
| Stainless 316L | 8000 | 16 | 193 | 290 | 870 |
| Rhenium | 21020 | 48 | 460 | 1070 | 2200 |

---

## 6. References

1. Anderson, J.D. (2003). *Modern Compressible Flow*. McGraw-Hill.
2. Shapiro, A.H. (1953). *The Dynamics and Thermodynamics of Compressible Fluid Flow*. Ronald Press.
3. Bartz, D.R. (1957). "A Simple Equation for Rapid Estimation of Rocket Nozzle Convective Heat Transfer Coefficients." *Jet Propulsion*, 27(1).
4. McBride, B.J. & Gordon, S. (1994). *Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications*. NASA RP-1311.
5. Rao, G.V.R. (1958). "Exhaust Nozzle Contour for Optimum Thrust." *Jet Propulsion*, 28(6).
6. Lee, C.C. (1963). "Numerical Analysis of Plug Nozzles." NASA TR.
7. Thwaites, B. (1949). "Approximate Calculation of the Laminar Boundary Layer." *Aeronautical Quarterly*, 1.
8. Goodman, J. (1899). *Mechanics Applied to Engineering*. Longmans Green.
