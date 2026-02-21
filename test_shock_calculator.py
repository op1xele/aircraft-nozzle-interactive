"""
Test Script for Oblique Shock Calculator - ULTIMATE VERSION
Demonstrates all capabilities including shock diamonds
"""

from core.shock_calculator import ShockCalculator
import numpy as np

print("=" * 70)
print("PICKLE BROTHERS - OBLIQUE SHOCK CALCULATOR TEST")
print("=" * 70)

# Initialize calculator
calc = ShockCalculator(gamma=1.4)

# ============================================================================
# TEST 1: NORMAL SHOCK
# ============================================================================
print("\n" + "=" * 70)
print("TEST 1: NORMAL SHOCK")
print("=" * 70)

M1 = 2.0
result = calc.normal_shock(M1)

print(f"\nUpstream Mach: M1 = {M1}")
print(f"Downstream Mach: M2 = {result['M2']:.4f}")
print(f"Pressure Ratio: p2/p1 = {result['p2_p1']:.4f}")
print(f"Density Ratio: ρ2/ρ1 = {result['rho2_rho1']:.4f}")
print(f"Temperature Ratio: T2/T1 = {result['T2_T1']:.4f}")
print(f"Stagnation Pressure Ratio: p02/p01 = {result['p02_p01']:.4f}")

# ============================================================================
# TEST 2: OBLIQUE SHOCK (GIVEN SHOCK ANGLE)
# ============================================================================
print("\n" + "=" * 70)
print("TEST 2: OBLIQUE SHOCK (Given Shock Angle)")
print("=" * 70)

M1 = 3.0
beta = 35.0  # degrees

result = calc.oblique_shock(M1, beta)

print(f"\nUpstream Mach: M1 = {M1}")
print(f"Shock Angle: β = {beta}°")
print(f"Deflection Angle: θ = {result['theta']:.4f}°")
print(f"Downstream Mach: M2 = {result['M2']:.4f}")
print(f"Pressure Ratio: p2/p1 = {result['p2_p1']:.4f}")

# ============================================================================
# TEST 3: OBLIQUE SHOCK (GIVEN DEFLECTION ANGLE)
# ============================================================================
print("\n" + "=" * 70)
print("TEST 3: OBLIQUE SHOCK (Given Deflection Angle)")
print("=" * 70)

M1 = 2.5
theta = 20.0  # degrees

solutions = calc.theta_beta_M(M1, theta)

print(f"\nUpstream Mach: M1 = {M1}")
print(f"Deflection Angle: θ = {theta}°")
print(f"Maximum Deflection: θ_max = {solutions['theta_max']:.4f}°")

if solutions['weak']:
    print(f"\n--- WEAK SHOCK SOLUTION ---")
    print(f"Shock Angle: β = {solutions['weak']['beta']:.4f}°")
    print(f"Downstream Mach: M2 = {solutions['weak']['M2']:.4f}")

# ============================================================================
# TEST 4: MAXIMUM DEFLECTION ANGLES
# ============================================================================
print("\n" + "=" * 70)
print("TEST 4: MAXIMUM DEFLECTION ANGLES")
print("=" * 70)

print("\n  M1    |  θ_max  |  μ")
print("-" * 30)

for M in [1.5, 2.0, 2.5, 3.0, 4.0]:
    theta_max = calc.get_theta_max(M)
    mu = calc.mach_angle(M)
    print(f"  {M:.1f}  |  {theta_max:.2f}°  |  {mu:.2f}°")

# ============================================================================
# TEST 5: SHOCK DIAMONDS
# ============================================================================
print("\n" + "=" * 70)
print("TEST 5: SHOCK DIAMOND PATTERN")
print("=" * 70)

M_exit = 2.5
p_exit = 30.0  # psia
p_ambient = 14.7  # psia
nozzle_radius = 2.0  # inches

diamond_result = calc.shock_diamond_pattern(M_exit, p_exit, p_ambient, nozzle_radius)

print(f"\nExit Conditions:")
print(f"  M_exit = {M_exit}")
print(f"  p_exit = {p_exit} psia")
print(f"  p_ambient = {p_ambient} psia")
print(f"  Pressure Ratio = {p_exit/p_ambient:.2f}")

print(f"\nShock Diamond Pattern:")
print(f"  Type: {diamond_result['type']}")
print(f"  Number of Diamonds: {len(diamond_result['diamonds'])}")

if diamond_result['diamonds']:
    print(f"\n  Diamond Details:")
    for i, diamond in enumerate(diamond_result['diamonds'][:3]):  # Show first 3
        print(f"    Diamond {i+1}:")
        print(f"      x_center: {diamond['x_center']:.2f} in")
        print(f"      radius: {diamond['radius']:.2f} in")

print("\n" + "=" * 70)
print("✅ ALL TESTS COMPLETED!")
print("=" * 70)