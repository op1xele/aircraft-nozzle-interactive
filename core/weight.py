"""
Weight Estimation Module
Nozzle mass estimation based on geometry and material.
"""

import numpy as np


class WeightEstimation:
    """
    Nozzle weight estimation using shell/frustum geometry
    and material density.
    """

    def __init__(self, material='Inconel 718', wall_thickness_mm=3.0):
        from core.thermal import MATERIALS
        self.mat = MATERIALS.get(material, MATERIALS['Inconel 718'])
        self.mat_name = material
        self.t = wall_thickness_mm / 1000.0  # m

    def nozzle_mass(self, wall_x_in, wall_y_in):
        """
        Compute nozzle wall mass using frustum shell integration.
        Revolves the 2D wall profile around the axis.

        Returns: mass_kg, mass_lbm, volume_m3, surface_area_m2
        """
        x = np.array(wall_x_in) * 0.0254  # m
        r = np.array(wall_y_in) * 0.0254  # m

        surface_area = 0.0
        for i in range(1, len(x)):
            r1, r2 = r[i-1], r[i]
            dx = x[i] - x[i-1]
            dr = r2 - r1
            slant = np.sqrt(dx**2 + dr**2)
            # Frustum lateral surface area
            surface_area += np.pi * (r1 + r2) * slant

        volume = surface_area * self.t
        mass_kg = volume * self.mat['density']
        mass_lbm = mass_kg * 2.20462

        return {
            'mass_kg':        mass_kg,
            'mass_lbm':       mass_lbm,
            'volume_m3':      volume,
            'surface_area_m2': surface_area,
            'material':       self.mat_name,
            'density_kg_m3':  self.mat['density'],
            'wall_thickness_mm': self.t * 1000,
        }

    def mass_breakdown(self, wall_x_in, wall_y_in):
        """
        Mass breakdown by nozzle section:
        convergent, throat, divergent.
        """
        x = np.array(wall_x_in)
        y = np.array(wall_y_in)

        throat_idx = int(np.argmin(y))
        mid_idx = (throat_idx + len(x)) // 2

        sections = {
            'Convergent': (0, throat_idx),
            'Throat':     (max(0, throat_idx-2), min(len(x)-1, throat_idx+2)),
            'Divergent':  (throat_idx, len(x)-1),
        }

        breakdown = {}
        for name, (i0, i1) in sections.items():
            if i1 > i0:
                result = self.nozzle_mass(x[i0:i1+1], y[i0:i1+1])
                breakdown[name] = result['mass_kg']

        total = sum(breakdown.values())
        breakdown['Total'] = total
        return breakdown

    def compare_materials(self, wall_x_in, wall_y_in):
        """Compare mass across all available materials."""
        from core.thermal import MATERIALS
        results = {}
        for mat_name in MATERIALS:
            we = WeightEstimation(mat_name, self.t * 1000)
            r = we.nozzle_mass(wall_x_in, wall_y_in)
            results[mat_name] = {
                'mass_kg':   r['mass_kg'],
                'mass_lbm':  r['mass_lbm'],
                'density':   MATERIALS[mat_name]['density'],
            }
        return results

    def thrust_to_weight(self, mass_kg, thrust_N):
        """Nozzle T/W ratio."""
        weight_N = mass_kg * 9.81
        return thrust_N / weight_N if weight_N > 0 else 0.0