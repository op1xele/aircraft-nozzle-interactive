"""
Structural Analysis Module
Nozzle wall stress analysis under pressure and thermal loads.
"""

import numpy as np


class StructuralAnalysis:

    def __init__(self, material='Inconel 718', wall_thickness_mm=3.0, safety_factor=1.5):
        from core.thermal import MATERIALS
        self.mat = MATERIALS.get(material, MATERIALS['Inconel 718'])
        self.mat_name = material
        self.t = wall_thickness_mm / 1000.0
        self.SF = safety_factor

    def pressure_stresses(self, wall_y_in, p_internal_Pa, p_external_Pa=0.0):
        r = np.array(wall_y_in) * 0.0254
        dp = p_internal_Pa - p_external_Pa
        sigma_hoop  = dp * r / self.t
        sigma_axial = sigma_hoop / 2.0
        return sigma_hoop, sigma_axial

    def thermal_stress(self, delta_T_K):
        E = self.mat['E']; alpha = self.mat['alpha']; nu = 0.3
        return E * alpha * np.array(delta_T_K) / (1.0 - nu)

    def von_mises(self, sigma_hoop, sigma_axial, sigma_thermal=None):
        s1 = sigma_hoop + (sigma_thermal if sigma_thermal is not None else 0)
        s2 = sigma_axial + (sigma_thermal * 0.5 if sigma_thermal is not None else 0)
        return np.sqrt(s1**2 - s1*s2 + s2**2)

    def safety_factors(self, sigma_vm):
        return self.mat['yield_strength'] / np.maximum(sigma_vm, 1.0)

    def fatigue_life(self, sigma_vm, stress_ratio=0.1, n_cycles_design=1000):
        sy = self.mat['yield_strength']
        su = sy * 1.3
        se = 0.5 * su
        sigma_mean = sigma_vm * (1 + stress_ratio) / 2.0
        sigma_alt  = sigma_vm * (1 - stress_ratio) / 2.0
        goodman = np.clip(sigma_alt/se + sigma_mean/su, 1e-6, 10.0)
        N_fail = np.where(goodman < 1.0, 1e9, 1.0/goodman * 1e6)
        utilization = n_cycles_design / np.maximum(N_fail, 1.0) * 100.0
        return N_fail, utilization

    def full_analysis(self, wall_x_in, wall_y_in, p_internal_psia,
                      delta_T_dist=None, p_external_psia=14.7):
        p_int = np.array(p_internal_psia) * 6894.76
        if np.isscalar(p_int):
            p_int = np.full(len(wall_x_in), float(p_int))
        p_ext = p_external_psia * 6894.76

        sigma_h, sigma_a = self.pressure_stresses(wall_y_in, p_int, p_ext)
        sigma_t = self.thermal_stress(delta_T_dist) if delta_T_dist is not None else None
        sigma_vm = self.von_mises(sigma_h, sigma_a, sigma_t)
        sf = self.safety_factors(sigma_vm)
        N_fail, utilization = self.fatigue_life(sigma_vm)

        critical_idx = int(np.argmax(sigma_vm))
        return {
            'sigma_hoop_MPa':     sigma_h / 1e6,
            'sigma_axial_MPa':    sigma_a / 1e6,
            'sigma_thermal_MPa':  sigma_t / 1e6 if sigma_t is not None else np.zeros(len(wall_x_in)),
            'sigma_vm_MPa':       sigma_vm / 1e6,
            'safety_factor':      sf,
            'fatigue_N_fail':     N_fail,
            'fatigue_util_pct':   utilization,
            'yield_strength_MPa': self.mat['yield_strength'] / 1e6,
            'max_vm_MPa':         float(np.max(sigma_vm) / 1e6),
            'min_sf':             float(np.min(sf)),
            'critical_x':         float(wall_x_in[critical_idx]),
            'material_ok':        bool(np.min(sf) >= self.SF),
            'design_sf':          self.SF,
        }