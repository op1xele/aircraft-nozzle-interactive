"""
Isentropic Relations Module
"""

import numpy as np
from config import CONV_DEG_RAD, CONV_RAD_DEG, TOLERANCE, MAX_ITERATIONS


class IsentropicRelations:
    """Isentropic relations for compressible flow"""
    
    def __init__(self, gamma=1.4):
        self.gamma = gamma
        self.gamma_m1 = gamma - 1.0
        self.gamma_p1 = gamma + 1.0
        self._cache = {}
    
    def get_all_properties(self, mach):
        """Calculate all isentropic properties"""
        if mach <= 0:
            raise ValueError("Mach must be positive")
        
        cache_key = f"mach_{mach:.6f}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        mach2 = mach * mach
        mach2_m1 = mach2 - 1.0
        
        # Temperature ratio
        factor = 1.0 + 0.5 * self.gamma_m1 * mach2
        T_Tt = 1.0 / factor
        
        # Pressure ratio
        p_pt = np.power(T_Tt, self.gamma / self.gamma_m1)
        
        # Density ratio
        rho_rhot = np.power(T_Tt, 1.0 / self.gamma_m1)
        
        # Area ratio
        factor2 = self.gamma_p1 / (2.0 * self.gamma_m1)
        A_Astar = (1.0 / mach) * np.power(factor / (self.gamma_p1 / 2.0), factor2)
        
        # Mach angle
        if mach >= 1.0:
            mu = np.arcsin(1.0 / mach) * CONV_RAD_DEG
        else:
            mu = 90.0
        
        # Prandtl-Meyer angle
        if mach > 1.0:
            nu = self._calculate_prandtl_meyer(mach)
        else:
            nu = 0.0
        
        result = {
            'p_pt': p_pt,
            'T_Tt': T_Tt,
            'rho_rhot': rho_rhot,
            'A_Astar': A_Astar,
            'mu': mu,
            'nu': nu,
            'mach': mach
        }
        
        self._cache[cache_key] = result
        return result
    
    def _calculate_prandtl_meyer(self, mach):
        """Calculate Prandtl-Meyer angle"""
        if mach <= 1.0:
            return 0.0
        
        mach2 = mach * mach
        mach2_m1 = mach2 - 1.0
        
        sqrt_factor = np.sqrt(self.gamma_p1 / self.gamma_m1)
        inner_sqrt = np.sqrt(self.gamma_m1 * mach2_m1 / self.gamma_p1)
        
        nu_rad = sqrt_factor * np.arctan(inner_sqrt) - np.arctan(np.sqrt(mach2_m1))
        nu_deg = nu_rad * CONV_RAD_DEG
        
        return nu_deg
    
    def get_mach_from_pm(self, nu_deg):
        """Get Mach from Prandtl-Meyer angle"""
        if nu_deg <= 0:
            return 1.0
        
        nu_max = (np.sqrt(self.gamma_p1 / self.gamma_m1) - 1.0) * 90.0
        if nu_deg >= nu_max:
            return 100.0
        
        nu_rad = nu_deg * CONV_DEG_RAD
        msm1 = 1.0
        
        for iteration in range(MAX_ITERATIONS):
            sqrt_factor = np.sqrt(self.gamma_p1 / self.gamma_m1)
            inner_sqrt = np.sqrt(self.gamma_m1 * msm1 / self.gamma_p1)
            
            nu_current = sqrt_factor * np.arctan(inner_sqrt) - np.arctan(np.sqrt(msm1))
            
            error = nu_rad - nu_current
            if abs(error) < TOLERANCE:
                break
            
            dnu_dmsm1 = (sqrt_factor / (1.0 + self.gamma_m1 * msm1 / self.gamma_p1) * 
                        self.gamma_m1 / (2.0 * self.gamma_p1 * np.sqrt(msm1)) -
                        1.0 / (2.0 * msm1 * (1.0 + msm1)))
            
            msm1_new = msm1 + error / dnu_dmsm1
            if msm1_new < 0:
                msm1_new = msm1 / 2.0
            
            msm1 = msm1_new
        
        return np.sqrt(msm1 + 1.0)
    
    def clear_cache(self):
        self._cache.clear()


def get_isentropic_properties(mach, gamma=1.4):
    calc = IsentropicRelations(gamma)
    return calc.get_all_properties(mach)


def get_mach_from_prandtl_meyer(nu, gamma=1.4):
    calc = IsentropicRelations(gamma)
    return calc.get_mach_from_pm(nu)