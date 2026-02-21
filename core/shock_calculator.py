"""
Oblique Shock Calculator Module - ULTIMATE VERSION
Complete oblique and normal shock relations + Shock Diamonds
Based on Anderson's Compressible Flow
"""

import numpy as np
from typing import Tuple, Dict, Optional

# Try to import from config, fallback to constants if not available
try:
    from config import CONV_DEG_RAD, CONV_RAD_DEG, TOLERANCE, MAX_ITERATIONS
except:
    CONV_DEG_RAD = np.pi / 180.0
    CONV_RAD_DEG = 180.0 / np.pi
    TOLERANCE = 1e-6
    MAX_ITERATIONS = 100


class ShockCalculator:
    """
    Comprehensive shock wave calculator
    - Normal shock relations
    - Oblique shock relations
    - Shock angle solver (beta-theta-M)
    - Weak/strong shock solutions
    - Detachment angle
    - Shock polar calculations
    - Shock diamond patterns (plume)
    """
    
    def __init__(self, gamma: float = 1.4):
        """
        Initialize shock calculator
        
        Args:
            gamma: Ratio of specific heats (default 1.4 for air)
        """
        self.gamma = gamma
        self.gamma_m1 = gamma - 1.0
        self.gamma_p1 = gamma + 1.0
        self.gm1_gp1 = self.gamma_m1 / self.gamma_p1
        
    def normal_shock(self, M1: float) -> Dict[str, float]:
        """
        Calculate normal shock relations
        
        Args:
            M1: Upstream Mach number
            
        Returns:
            Dictionary with all normal shock properties:
                M2: Downstream Mach number
                p2_p1: Pressure ratio
                rho2_rho1: Density ratio
                T2_T1: Temperature ratio
                p02_p01: Stagnation pressure ratio
                pt2_pt1: Same as p02_p01 (compatibility)
                u2_u1: Velocity ratio
        """
        if M1 <= 1.0:
            raise ValueError(f"Normal shock requires M1 > 1.0, got M1={M1}")
        
        M1_sq = M1 * M1
        
        # Downstream Mach number
        numerator = 1.0 + 0.5 * self.gamma_m1 * M1_sq
        denominator = self.gamma * M1_sq - 0.5 * self.gamma_m1
        M2_sq = numerator / denominator
        M2 = np.sqrt(M2_sq)
        
        # Pressure ratio
        p2_p1 = 1.0 + (2.0 * self.gamma / self.gamma_p1) * (M1_sq - 1.0)
        
        # Density ratio
        rho2_rho1 = self.gamma_p1 * M1_sq / (self.gamma_m1 * M1_sq + 2.0)
        
        # Temperature ratio
        T2_T1 = p2_p1 / rho2_rho1
        
        # Stagnation pressure ratio (loss)
        term1 = (self.gamma_p1 * M1_sq / (self.gamma_m1 * M1_sq + 2.0)) ** (self.gamma / self.gamma_m1)
        term2 = (self.gamma_p1 / (2.0 * self.gamma * M1_sq - self.gamma_m1)) ** (1.0 / self.gamma_m1)
        p02_p01 = term1 * term2
        
        # Velocity ratio
        u2_u1 = rho2_rho1 ** (-1.0)
        
        return {
            'M2': M2,
            'p2_p1': p2_p1,
            'rho2_rho1': rho2_rho1,
            'T2_T1': T2_T1,
            'p02_p01': p02_p01,
            'pt2_pt1': p02_p01,  # Compatibility with old code
            'u2_u1': u2_u1,
            'type': 'normal'
        }
    
    def oblique_shock(self, M1: float, beta: float, return_both: bool = False) -> Dict[str, float]:
        """
        Calculate oblique shock relations given shock angle
        
        Args:
            M1: Upstream Mach number
            beta: Shock angle in degrees
            return_both: If True, return both weak and strong solutions
            
        Returns:
            Dictionary with oblique shock properties
        """
        if M1 <= 1.0:
            raise ValueError(f"Oblique shock requires M1 > 1.0, got M1={M1}")
        
        beta_rad = np.radians(beta)
        
        # Normal component of Mach number
        Mn1 = M1 * np.sin(beta_rad)
        
        if Mn1 <= 1.0:
            raise ValueError(f"Normal Mach component must be > 1.0, got Mn1={Mn1:.3f}")
        
        # Normal shock relations for normal component
        normal_results = self.normal_shock(Mn1)
        Mn2 = normal_results['M2']
        
        # Tangential component (unchanged across shock)
        Mt1 = M1 * np.cos(beta_rad)
        Mt2 = Mt1  # Tangential component unchanged
        
        # Downstream Mach number
        M2 = np.sqrt(Mn2**2 + Mt2**2)
        
        # Flow deflection angle (theta)
        tan_theta = 2.0 / np.tan(beta_rad) * (M1**2 * np.sin(beta_rad)**2 - 1.0)
        tan_theta /= (M1**2 * (self.gamma + np.cos(2.0 * beta_rad)) + 2.0)
        theta = np.degrees(np.arctan(tan_theta))
        
        return {
            'M1': M1,
            'M2': M2,
            'Mn1': Mn1,
            'Mn2': Mn2,
            'M1n': Mn1,  # Compatibility
            'M2n': Mn2,  # Compatibility
            'beta': beta,
            'theta': theta,
            'p2_p1': normal_results['p2_p1'],
            'rho2_rho1': normal_results['rho2_rho1'],
            'T2_T1': normal_results['T2_T1'],
            'p02_p01': normal_results['p02_p01'],
            'pt2_pt1': normal_results['p02_p01'],  # Compatibility
            'type': 'oblique'
        }
    
    def oblique_shock_properties(self, M1: float, beta_deg: float) -> Dict[str, float]:
        """
        Calculate oblique shock properties for given β (compatibility method)
        
        Args:
            M1: Upstream Mach number
            beta_deg: Shock angle β (degrees)
            
        Returns:
            dict: Flow properties across oblique shock
        """
        return self.oblique_shock(M1, beta_deg)
    
    def theta_beta_M(self, M1: float, theta: float) -> Dict[str, float]:
        """
        Solve for shock angle given deflection angle (theta-beta-M relation)
        Returns both weak and strong shock solutions
        
        Args:
            M1: Upstream Mach number
            theta: Flow deflection angle in degrees
            
        Returns:
            Dictionary with weak and strong shock solutions
        """
        if M1 <= 1.0:
            raise ValueError(f"Oblique shock requires M1 > 1.0, got M1={M1}")
        
        theta_rad = np.radians(theta)
        
        # Check if solution exists (theta < theta_max)
        theta_max = self.get_theta_max(M1)
        if abs(theta) > theta_max:
            return {
                'M1': M1,
                'theta': theta,
                'theta_max': theta_max,
                'mu': self.mach_angle(M1),
                'weak': None,
                'strong': None,
                'beta_weak': None,
                'beta_strong': None,
                'detached': True
            }
        
        # Solve theta-beta-M equation iteratively
        def theta_beta_eq(beta_deg):
            beta_rad = np.radians(beta_deg)
            M1_sq = M1 * M1
            sin_beta_sq = np.sin(beta_rad) ** 2
            
            numerator = 2.0 * (M1_sq * sin_beta_sq - 1.0) / np.tan(beta_rad)
            denominator = M1_sq * (self.gamma + np.cos(2.0 * beta_rad)) + 2.0
            
            return np.degrees(np.arctan(numerator / denominator)) - theta
        
        # Weak shock solution (smaller beta)
        mu = self.mach_angle(M1)  # Mach angle
        
        try:
            from scipy.optimize import brentq
            
            # Weak shock: between mu and ~70 degrees
            beta_weak = brentq(theta_beta_eq, mu + 0.1, 85.0)
            weak_solution = self.oblique_shock(M1, beta_weak)
        except:
            # Fallback: Newton-Raphson
            beta_weak = self._solve_theta_beta_mach(M1, theta, weak=True)
            weak_solution = self.oblique_shock(M1, beta_weak) if beta_weak else None
        
        try:
            from scipy.optimize import brentq
            
            # Strong shock: between ~70 and 90 degrees
            beta_strong = brentq(theta_beta_eq, 70.0, 89.9)
            strong_solution = self.oblique_shock(M1, beta_strong)
        except:
            # Fallback: Newton-Raphson
            beta_strong = self._solve_theta_beta_mach(M1, theta, weak=False)
            strong_solution = self.oblique_shock(M1, beta_strong) if beta_strong else None
        
        return {
            'M1': M1,
            'theta': theta,
            'theta_max': theta_max,
            'mu': mu,
            'weak': weak_solution,
            'strong': strong_solution,
            'beta_weak': beta_weak if weak_solution else None,
            'beta_strong': beta_strong if strong_solution else None,
            'detached': False
        }
    
    def oblique_shock_beta(self, M1: float, theta_deg: float) -> Dict[str, float]:
        """
        Calculate oblique shock angle β for given M1 and deflection angle θ
        (Compatibility method with old interface)
        
        Args:
            M1: Upstream Mach number
            theta_deg: Flow deflection angle (degrees)
            
        Returns:
            dict: Beta angles (weak and strong solutions)
        """
        result = self.theta_beta_M(M1, theta_deg)
        
        return {
            'beta_weak': result['beta_weak'],
            'beta_strong': result['beta_strong'],
            'detached': result['detached'],
            'theta_max': result['theta_max']
        }
    
    def _solve_theta_beta_mach(self, M1: float, theta_deg: float, weak: bool = True) -> Optional[float]:
        """
        Solve theta-beta-Mach relation for β using Newton-Raphson
        
        Args:
            M1: Upstream Mach number
            theta_deg: Deflection angle in degrees
            weak: True for weak shock, False for strong shock
            
        Returns:
            Beta angle in degrees, or None if failed
        """
        theta_rad = theta_deg * CONV_DEG_RAD
        M1_sq = M1 * M1
        
        # Initial guess
        if weak:
            beta = np.arcsin(1.0 / M1) + 0.1  # Start just above Mach angle
        else:
            beta = np.pi / 2.0 - 0.1  # Start just below 90 degrees
        
        for iteration in range(MAX_ITERATIONS):
            sin_beta = np.sin(beta)
            cos_beta = np.cos(beta)
            tan_beta = np.tan(beta)
            
            # Theta-beta-Mach relation
            numerator = M1_sq * sin_beta**2 - 1.0
            denominator = M1_sq * (self.gamma + np.cos(2.0 * beta)) + 2.0
            tan_theta_calc = 2.0 / tan_beta * numerator / denominator
            theta_calc = np.arctan(tan_theta_calc)
            
            error = theta_calc - theta_rad
            
            if abs(error) < TOLERANCE:
                return beta * CONV_RAD_DEG
            
            # Numerical derivative
            delta = 0.0001
            beta_plus = beta + delta
            numerator_plus = M1_sq * np.sin(beta_plus)**2 - 1.0
            denominator_plus = M1_sq * (self.gamma + np.cos(2.0 * beta_plus)) + 2.0
            tan_theta_plus = 2.0 / np.tan(beta_plus) * numerator_plus / denominator_plus
            theta_plus = np.arctan(tan_theta_plus)
            
            derivative = (theta_plus - theta_calc) / delta
            
            if abs(derivative) < 1e-10:
                break
            
            # Update beta
            beta = beta - error / derivative
            
            # Keep in valid range
            mu = np.arcsin(1.0 / M1)  # Mach angle
            if weak:
                beta = np.clip(beta, mu + 0.01, np.pi / 2.0 - 0.01)
            else:
                beta = np.clip(beta, mu + 0.1, np.pi / 2.0 - 0.01)
        
        return None  # Failed to converge
    
    def get_theta_max(self, M1: float) -> float:
        """
        Calculate maximum deflection angle for attached shock
        
        Args:
            M1: Upstream Mach number
            
        Returns:
            Maximum deflection angle in degrees
        """
        if M1 <= 1.0:
            return 0.0
        
        # More accurate: numerical search
        M1_sq = M1 * M1
        theta_max = 0.0
        
        mu = np.arcsin(1.0 / M1)  # Mach angle in radians
        
        for beta_rad in np.linspace(mu + 0.01, np.pi / 2.0 - 0.01, 200):
            try:
                sin_beta_sq = np.sin(beta_rad) ** 2
                numerator = 2.0 * (M1_sq * sin_beta_sq - 1.0) / np.tan(beta_rad)
                denominator = M1_sq * (self.gamma + np.cos(2.0 * beta_rad)) + 2.0
                
                theta_rad = np.arctan(numerator / denominator)
                theta_deg = np.degrees(theta_rad)
                
                if theta_deg > theta_max:
                    theta_max = theta_deg
            except:
                continue
        
        return theta_max
    
    def _get_max_deflection_angle(self, M1: float) -> float:
        """Compatibility wrapper"""
        return self.get_theta_max(M1)
    
    def _get_theta_from_beta(self, M1: float, beta_deg: float) -> float:
        """Helper: calculate theta from beta"""
        beta_rad = np.radians(beta_deg)
        M1_sq = M1 * M1
        sin_beta_sq = np.sin(beta_rad) ** 2
        
        numerator = 2.0 * (M1_sq * sin_beta_sq - 1.0) / np.tan(beta_rad)
        denominator = M1_sq * (self.gamma + np.cos(2.0 * beta_rad)) + 2.0
        
        theta = np.arctan(numerator / denominator)
        return np.degrees(theta)
    
    def detachment_angle(self, M1: float) -> float:
        """
        Calculate detachment (maximum) deflection angle
        Same as get_theta_max but with clearer name
        
        Args:
            M1: Upstream Mach number
            
        Returns:
            Detachment angle in degrees
        """
        return self.get_theta_max(M1)
    
    def shock_from_deflection(self, M1: float, theta: float, 
                             solution: str = 'weak') -> Dict[str, float]:
        """
        Complete shock solution from deflection angle
        
        Args:
            M1: Upstream Mach number
            theta: Deflection angle in degrees
            solution: 'weak' or 'strong' shock solution
            
        Returns:
            Complete shock properties
        """
        solutions = self.theta_beta_M(M1, theta)
        
        if solutions['detached']:
            raise ValueError(f"Shock is detached: theta={theta:.2f}° > theta_max={solutions['theta_max']:.2f}°")
        
        if solution == 'weak' and solutions['weak'] is not None:
            return solutions['weak']
        elif solution == 'strong' and solutions['strong'] is not None:
            return solutions['strong']
        else:
            raise ValueError(f"No {solution} shock solution available")
    
    def mach_angle(self, M: float) -> float:
        """
        Calculate Mach angle
        
        Args:
            M: Mach number
            
        Returns:
            Mach angle in degrees
        """
        if M < 1.0:
            raise ValueError(f"Mach angle undefined for M < 1.0, got M={M}")
        
        mu = np.degrees(np.arcsin(1.0 / M))
        return mu
    
    def shock_polar(self, M1: float, num_points: int = 100) -> Dict[str, np.ndarray]:
        """
        Generate shock polar curve (pressure vs deflection)
        
        Args:
            M1: Upstream Mach number
            num_points: Number of points on curve
            
        Returns:
            Dictionary with theta and pressure ratio arrays
        """
        mu = self.mach_angle(M1)
        theta_max = self.get_theta_max(M1)
        
        # Generate beta values from Mach angle to 90 degrees
        beta_vals = np.linspace(mu + 0.1, 89.9, num_points)
        
        theta_vals = []
        p_ratio_vals = []
        
        for beta in beta_vals:
            try:
                result = self.oblique_shock(M1, beta)
                theta_vals.append(result['theta'])
                p_ratio_vals.append(result['p2_p1'])
            except:
                continue
        
        return {
            'theta': np.array(theta_vals),
            'p2_p1': np.array(p_ratio_vals),
            'M1': M1,
            'theta_max': theta_max
        }
    
    def get_shock_properties(self, M1: float, p1: float, T1: float, 
                           theta: float = None, beta: float = None,
                           solution: str = 'weak') -> Dict[str, float]:
        """
        Get complete shock properties with dimensional quantities
        
        Args:
            M1: Upstream Mach number
            p1: Upstream static pressure
            T1: Upstream static temperature
            theta: Deflection angle (if known)
            beta: Shock angle (if known)
            solution: 'weak' or 'strong' for oblique shock
            
        Returns:
            Complete shock properties with dimensional values
        """
        if theta is not None:
            # Given deflection angle
            result = self.shock_from_deflection(M1, theta, solution)
        elif beta is not None:
            # Given shock angle
            result = self.oblique_shock(M1, beta)
        else:
            # Normal shock
            result = self.normal_shock(M1)
        
        # Add dimensional properties
        result['p1'] = p1
        result['p2'] = p1 * result['p2_p1']
        result['T1'] = T1
        result['T2'] = T1 * result['T2_T1']
        
        return result
    
    # ========================================================================
    # SHOCK DIAMOND PATTERNS (from old version)
    # ========================================================================
    
    def shock_diamond_pattern(self, M_exit: float, p_exit: float, 
                             p_ambient: float, nozzle_radius: float) -> Dict:
        """
        Calculate shock diamond pattern in plume
        For under-expanded or over-expanded flows
        
        Args:
            M_exit: Exit Mach number
            p_exit: Exit pressure
            p_ambient: Ambient pressure
            nozzle_radius: Nozzle exit radius
            
        Returns:
            dict: Shock diamond geometry
        """
        # Pressure ratio determines expansion type
        pressure_ratio = p_exit / p_ambient
        
        if abs(pressure_ratio - 1.0) < 0.05:
            # Perfectly expanded - no shocks
            return {
                'type': 'perfectly_expanded',
                'diamonds': []
            }
        
        elif pressure_ratio > 1.0:
            # Under-expanded - expansion waves and oblique shocks
            return self._calculate_underexpanded_diamonds(
                M_exit, pressure_ratio, nozzle_radius
            )
        
        else:
            # Over-expanded - oblique shocks at exit
            return self._calculate_overexpanded_diamonds(
                M_exit, pressure_ratio, nozzle_radius
            )
    
    def _calculate_underexpanded_diamonds(self, M_exit: float, PR: float, R: float) -> Dict:
        """Calculate shock diamonds for under-expanded flow"""
        # Simplified shock diamond calculation
        diamonds = []
        
        # First diamond length (approximate)
        L1 = R * M_exit * 2.0
        
        # Number of visible diamonds (decreases with distance)
        num_diamonds = min(5, int(M_exit))
        
        x_pos = 0.0
        for i in range(num_diamonds):
            diamond_length = L1 * (0.8 ** i)  # Each diamond gets shorter
            
            diamonds.append({
                'x_start': x_pos,
                'x_end': x_pos + diamond_length,
                'x_center': x_pos + diamond_length / 2.0,
                'radius': R * (1.0 + 0.2 * i)  # Slight expansion
            })
            
            x_pos += diamond_length
        
        return {
            'type': 'under_expanded',
            'diamonds': diamonds,
            'pressure_ratio': PR
        }
    
    def _calculate_overexpanded_diamonds(self, M_exit: float, PR: float, R: float) -> Dict:
        """Calculate shock pattern for over-expanded flow"""
        # Oblique shock at exit
        # Simplified - single shock reflection
        
        return {
            'type': 'over_expanded',
            'diamonds': [],
            'shock_angle': 30.0,  # Approximate
            'pressure_ratio': PR
        }


# ============================================================================
# UTILITY FUNCTIONS FOR MOC ENGINE INTEGRATION
# ============================================================================

def calculate_oblique_shock(M1: float, theta: float, gamma: float = 1.4) -> Tuple[float, float, float]:
    """
    Quick oblique shock calculation (weak shock solution)
    Returns (M2, beta, pressure_ratio)
    """
    calc = ShockCalculator(gamma)
    result = calc.shock_from_deflection(M1, theta, solution='weak')
    return result['M2'], result['beta'], result['p2_p1']


def calculate_normal_shock(M1: float, gamma: float = 1.4) -> Tuple[float, float]:
    """
    Quick normal shock calculation
    Returns (M2, pressure_ratio)
    """
    calc = ShockCalculator(gamma)
    result = calc.normal_shock(M1)
    return result['M2'], result['p2_p1']


def check_shock_detachment(M1: float, theta: float, gamma: float = 1.4) -> bool:
    """
    Check if shock is detached for given conditions
    Returns True if detached
    """
    calc = ShockCalculator(gamma)
    theta_max = calc.get_theta_max(M1)
    return abs(theta) > theta_max