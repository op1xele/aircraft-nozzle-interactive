"""
External Flow Calculator Module
External expansion, cowl, afterbody calculations
"""

import numpy as np
from config import CONV_DEG_RAD, CONV_RAD_DEG


class ExternalFlowCalculator:
    """
    Calculate external flow fields around nozzle
    - External expansion beyond nozzle exit
    - Cowl integration and design
    - Afterbody pressure distribution
    - External streamlines
    - Drag calculations
    """
    
    def __init__(self, gamma=1.4):
        self.gamma = gamma
        self.gamma_m1 = gamma - 1.0
        self.gamma_p1 = gamma + 1.0
    
    def calculate_external_expansion(self, exit_mach, exit_radius, exit_angle, 
                                    p_exit, p_ambient, num_streamlines=10):
        """
        Calculate external expansion field beyond nozzle exit
        
        Args:
            exit_mach: Exit Mach number
            exit_radius: Exit radius (inches)
            exit_angle: Exit flow angle (degrees)
            p_exit: Exit static pressure (psia)
            p_ambient: Ambient pressure (psia)
            num_streamlines: Number of external streamlines
        
        Returns:
            dict: External expansion field data
        """
        expansion_type = self._determine_expansion_type(p_exit, p_ambient)
        
        if expansion_type == 'over_expanded':
            # Over-expanded: external compression
            return self._calculate_external_compression(
                exit_mach, exit_radius, exit_angle, p_exit, p_ambient, num_streamlines
            )
        elif expansion_type == 'under_expanded':
            # Under-expanded: external expansion
            return self._calculate_external_prandtl_meyer(
                exit_mach, exit_radius, exit_angle, p_exit, p_ambient, num_streamlines
            )
        else:
            # Perfectly expanded: minimal external flow
            return self._calculate_minimal_external(
                exit_mach, exit_radius, exit_angle, num_streamlines
            )
    
    def _determine_expansion_type(self, p_exit, p_ambient):
        """Determine if flow is over/under/perfectly expanded"""
        pressure_ratio = p_exit / p_ambient
        
        if pressure_ratio > 1.05:
            return 'under_expanded'
        elif pressure_ratio < 0.95:
            return 'over_expanded'
        else:
            return 'perfectly_expanded'
    
    def _calculate_external_prandtl_meyer(self, M_exit, r_exit, theta_exit, 
                                         p_exit, p_ambient, num_streamlines):
        """Calculate external Prandtl-Meyer expansion"""
        # Expansion fan from exit
        streamlines = []
        
        # Prandtl-Meyer function
        nu_exit = self._prandtl_meyer_angle(M_exit)
        
        # Additional turning needed
        pressure_ratio = p_exit / p_ambient
        # From isentropic relation
        additional_mach = M_exit * 1.2  # Simplified
        nu_additional = self._prandtl_meyer_angle(additional_mach) - nu_exit
        
        # Generate external streamlines
        for i in range(num_streamlines):
            fraction = float(i) / (num_streamlines - 1)
            r_start = r_exit * (1.0 - 0.5 * fraction)
            
            # Streamline points
            x_points = []
            y_points = []
            
            # Start at exit
            x_points.append(0.0)
            y_points.append(r_start)
            
            # Expansion region
            for j in range(10):
                x = j * 0.5
                turn = theta_exit + nu_additional * (j / 10.0) * fraction
                dr = x * np.tan(turn * CONV_DEG_RAD)
                
                x_points.append(x)
                y_points.append(r_start + dr)
            
            streamlines.append({
                'x': x_points,
                'y': y_points,
                'type': 'expansion'
            })
        
        return {
            'streamlines': streamlines,
            'expansion_type': 'under_expanded',
            'expansion_angle': nu_additional
        }
    
    def _calculate_external_compression(self, M_exit, r_exit, theta_exit, 
                                       p_exit, p_ambient, num_streamlines):
        """Calculate external compression (oblique shocks)"""
        # External oblique shock
        streamlines = []
        
        for i in range(num_streamlines):
            fraction = float(i) / (num_streamlines - 1)
            r_start = r_exit * (1.0 - 0.5 * fraction)
            
            x_points = [0.0]
            y_points = [r_start]
            
            # Compression/shock region
            for j in range(10):
                x = j * 0.3
                # Compression turns inward
                turn = theta_exit * (1.0 - j / 10.0)
                dr = x * np.tan(turn * CONV_DEG_RAD)
                
                x_points.append(x)
                y_points.append(r_start + dr * 0.5)  # Gradual compression
            
            streamlines.append({
                'x': x_points,
                'y': y_points,
                'type': 'compression'
            })
        
        return {
            'streamlines': streamlines,
            'expansion_type': 'over_expanded',
            'shock_present': True
        }
    
    def _calculate_minimal_external(self, M_exit, r_exit, theta_exit, num_streamlines):
        """Minimal external flow for perfectly expanded"""
        streamlines = []
        
        for i in range(num_streamlines):
            fraction = float(i) / (num_streamlines - 1)
            r_start = r_exit * (1.0 - 0.5 * fraction)
            
            # Nearly straight lines
            x_points = [0.0, 5.0]
            y_points = [r_start, r_start + 5.0 * np.tan(theta_exit * CONV_DEG_RAD * 0.5)]
            
            streamlines.append({
                'x': x_points,
                'y': y_points,
                'type': 'parallel'
            })
        
        return {
            'streamlines': streamlines,
            'expansion_type': 'perfectly_expanded'
        }
    
    def calculate_cowl_design(self, nozzle_length, exit_radius, exit_angle):
        """
        Calculate cowl (external shroud) design
        
        Args:
            nozzle_length: Nozzle length
            exit_radius: Exit radius
            exit_angle: Exit angle
        
        Returns:
            dict: Cowl geometry
        """
        # Cowl starts at throat region
        cowl_start_x = -nozzle_length * 0.1
        cowl_start_r = exit_radius * 1.5
        
        # Cowl extends beyond exit
        cowl_end_x = nozzle_length * 1.2
        cowl_end_r = exit_radius * 1.3
        
        # Simple linear cowl
        x_cowl = [cowl_start_x, 0.0, nozzle_length, cowl_end_x]
        y_cowl = [cowl_start_r, cowl_start_r, exit_radius * 1.2, cowl_end_r]
        
        return {
            'x': x_cowl,
            'y': y_cowl,
            'length': cowl_end_x - cowl_start_x
        }
    
    def calculate_afterbody_pressure(self, nozzle_results, p_ambient, afterbody_length=5.0):
        """
        Calculate pressure distribution on afterbody
        
        Args:
            nozzle_results: Nozzle MOC results
            p_ambient: Ambient pressure
            afterbody_length: Afterbody length (inches)
        
        Returns:
            dict: Afterbody pressure distribution
        """
        exit_x = np.max(nozzle_results['wall_x'])
        exit_r = np.max(nozzle_results['wall_y'])
        exit_p = nozzle_results.get('exit_p_pt', 0.5) * 50.0  # Estimate
        
        # Afterbody coordinates
        x_afterbody = np.linspace(exit_x, exit_x + afterbody_length, 20)
        
        # Pressure recovery along afterbody
        pressure_ratio = []
        for i, x in enumerate(x_afterbody):
            # Gradual pressure recovery
            progress = i / len(x_afterbody)
            p_local = exit_p + (p_ambient - exit_p) * progress ** 0.5
            pressure_ratio.append(p_local / p_ambient)
        
        return {
            'x': x_afterbody.tolist(),
            'pressure_ratio': pressure_ratio,
            'afterbody_length': afterbody_length
        }
    
    def calculate_drag(self, nozzle_results, p_ambient, afterbody_results):
        """
        Calculate external drag
        
        Returns:
            dict: Drag components
        """
        # Simplified drag calculation
        exit_area = np.pi * (np.max(nozzle_results['wall_y']) / 12.0) ** 2  # ft²
        exit_p = nozzle_results.get('exit_p_pt', 0.5) * 50.0 * 144.0  # psf
        
        # Pressure drag on afterbody
        pressure_drag = (p_ambient * 144.0 - exit_p) * exit_area * 0.3
        
        return {
            'pressure_drag': abs(pressure_drag),
            'total_drag': abs(pressure_drag)
        }
    
    def _prandtl_meyer_angle(self, mach):
        """Calculate Prandtl-Meyer angle"""
        if mach <= 1.0:
            return 0.0
        
        term1 = np.sqrt((self.gamma_p1) / (self.gamma_m1))
        term2 = np.arctan(np.sqrt((self.gamma_m1 / self.gamma_p1) * (mach**2 - 1.0)))
        term3 = np.arctan(np.sqrt(mach**2 - 1.0))
        
        nu = term1 * term2 - term3
        return nu * CONV_RAD_DEG