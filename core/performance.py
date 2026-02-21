"""
Performance Calculator Module
Complete nozzle performance analysis
"""

import numpy as np
from config import R_GAS, GC


class PerformanceCalculator:
    """Advanced nozzle performance calculator"""
    
    def __init__(self, gamma=1.4):
        self.gamma = gamma
        self.gamma_m1 = gamma - 1.0
        self.gamma_p1 = gamma + 1.0
        self.R = R_GAS
        self.gc = GC
    
    def calculate_all_performance(self, moc_results, p_total, T_total, p_ambient, 
                                  throat_area, exit_area):
        """
        Calculate complete performance metrics
        
        Args:
            moc_results: MOC calculation results
            p_total: Total pressure (lbf/ft²)
            T_total: Total temperature (Rankine)
            p_ambient: Ambient pressure (lbf/ft²)
            throat_area: Throat area (ft²)
            exit_area: Exit area (ft²)
            
        Returns:
            dict: Complete performance data
        """
        # Exit conditions from MOC
        M_exit = moc_results['exit_mach']
        p_ratio_exit = moc_results['exit_p_pt']
        T_ratio_exit = moc_results['exit_T_Tt']
        
        # Static conditions at exit
        p_exit = p_total * p_ratio_exit
        T_exit = T_total * T_ratio_exit
        
        # Mass flow rate (choked flow through throat)
        mass_flow = self.calculate_mass_flow_choked(p_total, T_total, throat_area)
        
        # Exit velocity
        V_exit = self.calculate_exit_velocity(M_exit, T_exit)
        
        # Thrust components
        momentum_thrust = (mass_flow / self.gc) * V_exit
        pressure_thrust = (p_exit - p_ambient) * exit_area
        total_thrust = momentum_thrust + pressure_thrust
        
        # Nozzle Pressure Ratio
        NPR = p_total / p_ambient
        
        # Ideal thrust (perfectly expanded)
        p_exit_ideal = p_ambient
        # Calculate ideal exit Mach for perfectly expanded
        M_exit_ideal = self._get_mach_for_pressure_ratio(p_exit_ideal / p_total)
        T_exit_ideal = T_total / (1.0 + 0.5 * self.gamma_m1 * M_exit_ideal**2)
        V_exit_ideal = self.calculate_exit_velocity(M_exit_ideal, T_exit_ideal)
        ideal_thrust = (mass_flow / self.gc) * V_exit_ideal
        
        # Nozzle efficiency
        efficiency = (total_thrust / ideal_thrust) * 100.0 if ideal_thrust > 0 else 0.0
        
        # Specific impulse
        Isp = total_thrust / (mass_flow / self.gc)
        
        # Characteristic velocity (C*)
        C_star = self._calculate_characteristic_velocity(T_total)
        
        # Thrust coefficient
        C_F = total_thrust / (p_total * throat_area)
        
        # Expansion ratio
        expansion_ratio = exit_area / throat_area
        
        # Determine expansion condition
        expansion_condition = self._determine_expansion_condition(p_exit, p_ambient)
        
        return {
            'thrust': total_thrust,
            'momentum_thrust': momentum_thrust,
            'pressure_thrust': pressure_thrust,
            'mass_flow': mass_flow,
            'NPR': NPR,
            'Isp': Isp,
            'efficiency': efficiency,
            'C_star': C_star,
            'C_F': C_F,
            'expansion_ratio': expansion_ratio,
            'exit_velocity': V_exit,
            'exit_pressure': p_exit,
            'exit_temperature': T_exit,
            'expansion_condition': expansion_condition
        }
    
    def calculate_mass_flow_choked(self, p_total, T_total, throat_area):
        """
        Calculate mass flow rate for choked flow (M=1 at throat)
        
        ṁ = (p_t * A*) / √(T_t) * √(γ/R) * [(2/(γ+1))]^[(γ+1)/(2(γ-1))]
        """
        factor = np.sqrt(self.gamma / self.R)
        exponent = self.gamma_p1 / (2.0 * self.gamma_m1)
        flow_factor = np.power(2.0 / self.gamma_p1, exponent)
        
        mass_flow = (p_total * throat_area * factor * flow_factor) / np.sqrt(T_total)
        
        return mass_flow
    
    def calculate_exit_velocity(self, mach, temperature):
        """
        Calculate exit velocity: V = M * √(γ * R * T)
        """
        a = np.sqrt(self.gamma * self.R * temperature)
        velocity = mach * a
        return velocity
    
    def _calculate_characteristic_velocity(self, T_total):
        """
        Calculate characteristic velocity (C*)
        
        C* = √(γ * R * T_t) / (γ * √[(2/(γ+1))^((γ+1)/(γ-1))])
        """
        numerator = np.sqrt(self.gamma * self.R * T_total)
        exponent = self.gamma_p1 / self.gamma_m1
        denominator = self.gamma * np.sqrt(np.power(2.0 / self.gamma_p1, exponent))
        
        C_star = numerator / denominator
        return C_star
    
    def _get_mach_for_pressure_ratio(self, p_p0):
        """
        Get Mach number from pressure ratio using isentropic relation
        """
        exponent = self.gamma_m1 / self.gamma
        mach_squared = (2.0 / self.gamma_m1) * (np.power(p_p0, -exponent) - 1.0)
        
        if mach_squared < 0:
            return 1.0
        
        return np.sqrt(mach_squared)
    
    def _determine_expansion_condition(self, p_exit, p_ambient):
        """Determine if nozzle is over/under/perfectly expanded"""
        ratio = p_exit / p_ambient
        
        if ratio > 1.05:
            return "Under-Expanded"
        elif ratio < 0.95:
            return "Over-Expanded"
        else:
            return "Perfectly Expanded"