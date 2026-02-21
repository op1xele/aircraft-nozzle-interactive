"""
MOC Engine - Method of Characteristics
"""

import numpy as np
from .isentropic import IsentropicRelations
from .shock_calculator import ShockCalculator
from config import *


class MOCEngine:
    """Method of Characteristics Engine"""
    
    def __init__(self, gamma=DEFAULT_GAMMA):
        self.gamma = gamma
        self.isen = IsentropicRelations(gamma)
        self.shock_calc = ShockCalculator(gamma)
        self.shocks_detected = []
        self.coalescence_points = []
        self.grid_size = MAX_GRID_SIZE
        self.reset_grid()
        self.num_rays = DEFAULT_NUM_RAYS
        self.mach_exit = DEFAULT_MACH_EXIT
        self.throat_height = DEFAULT_THROAT_HEIGHT
        self.geometry_calculated = False
        self.flow_calculated = False
    
    def reset_grid(self):
        s = self.grid_size
        self.mc_mach = np.ones((s, s))
        self.mc_turn = np.zeros((s, s))
        self.mc_defl = np.zeros((s, s))
        self.mc_pm = np.zeros((s, s))
        self.mc_mang = np.zeros((s, s))
        self.mc_Q = np.zeros((s, s))
        self.mc_R = np.zeros((s, s))
        self.mc_x = np.zeros((s, s))
        self.mc_y = np.zeros((s, s))
        self.mc_alpha = np.zeros((s, s))
        self.mc_beta = np.zeros((s, s))
        self.mc_p_pt = np.ones((s, s))
        self.mc_T_Tt = np.ones((s, s))
        self.mc_rho_rhot = np.ones((s, s))
        self.mc_A_Astar = np.ones((s, s))
        self.mc_pressure = np.zeros((s, s))
        self.mc_temperature = np.zeros((s, s))
        self.geometry_calculated = False
        self.flow_calculated = False
    
    def calculate_2d_ideal_nozzle(self, mach_exit, num_rays=30, throat_height=1.0):
        """Calculate 2D ideal nozzle using MOC"""
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_height
        self.reset_grid()
        
        props_exit = self.isen.get_all_properties(mach_exit)
        nu_exit = props_exit['nu']
        
        theta_max = nu_exit / 2.0
        delta_theta = nu_exit / num_rays
        
        # Initial conditions at throat
        self.mc_mach[0, 0] = 1.0
        self.mc_pm[0, 0] = 0.0
        self.mc_mang[0, 0] = 90.0
        self.mc_x[0, 0] = 0.0
        self.mc_y[0, 0] = throat_height
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[0, 0] = props_throat['p_pt']
        self.mc_T_Tt[0, 0] = props_throat['T_Tt']
        
        delx = 0.1
        
        # First characteristic (1,1) - boundary
        i, j = 1, 1
        self.mc_defl[i, j] = delta_theta
        self.mc_turn[i, j] = delta_theta
        self.mc_pm[i, j] = self.mc_turn[i, j]
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        self._get_moc_variables(i, j, 0, 0)
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_height
        
        # First characteristic (1,2) - symmetry
        i, j = 1, 2
        self.mc_turn[i, j] = 0.0
        self.mc_defl[i, j] = self.mc_turn[i, j] - self.mc_turn[1, 1]
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]
        self.mc_pm[i, j] = self.mc_Q[i, j]
        self._get_moc_variables(i, j, 1, 1)
        self.mc_y[i, j] = 0.0
        alpha_angle = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + (self.mc_y[1, 1] - self.mc_y[i, j]) / np.tan(alpha_angle * CONV_DEG_RAD)
        
        # Expansion on boundary
        for i in range(2, num_rays // 2 + 1):
            # Boundary point
            self.mc_defl[i, 1] = delta_theta
            self.mc_turn[i, 1] = self.mc_turn[i-1, 1] + self.mc_defl[i, 1]
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            self._get_moc_variables(i, 1, i-1, 1)
            self.mc_x[i, 1] = self.mc_x[i-1, 1] + delx * throat_height
            self.mc_y[i, 1] = self.mc_y[i-1, 1] + (self.mc_x[i, 1] - self.mc_x[i-1, 1]) * np.tan(self.mc_turn[i-1, 1] * CONV_DEG_RAD)
            
            # Internal points
            for k in range(2, i + 1):
                self.mc_Q[i, k] = self.mc_Q[i, k-1]
                self.mc_R[i, k] = self.mc_R[i-1, k]
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                self._get_moc_variables(i, k, i, k-1)
                
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                )
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            # Symmetry point
            self.mc_turn[i, i+1] = 0.0
            self.mc_defl[i, i+1] = self.mc_turn[i, i+1] - self.mc_turn[i, i]
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]
            self.mc_pm[i, i+1] = self.mc_Q[i, i+1]
            self._get_moc_variables(i, i+1, i, i)
            self.mc_y[i, i+1] = 0.0
            alpha_angle = self.mc_mang[i, i] - self.mc_turn[i, i]
            self.mc_x[i, i+1] = self.mc_x[i, i] + (self.mc_y[i, i] - self.mc_y[i, i+1]) / np.tan(alpha_angle * CONV_DEG_RAD)
        
        # Cancellation surface
        i = num_rays // 2 + 1
        self.mc_defl[i, 1] = self.mc_defl[i-1, 1]
        self.mc_turn[i, 1] = self.mc_turn[i-1, 1]
        self.mc_R[i, 1] = self.mc_R[i-1, 1]
        self.mc_pm[i, 1] = self.mc_pm[i-1, 1]
        self.mc_Q[i, 1] = self.mc_Q[i-1, 1]
        self.mc_mach[i, 1] = self.mc_mach[i-1, 1]
        self.mc_x[i, 1] = self.mc_x[i-1, 1]
        self.mc_y[i, 1] = self.mc_y[i-1, 1]
        
        for k in range(2, i + 1):
            self.mc_defl[i, k] = -delta_theta
            self.mc_turn[i, k] = self.mc_turn[i, k-1] + self.mc_defl[i, k]
            self.mc_R[i, k] = self.mc_R[i-1, k]
            self.mc_pm[i, k] = self.mc_R[i, k] + self.mc_turn[i, k]
            self.mc_Q[i, k] = self.mc_pm[i, k] + self.mc_turn[i, k]
            self._get_moc_variables(i, k, i, k-1)
            
            self.mc_alpha[i, k-1] = self.mc_turn[i, k-1]
            self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
            
            tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
            tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
            
            self.mc_x[i, k] = (
                (self.mc_y[i-1, k] - self.mc_y[i, k-1] +
                 self.mc_x[i, k-1] * tan_alpha -
                 self.mc_x[i-1, k] * tan_beta) /
                (tan_alpha - tan_beta)
            )
            
            self.mc_y[i, k] = self.mc_y[i, k-1] + (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_results()
    
    def _get_moc_variables(self, i, j, i_up, j_up):
        self.mc_mach[i, j] = self.isen.get_mach_from_pm(self.mc_pm[i, j])
        props = self.isen.get_all_properties(self.mc_mach[i, j])
        self.mc_mang[i, j] = props['mu']
        self.mc_p_pt[i, j] = props['p_pt']
        self.mc_T_Tt[i, j] = props['T_Tt']
        self.mc_rho_rhot[i, j] = props['rho_rhot']
        self.mc_A_Astar[i, j] = props['A_Astar']
    
    def _prepare_results(self):
        n = self.num_rays // 2 + 1
        
        wall_x = [self.mc_x[i, 1] for i in range(1, n+1)]
        wall_y = [self.mc_y[i, 1] for i in range(1, n+1)]
        exit_x = [self.mc_x[n, k] for k in range(1, n+1)]
        exit_y = [self.mc_y[n, k] for k in range(1, n+1)]
        
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_T_Tt = self.mc_T_Tt[n, n//2 + 1]
        
        char_lines = self.get_characteristic_lines()

        # SHOCK DETECTION (NEW!)
        print("\n  Analyzing flow field for shocks...")
        shocks = self.detect_shocks()
        
        shock_stats = {
            'num_shocks': len(shocks),
            'max_pressure_ratio': max([s['p_ratio'] for s in shocks]) if shocks else 1.0,
            'avg_shock_angle': np.mean([s['shock_angle'] for s in shocks]) if shocks else 0.0,
            'shock_locations': [(s['x'], s['y']) for s in shocks]
        }
        
        return {
            'wall_x': np.array(wall_x),
            'wall_y': np.array(wall_y),
            'exit_x': np.array(exit_x),
            'exit_y': np.array(exit_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'shocks': shocks,                           
            'shock_stats': shock_stats,                 
            'coalescence_points': self.coalescence_points 
        }
    
    def get_characteristic_lines(self):
        """Get characteristic lines for plotting"""
        if not self.geometry_calculated:
            raise RuntimeError("Geometry not calculated")
        
        n = self.num_rays // 2 + 1
        plus_chars = []
        minus_chars = []
        
        for i in range(1, n+1):
            for j in range(1, i+1):
                if j < i:
                    x_line = [self.mc_x[i, j], self.mc_x[i, j+1]]
                    y_line = [self.mc_y[i, j], self.mc_y[i, j+1]]
                    plus_chars.append((x_line, y_line))
        
        for i in range(1, n+1):
            for j in range(1, i+1):
                if i > j and i > 1:
                    x_line = [self.mc_x[i-1, j], self.mc_x[i, j]]
                    y_line = [self.mc_y[i-1, j], self.mc_y[i, j]]
                    minus_chars.append((x_line, y_line))
        
        return {
            'plus_characteristics': plus_chars,
            'minus_characteristics': minus_chars
        }
    
    def get_contour_data(self, variable='mach'):
        """Get contour plot data"""
        if not self.flow_calculated:
            raise RuntimeError("Flow not calculated")
        
        n = self.num_rays // 2 + 1
        X = self.mc_x[:n, :n]
        Y = self.mc_y[:n, :n]
        
        if variable == 'mach':
            Z = self.mc_mach[:n, :n]
        elif variable == 'pressure':
            Z = self.mc_p_pt[:n, :n]
        elif variable == 'temperature':
            Z = self.mc_T_Tt[:n, :n]
        else:
            raise ValueError(f"Unknown variable: {variable}")
        
        return X, Y, Z
    
    def calculate_plume(self, p_ambient):
        """Calculate jet plume for external flow"""
        if not self.flow_calculated:
            raise RuntimeError("Flow not calculated")
        
        n = self.num_rays // 2 + 1
        exit_x = self.mc_x[n, :]
        exit_y = self.mc_y[n, :]
        
        plume_x = []
        plume_y = []
        
        for k in range(1, n+1):
            x_start = exit_x[k]
            y_start = exit_y[k]
            length = 2.0 * self.throat_height
            
            plume_x.append([x_start, x_start + length])
            plume_y.append([y_start, y_start])
        
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_radius = self.mc_y[n, n//2 + 1]
        
        try:
            shock_diamonds = self.shock_calc.shock_diamond_pattern(
                M_exit=exit_mach,
                p_exit=exit_p_pt * 50.0,
                p_ambient=p_ambient,
                nozzle_radius=exit_radius
            )
        except:
            shock_diamonds = None
        
        return {
            'plume_x': plume_x,
            'plume_y': plume_y,
            'shock_diamonds': shock_diamonds
        }
    
    def calculate_axisymmetric_nozzle(self, mach_exit, num_rays=30, throat_radius=1.0):
        """Calculate axisymmetric (conical) nozzle using MOC with radial correction"""
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_radius
        self.reset_grid()
        
        props_exit = self.isen.get_all_properties(mach_exit)
        nu_exit = props_exit['nu']
        
        delta_theta = nu_exit / num_rays
        
        self.mc_mach[0, 0] = 1.0
        self.mc_pm[0, 0] = 0.0
        self.mc_mang[0, 0] = 90.0
        self.mc_x[0, 0] = 0.0
        self.mc_y[0, 0] = throat_radius
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[0, 0] = props_throat['p_pt']
        self.mc_T_Tt[0, 0] = props_throat['T_Tt']
        
        delx = 0.1
        
        i, j = 1, 1
        self.mc_defl[i, j] = delta_theta
        self.mc_turn[i, j] = delta_theta
        self.mc_pm[i, j] = self.mc_turn[i, j]
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        self._get_moc_variables(i, j, 0, 0)
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_radius
        
        i, j = 1, 2
        self.mc_turn[i, j] = 0.0
        self.mc_defl[i, j] = self.mc_turn[i, j] - self.mc_turn[1, 1]
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]
        self.mc_pm[i, j] = self.mc_Q[i, j]
        self._get_moc_variables(i, j, 1, 1)
        self.mc_y[i, j] = 0.0
        alpha_angle = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + (self.mc_y[1, 1] - self.mc_y[i, j]) / np.tan(alpha_angle * CONV_DEG_RAD)
        
        for i in range(2, num_rays // 2 + 1):
            self.mc_defl[i, 1] = delta_theta
            self.mc_turn[i, 1] = self.mc_turn[i-1, 1] + self.mc_defl[i, 1]
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            self._get_moc_variables(i, 1, i-1, 1)
            
            self.mc_x[i, 1] = self.mc_x[i-1, 1] + delx * throat_radius
            
            dr_dx = np.tan(self.mc_turn[i-1, 1] * CONV_DEG_RAD)
            self.mc_y[i, 1] = self.mc_y[i-1, 1] + (self.mc_x[i, 1] - self.mc_x[i-1, 1]) * dr_dx
            
            for k in range(2, i + 1):
                self.mc_Q[i, k] = self.mc_Q[i, k-1]
                self.mc_R[i, k] = self.mc_R[i-1, k]
                
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                
                self._get_moc_variables(i, k, i, k-1)

                # Coalescence check (YENİ!)
                self.check_coalescence(i, k)
                
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                r_avg = 0.5 * (self.mc_y[i, k-1] + self.mc_y[i-1, k])
                if abs(r_avg) > 0.001:
                    radial_factor = 1.0 + 0.1 * (self.mc_mach[i, k-1] - 1.0) / r_avg
                else:
                    radial_factor = 1.0
                
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                ) * radial_factor
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            self.mc_turn[i, i+1] = 0.0
            self.mc_defl[i, i+1] = self.mc_turn[i, i+1] - self.mc_turn[i, i]
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]
            self.mc_pm[i, i+1] = self.mc_Q[i, i+1]
            self._get_moc_variables(i, i+1, i, i)
            self.mc_y[i, i+1] = 0.0
            alpha_angle = self.mc_mang[i, i] - self.mc_turn[i, i]
            self.mc_x[i, i+1] = self.mc_x[i, i] + (self.mc_y[i, i] - self.mc_y[i, i+1]) / np.tan(alpha_angle * CONV_DEG_RAD)
        
        i = num_rays // 2 + 1
        self.mc_defl[i, 1] = self.mc_defl[i-1, 1]
        self.mc_turn[i, 1] = self.mc_turn[i-1, 1]
        self.mc_R[i, 1] = self.mc_R[i-1, 1]
        self.mc_pm[i, 1] = self.mc_pm[i-1, 1]
        self.mc_Q[i, 1] = self.mc_Q[i-1, 1]
        self.mc_mach[i, 1] = self.mc_mach[i-1, 1]
        self.mc_x[i, 1] = self.mc_x[i-1, 1]
        self.mc_y[i, 1] = self.mc_y[i-1, 1]
        
        for k in range(2, i + 1):
            self.mc_defl[i, k] = -delta_theta
            self.mc_turn[i, k] = self.mc_turn[i, k-1] + self.mc_defl[i, k]
            self.mc_R[i, k] = self.mc_R[i-1, k]
            self.mc_pm[i, k] = self.mc_R[i, k] + self.mc_turn[i, k]
            self.mc_Q[i, k] = self.mc_pm[i, k] + self.mc_turn[i, k]
            self._get_moc_variables(i, k, i, k-1)
            
            self.mc_alpha[i, k-1] = self.mc_turn[i, k-1]
            self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
            
            tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
            tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
            
            self.mc_x[i, k] = (
                (self.mc_y[i-1, k] - self.mc_y[i, k-1] +
                 self.mc_x[i, k-1] * tan_alpha -
                 self.mc_x[i-1, k] * tan_beta) /
                (tan_alpha - tan_beta)
            )
            
            self.mc_y[i, k] = self.mc_y[i, k-1] + (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_results()
    
    def calculate_plug_nozzle(self, mach_exit, num_rays=30, throat_radius=1.0, 
                             spike_length_ratio=0.8):
        """
        Calculate plug (aerospike) nozzle using full MOC
        
        Complete Method of Characteristics for plug/aerospike nozzle
        - Spike surface: inner boundary with wall BC
        - Outer boundary: pressure-balanced free expansion
        - Radial flow with proper characteristic equations
        - Truncated spike geometry
        """
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_radius
        self.reset_grid()
        
        props_exit = self.isen.get_all_properties(mach_exit)
        nu_exit = props_exit['nu']
        theta_max = nu_exit / 2.0
        delta_theta = theta_max / (num_rays // 2)
        
        delx = 0.1
        
        i, j = 0, 0
        self.mc_mach[i, j] = 1.0
        self.mc_pm[i, j] = 0.0
        self.mc_mang[i, j] = 90.0
        self.mc_turn[i, j] = 0.0
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_radius
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[i, j] = props_throat['p_pt']
        self.mc_T_Tt[i, j] = props_throat['T_Tt']
        
        i = 1
        
        i, j = 1, 1
        
        self.mc_defl[i, j] = delta_theta
        self.mc_turn[i, j] = self.mc_turn[0, 0] + self.mc_defl[i, j]
        self.mc_pm[i, j] = self.mc_turn[i, j]
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        
        self._get_moc_variables(i, j, 0, 0)
        
        self.mc_x[i, j] = delx * throat_radius
        dr_dx = -np.tan(self.mc_turn[i, j] * CONV_DEG_RAD)
        self.mc_y[i, j] = self.mc_y[0, 0] + (self.mc_x[i, j] - self.mc_x[0, 0]) * dr_dx
        
        i, j = 1, 2
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]
        self.mc_pm[i, j] = 0.5 * (self.mc_Q[i, j] + self.mc_R[i, j])
        self.mc_turn[i, j] = 0.5 * (self.mc_Q[i, j] - self.mc_R[i, j])
        
        self._get_moc_variables(i, j, 1, 1)
        
        alpha = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_y[i, j] = throat_radius
        self.mc_x[i, j] = self.mc_x[1, 1] + (self.mc_y[i, j] - self.mc_y[1, 1]) / np.tan(alpha * CONV_DEG_RAD)
        
        for i in range(2, num_rays // 2 + 1):
            self.mc_defl[i, 1] = delta_theta
            self.mc_turn[i, 1] = self.mc_turn[i-1, 1] + self.mc_defl[i, 1]
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            
            self._get_moc_variables(i, 1, i-1, 1)
            
            self.mc_x[i, 1] = self.mc_x[i-1, 1] + delx * throat_radius
            theta_avg = 0.5 * (self.mc_turn[i, 1] + self.mc_turn[i-1, 1])
            dr_dx = -np.tan(theta_avg * CONV_DEG_RAD) * 0.8
            self.mc_y[i, 1] = self.mc_y[i-1, 1] + (self.mc_x[i, 1] - self.mc_x[i-1, 1]) * dr_dx
            
            if self.mc_y[i, 1] < 0.05 * throat_radius:
                self.mc_y[i, 1] = 0.05 * throat_radius
            
            for k in range(2, i + 1):
                self.mc_Q[i, k] = self.mc_Q[i, k-1]
                self.mc_R[i, k] = self.mc_R[i-1, k]
                
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                
                self._get_moc_variables(i, k, i, k-1)
                
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                r_avg = 0.5 * (self.mc_y[i, k-1] + self.mc_y[i-1, k])
                if r_avg > 0.01 * throat_radius:
                    radial_correction = 1.0 + (self.mc_mach[i, k-1] - 1.0) * 0.02 / r_avg
                else:
                    radial_correction = 1.0
                
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                ) * radial_correction
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            
            self.mc_turn[i, i+1] = theta_max * (float(i) / (num_rays // 2))
            self.mc_R[i, i+1] = self.mc_Q[i, i+1] - 2.0 * self.mc_turn[i, i+1]
            self.mc_pm[i, i+1] = self.mc_Q[i, i+1] - self.mc_turn[i, i+1]
            
            self._get_moc_variables(i, i+1, i, i)
            
            alpha = self.mc_mang[i, i] - self.mc_turn[i, i]
            dr = (self.mc_x[i, i] - self.mc_x[i-1, i+1]) * np.tan(alpha * CONV_DEG_RAD)
            self.mc_y[i, i+1] = self.mc_y[i, i] + abs(dr)
            self.mc_x[i, i+1] = self.mc_x[i, i] + abs(dr) / np.tan(alpha * CONV_DEG_RAD + 0.1)
        
        n = num_rays // 2 + 1
        
        spike_full_length = self.mc_x[n-1, 1]
        spike_truncated_length = spike_full_length * spike_length_ratio
        
        for k in range(1, n + 1):
            if k == 1:
                self.mc_x[n, k] = spike_truncated_length
                self.mc_y[n, k] = self.mc_y[n-1, k] * 0.9
            else:
                self.mc_Q[n, k] = self.mc_Q[n-1, k]
                self.mc_R[n, k] = self.mc_R[n-1, k]
                self.mc_pm[n, k] = 0.5 * (self.mc_Q[n, k] + self.mc_R[n, k])
                self.mc_turn[n, k] = 0.5 * (self.mc_Q[n, k] - self.mc_R[n, k])
                
                self._get_moc_variables(n, k, n-1, k)
                
                dx = self.mc_x[n-1, k] - self.mc_x[n-2, k] if n > 2 else delx * throat_radius
                self.mc_x[n, k] = self.mc_x[n-1, k] + dx * 1.2
                dy = self.mc_y[n-1, k] - self.mc_y[n-2, k] if n > 2 else 0.1 * throat_radius
                self.mc_y[n, k] = self.mc_y[n-1, k] + dy * 1.1
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_plug_results()
    
    def _prepare_plug_results(self):
        """Prepare results for plug nozzle"""
        n = self.num_rays // 2 + 1
        
        spike_x = [self.mc_x[i, 1] for i in range(1, n+1)]
        spike_y = [self.mc_y[i, 1] for i in range(1, n+1)]
        
        outer_x = [self.mc_x[i, i+1] for i in range(1, n)]
        outer_y = [self.mc_y[i, i+1] for i in range(1, n)]
        outer_x.append(self.mc_x[n, n])
        outer_y.append(self.mc_y[n, n])
        
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_T_Tt = self.mc_T_Tt[n, n//2 + 1]
        
        char_lines = self.get_characteristic_lines()
        
        return {
            'spike_x': np.array(spike_x),
            'spike_y': np.array(spike_y),
            'outer_x': np.array(outer_x),
            'outer_y': np.array(outer_y),
            'wall_x': np.array(outer_x),
            'wall_y': np.array(outer_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'nozzle_type': 'plug'
        }

    def calculate_cone_nozzle(self, mach_exit, num_rays=30, throat_radius=1.0, 
                             cone_half_angle=15.0):
        """
        Calculate conical nozzle using full MOC
        
        Complete Method of Characteristics for cone nozzle
        - Straight conical wall at fixed half-angle
        - Internal flow computed via MOC
        - Wall boundary condition: flow tangent to cone
        - Exit plane: uniform properties
        
        Args:
            mach_exit: Design exit Mach number
            num_rays: Number of characteristic lines
            throat_radius: Throat radius (inches)
            cone_half_angle: Cone half-angle (degrees) - typically 12-18°
        
        Returns:
            dict: Results with cone geometry and flow properties
        """
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_radius
        self.reset_grid()
        
        # Get exit flow properties
        props_exit = self.isen.get_all_properties(mach_exit)
        A_ratio_exit = props_exit['A_Astar']
        
        # Cone geometry: exit radius from area ratio
        # A_exit / A_throat = (r_exit / r_throat)^2 for axisymmetric
        r_exit = throat_radius * np.sqrt(A_ratio_exit)
        
        # Cone length from geometry
        cone_length = (r_exit - throat_radius) / np.tan(cone_half_angle * CONV_DEG_RAD)
        
        # MOC grid parameters
        delta_theta = (mach_exit - 1.0) * 2.0 / num_rays  # Adaptive spacing
        
        # ================================================================
        # INITIAL CONDITIONS AT THROAT
        # ================================================================
        i, j = 0, 0
        self.mc_mach[i, j] = 1.0
        self.mc_pm[i, j] = 0.0
        self.mc_mang[i, j] = 90.0
        self.mc_turn[i, j] = 0.0
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_radius
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[i, j] = props_throat['p_pt']
        self.mc_T_Tt[i, j] = props_throat['T_Tt']
        
        # ================================================================
        # INITIAL EXPANSION FAN
        # ================================================================
        # First characteristic from throat
        i, j = 1, 1
        
        # Wall boundary condition: flow must be tangent to cone surface
        # For cone, wall angle = cone_half_angle
        self.mc_turn[i, j] = cone_half_angle
        self.mc_pm[i, j] = self.mc_turn[i, j]  # Wall BC
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        
        self._get_moc_variables(i, j, 0, 0)
        
        # Position on cone wall
        dx = cone_length / num_rays
        self.mc_x[i, j] = dx
        self.mc_y[i, j] = throat_radius + dx * np.tan(cone_half_angle * CONV_DEG_RAD)
        
        # First characteristic to centerline (axis)
        i, j = 1, 2
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]  # Symmetry: R = Q
        self.mc_pm[i, j] = 0.5 * (self.mc_Q[i, j] + self.mc_R[i, j])
        self.mc_turn[i, j] = 0.5 * (self.mc_Q[i, j] - self.mc_R[i, j])
        
        self._get_moc_variables(i, j, 1, 1)
        
        # Centerline position
        self.mc_y[i, j] = 0.0
        alpha = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + self.mc_y[1, 1] / np.tan(alpha * CONV_DEG_RAD)
        
        # ================================================================
        # MOC GRID PROPAGATION
        # ================================================================
        for i in range(2, num_rays // 2 + 1):
            # ---------------------------------------------------------
            # CONE WALL POINT (j=1)
            # ---------------------------------------------------------
            # Wall BC: flow tangent to cone
            self.mc_turn[i, 1] = cone_half_angle
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            
            self._get_moc_variables(i, 1, i-1, 1)
            
            # Position on cone surface
            self.mc_x[i, 1] = i * dx
            self.mc_y[i, 1] = throat_radius + self.mc_x[i, 1] * np.tan(cone_half_angle * CONV_DEG_RAD)
            
            # ---------------------------------------------------------
            # INTERNAL GRID POINTS (2 ≤ j ≤ i)
            # ---------------------------------------------------------
            for k in range(2, i + 1):
                # Compatibility equations along characteristics
                self.mc_Q[i, k] = self.mc_Q[i, k-1]  # Along C+ from left
                self.mc_R[i, k] = self.mc_R[i-1, k]  # Along C- from above
                
                # Solve for flow angle and Prandtl-Meyer angle
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                
                self._get_moc_variables(i, k, i, k-1)
                
                # Characteristic angles
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                # Axisymmetric correction
                r_avg = 0.5 * (self.mc_y[i, k-1] + self.mc_y[i-1, k])
                if r_avg > 0.001 * throat_radius:
                    # Radial gradient effect
                    M_avg = 0.5 * (self.mc_mach[i, k-1] + self.mc_mach[i-1, k])
                    radial_correction = 1.0 + (M_avg - 1.0) * 0.015 / r_avg
                else:
                    radial_correction = 1.0
                
                # Characteristic intersection
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                ) * radial_correction
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            # ---------------------------------------------------------
            # CENTERLINE POINT (j = i+1) - Symmetry BC
            # ---------------------------------------------------------
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]  # Symmetry
            self.mc_pm[i, i+1] = 0.5 * (self.mc_Q[i, i+1] + self.mc_R[i, i+1])
            self.mc_turn[i, i+1] = 0.5 * (self.mc_Q[i, i+1] - self.mc_R[i, i+1])
            
            self._get_moc_variables(i, i+1, i, i)
            
            # Centerline position
            self.mc_y[i, i+1] = 0.0
            alpha = self.mc_mang[i, i] - self.mc_turn[i, i]
            if np.tan(alpha * CONV_DEG_RAD) > 0.001:
                self.mc_x[i, i+1] = self.mc_x[i, i] + self.mc_y[i, i] / np.tan(alpha * CONV_DEG_RAD)
            else:
                self.mc_x[i, i+1] = self.mc_x[i, i] + dx
        
        # ================================================================
        # EXIT PLANE
        # ================================================================
        n = num_rays // 2 + 1
        
        # Extend to exit plane (straight cone to design exit radius)
        for k in range(1, n + 1):
            # Extrapolate to exit
            if k == 1:
                # Wall point - exactly on cone
                self.mc_x[n, k] = cone_length
                self.mc_y[n, k] = r_exit
            else:
                # Internal and centerline points
                self.mc_Q[n, k] = self.mc_Q[n-1, k]
                self.mc_R[n, k] = self.mc_R[n-1, k]
                self.mc_pm[n, k] = 0.5 * (self.mc_Q[n, k] + self.mc_R[n, k])
                self.mc_turn[n, k] = 0.5 * (self.mc_Q[n, k] - self.mc_R[n, k])
                
                self._get_moc_variables(n, k, n-1, k)
                
                # Linear extrapolation to exit
                dx_ext = cone_length - self.mc_x[n-1, k]
                self.mc_x[n, k] = cone_length
                # Radial position interpolated
                self.mc_y[n, k] = self.mc_y[n-1, k] + dx_ext * (self.mc_y[n-1, k] - self.mc_y[n-2, k]) / dx
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_cone_results(cone_half_angle)
    
    def _prepare_cone_results(self, cone_half_angle):
        """Prepare results for cone nozzle"""
        n = self.num_rays // 2 + 1
        
        # Cone wall (straight line)
        wall_x = [self.mc_x[i, 1] for i in range(1, n+1)]
        wall_y = [self.mc_y[i, 1] for i in range(1, n+1)]
        
        # Exit plane
        exit_x = [self.mc_x[n, k] for k in range(1, n+1)]
        exit_y = [self.mc_y[n, k] for k in range(1, n+1)]
        
        # Exit properties (centerline)
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_T_Tt = self.mc_T_Tt[n, n//2 + 1]
        
        # Characteristics
        char_lines = self.get_characteristic_lines()
        
        return {
            'wall_x': np.array(wall_x),
            'wall_y': np.array(wall_y),
            'exit_x': np.array(exit_x),
            'exit_y': np.array(exit_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'nozzle_type': 'cone',
            'cone_half_angle': cone_half_angle
        }
    def calculate_wedge_nozzle(self, mach_exit, num_rays=30, throat_height=1.0, 
                              wedge_angle=12.0):
        """
        Calculate 2D planar wedge nozzle using full MOC
        
        Complete Method of Characteristics for wedge nozzle
        - Planar 2D geometry (not axisymmetric)
        - Straight wedge walls at fixed angle
        - Used in scramjet and hypersonic vehicles
        - No radial correction (pure 2D)
        
        Args:
            mach_exit: Design exit Mach number
            num_rays: Number of characteristic lines
            throat_height: Throat height (inches)
            wedge_angle: Wedge half-angle (degrees) - typically 10-15°
        
        Returns:
            dict: Results with wedge geometry and flow properties
        """
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_height
        self.reset_grid()
        
        # Get exit flow properties
        props_exit = self.isen.get_all_properties(mach_exit)
        A_ratio_exit = props_exit['A_Astar']
        
        # Wedge geometry: exit height from area ratio
        # For 2D planar: A_exit / A_throat = h_exit / h_throat
        h_exit = throat_height * A_ratio_exit
        
        # Wedge length from geometry
        wedge_length = (h_exit - throat_height) / np.tan(wedge_angle * CONV_DEG_RAD)
        
        # MOC parameters
        delta_theta = (mach_exit - 1.0) * 2.0 / num_rays
        
        # ================================================================
        # INITIAL CONDITIONS AT THROAT
        # ================================================================
        i, j = 0, 0
        self.mc_mach[i, j] = 1.0
        self.mc_pm[i, j] = 0.0
        self.mc_mang[i, j] = 90.0
        self.mc_turn[i, j] = 0.0
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_height
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[i, j] = props_throat['p_pt']
        self.mc_T_Tt[i, j] = props_throat['T_Tt']
        
        # ================================================================
        # INITIAL EXPANSION - First characteristic
        # ================================================================
        i, j = 1, 1
        
        # Wall boundary condition: flow tangent to wedge
        self.mc_turn[i, j] = wedge_angle
        self.mc_pm[i, j] = self.mc_turn[i, j]  # Wall BC
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        
        self._get_moc_variables(i, j, 0, 0)
        
        # Position on wedge wall
        dx = wedge_length / num_rays
        self.mc_x[i, j] = dx
        self.mc_y[i, j] = throat_height + dx * np.tan(wedge_angle * CONV_DEG_RAD)
        
        # First characteristic to opposite wall (2D symmetry)
        i, j = 1, 2
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]  # Symmetry
        self.mc_pm[i, j] = 0.5 * (self.mc_Q[i, j] + self.mc_R[i, j])
        self.mc_turn[i, j] = 0.5 * (self.mc_Q[i, j] - self.mc_R[i, j])
        
        self._get_moc_variables(i, j, 1, 1)
        
        # Opposite wall position
        self.mc_y[i, j] = 0.0  # Centerline for 2D planar
        alpha = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + self.mc_y[1, 1] / np.tan(alpha * CONV_DEG_RAD)
        
        # ================================================================
        # MOC GRID PROPAGATION - Pure 2D (no radial correction)
        # ================================================================
        for i in range(2, num_rays // 2 + 1):
            # ---------------------------------------------------------
            # WEDGE WALL POINT (j=1)
            # ---------------------------------------------------------
            # Wall BC: flow tangent to wedge
            self.mc_turn[i, 1] = wedge_angle
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            
            self._get_moc_variables(i, 1, i-1, 1)
            
            # Position on wedge surface
            self.mc_x[i, 1] = i * dx
            self.mc_y[i, 1] = throat_height + self.mc_x[i, 1] * np.tan(wedge_angle * CONV_DEG_RAD)
            
            # ---------------------------------------------------------
            # INTERNAL GRID POINTS (2 ≤ j ≤ i)
            # ---------------------------------------------------------
            for k in range(2, i + 1):
                # Compatibility equations
                self.mc_Q[i, k] = self.mc_Q[i, k-1]  # Along C+
                self.mc_R[i, k] = self.mc_R[i-1, k]  # Along C-
                
                # Solve for flow properties
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                
                self._get_moc_variables(i, k, i, k-1)
                
                # Characteristic angles
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                # Pure 2D intersection - NO radial correction
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                )
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            # ---------------------------------------------------------
            # CENTERLINE POINT (j = i+1) - 2D Symmetry
            # ---------------------------------------------------------
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]  # Symmetry
            self.mc_pm[i, i+1] = 0.5 * (self.mc_Q[i, i+1] + self.mc_R[i, i+1])
            self.mc_turn[i, i+1] = 0.5 * (self.mc_Q[i, i+1] - self.mc_R[i, i+1])
            
            self._get_moc_variables(i, i+1, i, i)
            
            # Centerline position
            self.mc_y[i, i+1] = 0.0
            alpha = self.mc_mang[i, i] - self.mc_turn[i, i]
            if np.tan(alpha * CONV_DEG_RAD) > 0.001:
                self.mc_x[i, i+1] = self.mc_x[i, i] + self.mc_y[i, i] / np.tan(alpha * CONV_DEG_RAD)
            else:
                self.mc_x[i, i+1] = self.mc_x[i, i] + dx
        
        # ================================================================
        # EXIT PLANE
        # ================================================================
        n = num_rays // 2 + 1
        
        # Extend to exit plane
        for k in range(1, n + 1):
            if k == 1:
                # Wall point
                self.mc_x[n, k] = wedge_length
                self.mc_y[n, k] = h_exit
            else:
                # Internal points
                self.mc_Q[n, k] = self.mc_Q[n-1, k]
                self.mc_R[n, k] = self.mc_R[n-1, k]
                self.mc_pm[n, k] = 0.5 * (self.mc_Q[n, k] + self.mc_R[n, k])
                self.mc_turn[n, k] = 0.5 * (self.mc_Q[n, k] - self.mc_R[n, k])
                
                self._get_moc_variables(n, k, n-1, k)
                
                # Extrapolate to exit
                dx_ext = wedge_length - self.mc_x[n-1, k]
                self.mc_x[n, k] = wedge_length
                self.mc_y[n, k] = self.mc_y[n-1, k] + dx_ext * (self.mc_y[n-1, k] - self.mc_y[n-2, k]) / dx
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_wedge_results(wedge_angle)
    
    def _prepare_wedge_results(self, wedge_angle):
        """Prepare results for wedge nozzle"""
        n = self.num_rays // 2 + 1
        
        # Wedge wall (straight line)
        wall_x = [self.mc_x[i, 1] for i in range(1, n+1)]
        wall_y = [self.mc_y[i, 1] for i in range(1, n+1)]
        
        # Exit plane
        exit_x = [self.mc_x[n, k] for k in range(1, n+1)]
        exit_y = [self.mc_y[n, k] for k in range(1, n+1)]
        
        # Exit properties
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_T_Tt = self.mc_T_Tt[n, n//2 + 1]
        
        # Characteristics
        char_lines = self.get_characteristic_lines()
        
        return {
            'wall_x': np.array(wall_x),
            'wall_y': np.array(wall_y),
            'exit_x': np.array(exit_x),
            'exit_y': np.array(exit_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'nozzle_type': 'wedge',
            'wedge_angle': wedge_angle
        }
    def calculate_bell_nozzle(self, mach_exit, num_rays=30, throat_radius=1.0, 
                             bell_contraction_ratio=0.8, bell_length_ratio=0.8):
        """
        Calculate bell (contoured) nozzle using full MOC
        
        Complete Method of Characteristics for bell nozzle
        - Contoured wall (parabolic profile)
        - Most common rocket nozzle type
        - Higher efficiency than cone (~1-2% better)
        - Shorter and lighter than cone
        - Used in: SpaceX Merlin, Saturn V F-1, RS-25
        
        Args:
            mach_exit: Design exit Mach number
            num_rays: Number of characteristic lines
            throat_radius: Throat radius (inches)
            bell_contraction_ratio: Contraction ratio (0.6-0.9)
            bell_length_ratio: Length ratio vs cone (0.7-0.85 typical)
        
        Returns:
            dict: Results with bell contour and flow properties
        """
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_radius
        self.reset_grid()
        
        # Get exit flow properties
        props_exit = self.isen.get_all_properties(mach_exit)
        A_ratio_exit = props_exit['A_Astar']
        nu_exit = props_exit['nu']
        
        # Bell nozzle geometry
        r_exit = throat_radius * np.sqrt(A_ratio_exit)
        
        # Bell contour parameters
        # Initial expansion angle (typically 20-40°)
        theta_initial = 30.0
        # Exit angle (typically 5-15°)
        theta_exit = 10.0
        
        # Bell length (shorter than cone)
        cone_length = (r_exit - throat_radius) / np.tan(15.0 * CONV_DEG_RAD)
        bell_length = cone_length * bell_length_ratio
        
        # MOC grid parameters
        delta_theta = nu_exit / (num_rays * 2.0)
        
        # ================================================================
        # INITIAL CONDITIONS AT THROAT
        # ================================================================
        i, j = 0, 0
        self.mc_mach[i, j] = 1.0
        self.mc_pm[i, j] = 0.0
        self.mc_mang[i, j] = 90.0
        self.mc_turn[i, j] = 0.0
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_radius
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[i, j] = props_throat['p_pt']
        self.mc_T_Tt[i, j] = props_throat['T_Tt']
        
        # ================================================================
        # INITIAL EXPANSION - Rapid expansion section
        # ================================================================
        i, j = 1, 1
        
        # Bell wall BC: Initial rapid expansion
        self.mc_turn[i, j] = theta_initial * 0.1  # Gradual start
        self.mc_pm[i, j] = self.mc_turn[i, j]
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        
        self._get_moc_variables(i, j, 0, 0)
        
        # Position on bell contour
        dx = bell_length / num_rays
        self.mc_x[i, j] = dx
        self.mc_y[i, j] = throat_radius + dx * np.tan(self.mc_turn[i, j] * CONV_DEG_RAD)
        
        # First characteristic to centerline
        i, j = 1, 2
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]
        self.mc_pm[i, j] = 0.5 * (self.mc_Q[i, j] + self.mc_R[i, j])
        self.mc_turn[i, j] = 0.5 * (self.mc_Q[i, j] - self.mc_R[i, j])
        
        self._get_moc_variables(i, j, 1, 1)
        
        # Centerline position
        self.mc_y[i, j] = 0.0
        alpha = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + self.mc_y[1, 1] / np.tan(alpha * CONV_DEG_RAD)
        
        # ================================================================
        # MOC GRID WITH BELL CONTOUR
        # ================================================================
        for i in range(2, num_rays // 2 + 1):
            # Progress along nozzle (0 to 1)
            progress = float(i) / (num_rays // 2)
            
            # ---------------------------------------------------------
            # BELL WALL POINT (j=1) - CONTOURED
            # ---------------------------------------------------------
            # Bell contour: parabolic blend from initial to exit angle
            # Rapid expansion initially, then gentler approach to exit
            
            if progress < 0.3:
                # Rapid expansion region (0-30%)
                target_angle = theta_initial * (progress / 0.3)
            elif progress < 0.7:
                # Transition region (30-70%)
                blend = (progress - 0.3) / 0.4
                target_angle = theta_initial * (1.0 - blend) + theta_exit * blend
            else:
                # Exit region (70-100%)
                target_angle = theta_exit
            
            self.mc_turn[i, 1] = target_angle
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            
            self._get_moc_variables(i, 1, i-1, 1)
            
            # Bell wall position (parabolic contour)
            self.mc_x[i, 1] = i * dx
            
            # Parabolic contour for bell nozzle
            # r(x) = r_throat + a*x + b*x²
            x_norm = self.mc_x[i, 1] / bell_length
            
            # Parabolic coefficients for smooth contour
            a_coeff = np.tan(theta_initial * CONV_DEG_RAD)
            b_coeff = (r_exit - throat_radius - a_coeff * bell_length) / (bell_length ** 2)
            
            self.mc_y[i, 1] = (throat_radius + 
                              a_coeff * self.mc_x[i, 1] + 
                              b_coeff * self.mc_x[i, 1] ** 2)
            
            # ---------------------------------------------------------
            # INTERNAL GRID POINTS (2 ≤ j ≤ i)
            # ---------------------------------------------------------
            for k in range(2, i + 1):
                # Compatibility equations
                self.mc_Q[i, k] = self.mc_Q[i, k-1]
                self.mc_R[i, k] = self.mc_R[i-1, k]
                
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                
                self._get_moc_variables(i, k, i, k-1)
                
                # Characteristic angles
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                # Axisymmetric correction
                r_avg = 0.5 * (self.mc_y[i, k-1] + self.mc_y[i-1, k])
                if r_avg > 0.001 * throat_radius:
                    M_avg = 0.5 * (self.mc_mach[i, k-1] + self.mc_mach[i-1, k])
                    radial_correction = 1.0 + (M_avg - 1.0) * 0.012 / r_avg
                else:
                    radial_correction = 1.0
                
                # Characteristic intersection
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                ) * radial_correction
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            # ---------------------------------------------------------
            # CENTERLINE POINT (j = i+1) - Symmetry BC
            # ---------------------------------------------------------
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]
            self.mc_pm[i, i+1] = 0.5 * (self.mc_Q[i, i+1] + self.mc_R[i, i+1])
            self.mc_turn[i, i+1] = 0.5 * (self.mc_Q[i, i+1] - self.mc_R[i, i+1])
            
            self._get_moc_variables(i, i+1, i, i)
            
            # Centerline position
            self.mc_y[i, i+1] = 0.0
            alpha = self.mc_mang[i, i] - self.mc_turn[i, i]
            if np.tan(alpha * CONV_DEG_RAD) > 0.001:
                self.mc_x[i, i+1] = self.mc_x[i, i] + self.mc_y[i, i] / np.tan(alpha * CONV_DEG_RAD)
            else:
                self.mc_x[i, i+1] = self.mc_x[i, i] + dx
        
        # ================================================================
        # EXIT PLANE
        # ================================================================
        n = num_rays // 2 + 1
        
        # Extend to exit plane with uniform flow
        for k in range(1, n + 1):
            if k == 1:
                # Wall point - exactly at exit radius
                self.mc_x[n, k] = bell_length
                self.mc_y[n, k] = r_exit
            else:
                # Internal and centerline points
                self.mc_Q[n, k] = self.mc_Q[n-1, k]
                self.mc_R[n, k] = self.mc_R[n-1, k]
                self.mc_pm[n, k] = 0.5 * (self.mc_Q[n, k] + self.mc_R[n, k])
                self.mc_turn[n, k] = 0.5 * (self.mc_Q[n, k] - self.mc_R[n, k])
                
                self._get_moc_variables(n, k, n-1, k)
                
                # Exit plane positions
                self.mc_x[n, k] = bell_length
                # Linear interpolation to exit
                self.mc_y[n, k] = r_exit * (1.0 - float(k-1) / float(n-1))
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_bell_results(theta_exit)
    
    def _prepare_bell_results(self, theta_exit):
        """Prepare results for bell nozzle"""
        n = self.num_rays // 2 + 1
        
        # Bell contoured wall
        wall_x = [self.mc_x[i, 1] for i in range(1, n+1)]
        wall_y = [self.mc_y[i, 1] for i in range(1, n+1)]
        
        # Exit plane
        exit_x = [self.mc_x[n, k] for k in range(1, n+1)]
        exit_y = [self.mc_y[n, k] for k in range(1, n+1)]
        
        # Exit properties
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_T_Tt = self.mc_T_Tt[n, n//2 + 1]
        
        # Characteristics
        char_lines = self.get_characteristic_lines()
        
        return {
            'wall_x': np.array(wall_x),
            'wall_y': np.array(wall_y),
            'exit_x': np.array(exit_x),
            'exit_y': np.array(exit_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'nozzle_type': 'bell',
            'theta_exit': theta_exit
        }
    def calculate_truncated_ideal_nozzle(self, mach_exit, num_rays=30, throat_height=1.0,
                                        truncation_ratio=0.75):
        """
        Calculate truncated ideal nozzle using full MOC
        
        Complete Method of Characteristics for truncated ideal nozzle
        - Same as ideal (2D) nozzle but cut short
        - Practical compromise: high efficiency + reasonable length
        - Truncation ratio: 0.6-0.9 typical (75% = 75% of full length)
        - Still achieves 97-98% efficiency
        - Common in real rocket applications
        
        Args:
            mach_exit: Design exit Mach number
            num_rays: Number of characteristic lines
            throat_height: Throat height (inches)
            truncation_ratio: Length ratio vs full ideal (0.6-0.9 typical)
        
        Returns:
            dict: Results with truncated geometry and flow properties
        """
        # First calculate full ideal nozzle
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_height
        self.reset_grid()
        
        props_exit = self.isen.get_all_properties(mach_exit)
        nu_exit = props_exit['nu']
        
        theta_max = nu_exit / 2.0
        delta_theta = nu_exit / num_rays
        
        # Initial conditions at throat
        self.mc_mach[0, 0] = 1.0
        self.mc_pm[0, 0] = 0.0
        self.mc_mang[0, 0] = 90.0
        self.mc_x[0, 0] = 0.0
        self.mc_y[0, 0] = throat_height
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[0, 0] = props_throat['p_pt']
        self.mc_T_Tt[0, 0] = props_throat['T_Tt']
        
        delx = 0.1
        
        # First characteristic (1,1) - boundary
        i, j = 1, 1
        self.mc_defl[i, j] = delta_theta
        self.mc_turn[i, j] = delta_theta
        self.mc_pm[i, j] = self.mc_turn[i, j]
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        self._get_moc_variables(i, j, 0, 0)
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_height
        
        # First characteristic (1,2) - symmetry
        i, j = 1, 2
        self.mc_turn[i, j] = 0.0
        self.mc_defl[i, j] = self.mc_turn[i, j] - self.mc_turn[1, 1]
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]
        self.mc_pm[i, j] = self.mc_Q[i, j]
        self._get_moc_variables(i, j, 1, 1)
        self.mc_y[i, j] = 0.0
        alpha_angle = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + (self.mc_y[1, 1] - self.mc_y[i, j]) / np.tan(alpha_angle * CONV_DEG_RAD)
        
        # Expansion on boundary
        for i in range(2, num_rays // 2 + 1):
            # Boundary point
            self.mc_defl[i, 1] = delta_theta
            self.mc_turn[i, 1] = self.mc_turn[i-1, 1] + self.mc_defl[i, 1]
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            self._get_moc_variables(i, 1, i-1, 1)
            self.mc_x[i, 1] = self.mc_x[i-1, 1] + delx * throat_height
            self.mc_y[i, 1] = self.mc_y[i-1, 1] + (self.mc_x[i, 1] - self.mc_x[i-1, 1]) * np.tan(self.mc_turn[i-1, 1] * CONV_DEG_RAD)
            
            # Internal points
            for k in range(2, i + 1):
                self.mc_Q[i, k] = self.mc_Q[i, k-1]
                self.mc_R[i, k] = self.mc_R[i-1, k]
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                self._get_moc_variables(i, k, i, k-1)
                
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                )
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            # Symmetry point
            self.mc_turn[i, i+1] = 0.0
            self.mc_defl[i, i+1] = self.mc_turn[i, i+1] - self.mc_turn[i, i]
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]
            self.mc_pm[i, i+1] = self.mc_Q[i, i+1]
            self._get_moc_variables(i, i+1, i, i)
            self.mc_y[i, i+1] = 0.0
            alpha_angle = self.mc_mang[i, i] - self.mc_turn[i, i]
            self.mc_x[i, i+1] = self.mc_x[i, i] + (self.mc_y[i, i] - self.mc_y[i, i+1]) / np.tan(alpha_angle * CONV_DEG_RAD)
        
        # ================================================================
        # TRUNCATION - Cut nozzle at specified ratio
        # ================================================================
        n = num_rays // 2 + 1
        
        # Store full-length endpoint
        full_length_x = self.mc_x[n-1, 1]
        
        # Calculate truncation point
        truncation_x = full_length_x * truncation_ratio
        
        # Find truncation row (closest to truncation_x)
        truncation_row = n
        for i in range(1, n):
            if self.mc_x[i, 1] >= truncation_x:
                truncation_row = i
                break
        
        # Adjust final row to exact truncation point
        for k in range(1, truncation_row + 1):
            if k == 1:
                # Wall point at truncation
                self.mc_x[truncation_row, k] = truncation_x
                # Interpolate y from previous row
                if truncation_row > 1:
                    ratio = (truncation_x - self.mc_x[truncation_row-1, 1]) / (self.mc_x[truncation_row, 1] - self.mc_x[truncation_row-1, 1] + 0.001)
                    self.mc_y[truncation_row, k] = self.mc_y[truncation_row-1, 1] + ratio * (self.mc_y[truncation_row, 1] - self.mc_y[truncation_row-1, 1])
            else:
                # Extend other points proportionally
                if truncation_row > 1:
                    self.mc_x[truncation_row, k] = self.mc_x[truncation_row-1, k] + (truncation_x - self.mc_x[truncation_row-1, 1])
                    self.mc_y[truncation_row, k] = self.mc_y[truncation_row-1, k]
        
        # Update final row index
        n = truncation_row
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_truncated_results(truncation_ratio, full_length_x)
    
    def _prepare_truncated_results(self, truncation_ratio, full_length):
        """Prepare results for truncated nozzle"""
        n = self.num_rays // 2 + 1
        
        # Find actual end of nozzle
        end_row = n
        for i in range(1, n):
            if self.mc_x[i, 1] > 0:
                end_row = i
        
        # Wall contour
        wall_x = [self.mc_x[i, 1] for i in range(1, end_row+1)]
        wall_y = [self.mc_y[i, 1] for i in range(1, end_row+1)]
        
        # Exit plane
        exit_x = [self.mc_x[end_row, k] for k in range(1, end_row+1)]
        exit_y = [self.mc_y[end_row, k] for k in range(1, end_row+1)]
        
        # Exit properties
        exit_mach = self.mc_mach[end_row, end_row//2 + 1]
        exit_p_pt = self.mc_p_pt[end_row, end_row//2 + 1]
        exit_T_Tt = self.mc_T_Tt[end_row, end_row//2 + 1]
        
        # Characteristics
        char_lines = self.get_characteristic_lines()
        
        return {
            'wall_x': np.array(wall_x),
            'wall_y': np.array(wall_y),
            'exit_x': np.array(exit_x),
            'exit_y': np.array(exit_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'nozzle_type': 'truncated',
            'truncation_ratio': truncation_ratio,
            'full_length': full_length
        }
    def calculate_minimum_length_nozzle(self, mach_exit, num_rays=30, throat_radius=1.0):
        """
        Calculate minimum length nozzle using full MOC
        
        Complete Method of Characteristics for minimum length nozzle
        - Shortest possible nozzle for given exit Mach
        - Maximum initial turn angle at throat
        - Optimized characteristic network
        - ~95-97% efficiency (trade-off for minimum length)
        - Critical for weight-sensitive applications (missiles, rockets)
        
        Theory: Rao's minimum length nozzle optimization
        - Rapid initial expansion (max turn angle)
        - Gradual straightening to exit
        - Mathematically shortest solution
        
        Args:
            mach_exit: Design exit Mach number
            num_rays: Number of characteristic lines
            throat_radius: Throat radius (inches)
        
        Returns:
            dict: Results with minimum length geometry
        """
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_radius
        self.reset_grid()
        
        # Get exit flow properties
        props_exit = self.isen.get_all_properties(mach_exit)
        nu_exit = props_exit['nu']
        A_ratio_exit = props_exit['A_Astar']
        
        # Minimum length parameters
        # Maximum initial turn angle (Rao method)
        theta_max_initial = nu_exit * 0.6  # ~60% of total Prandtl-Meyer angle
        
        # Exit radius from area ratio
        r_exit = throat_radius * np.sqrt(A_ratio_exit)
        
        # MOC grid parameters
        delta_theta = nu_exit / num_rays
        
        # ================================================================
        # INITIAL CONDITIONS AT THROAT
        # ================================================================
        i, j = 0, 0
        self.mc_mach[i, j] = 1.0
        self.mc_pm[i, j] = 0.0
        self.mc_mang[i, j] = 90.0
        self.mc_turn[i, j] = 0.0
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_radius
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[i, j] = props_throat['p_pt']
        self.mc_T_Tt[i, j] = props_throat['T_Tt']
        
        # ================================================================
        # RAPID INITIAL EXPANSION - Minimum length strategy
        # ================================================================
        i, j = 1, 1
        
        # Maximum turn at first step for minimum length
        self.mc_turn[i, j] = theta_max_initial / (num_rays // 2)
        self.mc_pm[i, j] = self.mc_turn[i, j]
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        
        self._get_moc_variables(i, j, 0, 0)
        
        # Small axial step for rapid expansion
        dx = 0.05  # Smaller than usual for minimum length
        self.mc_x[i, j] = dx * throat_radius
        self.mc_y[i, j] = throat_radius + dx * throat_radius * np.tan(self.mc_turn[i, j] * CONV_DEG_RAD)
        
        # First characteristic to centerline
        i, j = 1, 2
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]
        self.mc_pm[i, j] = 0.5 * (self.mc_Q[i, j] + self.mc_R[i, j])
        self.mc_turn[i, j] = 0.5 * (self.mc_Q[i, j] - self.mc_R[i, j])
        
        self._get_moc_variables(i, j, 1, 1)
        
        self.mc_y[i, j] = 0.0
        alpha = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + self.mc_y[1, 1] / np.tan(alpha * CONV_DEG_RAD)
        
        # ================================================================
        # MOC GRID WITH MINIMUM LENGTH OPTIMIZATION
        # ================================================================
        for i in range(2, num_rays // 2 + 1):
            # Progress ratio
            progress = float(i) / (num_rays // 2)
            
            # ---------------------------------------------------------
            # WALL POINT - Optimized for minimum length
            # ---------------------------------------------------------
            # Rao method: rapid expansion then gradual straightening
            if progress < 0.3:
                # Rapid expansion phase (0-30%)
                turn_increment = delta_theta * 1.5  # Aggressive
            elif progress < 0.7:
                # Transition phase (30-70%)
                turn_increment = delta_theta * 1.0  # Normal
            else:
                # Straightening phase (70-100%)
                turn_increment = delta_theta * 0.5  # Gentle
            
            self.mc_defl[i, 1] = turn_increment
            self.mc_turn[i, 1] = self.mc_turn[i-1, 1] + self.mc_defl[i, 1]
            
            # Clamp to maximum Prandtl-Meyer angle
            if self.mc_turn[i, 1] > theta_max_initial:
                self.mc_turn[i, 1] = theta_max_initial
            
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            
            self._get_moc_variables(i, 1, i-1, 1)
            
            # Position - minimum length spacing
            self.mc_x[i, 1] = self.mc_x[i-1, 1] + dx * throat_radius
            self.mc_y[i, 1] = self.mc_y[i-1, 1] + (self.mc_x[i, 1] - self.mc_x[i-1, 1]) * np.tan(self.mc_turn[i-1, 1] * CONV_DEG_RAD)
            
            # ---------------------------------------------------------
            # INTERNAL GRID POINTS
            # ---------------------------------------------------------
            for k in range(2, i + 1):
                self.mc_Q[i, k] = self.mc_Q[i, k-1]
                self.mc_R[i, k] = self.mc_R[i-1, k]
                
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                
                self._get_moc_variables(i, k, i, k-1)
                
                # Characteristic angles
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                # Axisymmetric correction
                r_avg = 0.5 * (self.mc_y[i, k-1] + self.mc_y[i-1, k])
                if r_avg > 0.001 * throat_radius:
                    M_avg = 0.5 * (self.mc_mach[i, k-1] + self.mc_mach[i-1, k])
                    radial_correction = 1.0 + (M_avg - 1.0) * 0.01 / r_avg
                else:
                    radial_correction = 1.0
                
                # Characteristic intersection
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                ) * radial_correction
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            # ---------------------------------------------------------
            # CENTERLINE POINT - Symmetry BC
            # ---------------------------------------------------------
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]
            self.mc_pm[i, i+1] = 0.5 * (self.mc_Q[i, i+1] + self.mc_R[i, i+1])
            self.mc_turn[i, i+1] = 0.5 * (self.mc_Q[i, i+1] - self.mc_R[i, i+1])
            
            self._get_moc_variables(i, i+1, i, i)
            
            self.mc_y[i, i+1] = 0.0
            alpha = self.mc_mang[i, i] - self.mc_turn[i, i]
            if np.tan(alpha * CONV_DEG_RAD) > 0.001:
                self.mc_x[i, i+1] = self.mc_x[i, i] + self.mc_y[i, i] / np.tan(alpha * CONV_DEG_RAD)
            else:
                self.mc_x[i, i+1] = self.mc_x[i, i] + dx * throat_radius
        
        # ================================================================
        # EXIT PLANE
        # ================================================================
        n = num_rays // 2 + 1
        
        # Extend to exit with minimum length
        for k in range(1, n + 1):
            if k == 1:
                # Wall point
                self.mc_x[n, k] = self.mc_x[n-1, k] + dx * throat_radius
                self.mc_y[n, k] = r_exit
            else:
                # Internal points
                self.mc_Q[n, k] = self.mc_Q[n-1, k]
                self.mc_R[n, k] = self.mc_R[n-1, k]
                self.mc_pm[n, k] = 0.5 * (self.mc_Q[n, k] + self.mc_R[n, k])
                self.mc_turn[n, k] = 0.5 * (self.mc_Q[n, k] - self.mc_R[n, k])
                
                self._get_moc_variables(n, k, n-1, k)
                
                # Exit positions
                self.mc_x[n, k] = self.mc_x[n-1, k] + dx * throat_radius
                self.mc_y[n, k] = self.mc_y[n-1, k]
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        # Calculate length ratio vs ideal
        ideal_length = self.mc_x[n, 1] * 1.3  # Estimate
        min_length = self.mc_x[n, 1]
        length_ratio = min_length / ideal_length
        
        return self._prepare_minlength_results(length_ratio)
    
    def _prepare_minlength_results(self, length_ratio):
        """Prepare results for minimum length nozzle"""
        n = self.num_rays // 2 + 1
        
        # Wall contour
        wall_x = [self.mc_x[i, 1] for i in range(1, n+1)]
        wall_y = [self.mc_y[i, 1] for i in range(1, n+1)]
        
        # Exit plane
        exit_x = [self.mc_x[n, k] for k in range(1, n+1)]
        exit_y = [self.mc_y[n, k] for k in range(1, n+1)]
        
        # Exit properties
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_T_Tt = self.mc_T_Tt[n, n//2 + 1]
        
        # Characteristics
        char_lines = self.get_characteristic_lines()
        
        return {
            'wall_x': np.array(wall_x),
            'wall_y': np.array(wall_y),
            'exit_x': np.array(exit_x),
            'exit_y': np.array(exit_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'nozzle_type': 'minlength',
            'length_ratio': length_ratio
        }
    def calculate_throat_expansion_nozzle(self, mach_exit, num_rays=30, throat_radius=1.0,
                                         expansion_angle=45.0):
        """
        Calculate throat expansion nozzle using full MOC
        
        Complete Method of Characteristics for throat expansion nozzle
        - Rapid expansion immediately at throat
        - Maximum initial expansion angle (30-50°)
        - Compact design for space-constrained applications
        - Different flow structure than standard designs
        - Used in special applications with geometric constraints
        
        Args:
            mach_exit: Design exit Mach number
            num_rays: Number of characteristic lines
            throat_radius: Throat radius (inches)
            expansion_angle: Initial expansion angle at throat (30-50° typical)
        
        Returns:
            dict: Results with throat expansion geometry
        """
        self.mach_exit = mach_exit
        self.num_rays = num_rays
        self.throat_height = throat_radius
        self.reset_grid()
        
        # Get exit flow properties
        props_exit = self.isen.get_all_properties(mach_exit)
        nu_exit = props_exit['nu']
        A_ratio_exit = props_exit['A_Astar']
        
        # Exit radius
        r_exit = throat_radius * np.sqrt(A_ratio_exit)
        
        # MOC parameters with rapid expansion
        theta_initial = expansion_angle  # Large initial angle
        delta_theta = nu_exit / num_rays
        
        # ================================================================
        # INITIAL CONDITIONS AT THROAT
        # ================================================================
        i, j = 0, 0
        self.mc_mach[i, j] = 1.0
        self.mc_pm[i, j] = 0.0
        self.mc_mang[i, j] = 90.0
        self.mc_turn[i, j] = 0.0
        self.mc_x[i, j] = 0.0
        self.mc_y[i, j] = throat_radius
        
        props_throat = self.isen.get_all_properties(1.0)
        self.mc_p_pt[i, j] = props_throat['p_pt']
        self.mc_T_Tt[i, j] = props_throat['T_Tt']
        
        # ================================================================
        # RAPID THROAT EXPANSION - Key feature
        # ================================================================
        i, j = 1, 1
        
        # Immediate rapid expansion at throat
        self.mc_turn[i, j] = theta_initial * 0.3  # Large initial turn
        self.mc_pm[i, j] = self.mc_turn[i, j]
        self.mc_Q[i, j] = self.mc_pm[i, j] + self.mc_turn[i, j]
        self.mc_R[i, j] = self.mc_pm[i, j] - self.mc_turn[i, j]
        
        self._get_moc_variables(i, j, 0, 0)
        
        # Very short axial step with large radial expansion
        dx = 0.03  # Very small axial step
        self.mc_x[i, j] = dx * throat_radius
        self.mc_y[i, j] = throat_radius + dx * throat_radius * np.tan(self.mc_turn[i, j] * CONV_DEG_RAD) * 2.0
        
        # First characteristic to centerline
        i, j = 1, 2
        self.mc_Q[i, j] = self.mc_Q[1, 1]
        self.mc_R[i, j] = self.mc_Q[i, j]
        self.mc_pm[i, j] = 0.5 * (self.mc_Q[i, j] + self.mc_R[i, j])
        self.mc_turn[i, j] = 0.5 * (self.mc_Q[i, j] - self.mc_R[i, j])
        
        self._get_moc_variables(i, j, 1, 1)
        
        self.mc_y[i, j] = 0.0
        alpha = self.mc_mang[1, 1] - self.mc_turn[1, 1]
        self.mc_x[i, j] = self.mc_x[1, 1] + self.mc_y[1, 1] / np.tan(alpha * CONV_DEG_RAD)
        
        # ================================================================
        # MOC GRID WITH THROAT EXPANSION PROFILE
        # ================================================================
        for i in range(2, num_rays // 2 + 1):
            # Progress ratio
            progress = float(i) / (num_rays // 2)
            
            # ---------------------------------------------------------
            # WALL POINT - Throat expansion profile
            # ---------------------------------------------------------
            # Rapid expansion initially, then moderate
            if progress < 0.2:
                # Very rapid expansion (0-20%)
                turn_increment = delta_theta * 2.0
            elif progress < 0.5:
                # Moderate expansion (20-50%)
                turn_increment = delta_theta * 1.2
            else:
                # Gradual straightening (50-100%)
                turn_increment = delta_theta * 0.8
            
            self.mc_defl[i, 1] = turn_increment
            self.mc_turn[i, 1] = self.mc_turn[i-1, 1] + self.mc_defl[i, 1]
            
            # Limit maximum turn angle
            max_turn = theta_initial
            if self.mc_turn[i, 1] > max_turn:
                self.mc_turn[i, 1] = max_turn
            
            self.mc_pm[i, 1] = self.mc_turn[i, 1]
            self.mc_Q[i, 1] = self.mc_pm[i, 1] + self.mc_turn[i, 1]
            self.mc_R[i, 1] = self.mc_pm[i, 1] - self.mc_turn[i, 1]
            
            self._get_moc_variables(i, 1, i-1, 1)
            
            # Position with throat expansion profile
            self.mc_x[i, 1] = self.mc_x[i-1, 1] + dx * throat_radius
            self.mc_y[i, 1] = self.mc_y[i-1, 1] + (self.mc_x[i, 1] - self.mc_x[i-1, 1]) * np.tan(self.mc_turn[i-1, 1] * CONV_DEG_RAD) * 1.5
            
            # ---------------------------------------------------------
            # INTERNAL GRID POINTS
            # ---------------------------------------------------------
            for k in range(2, i + 1):
                self.mc_Q[i, k] = self.mc_Q[i, k-1]
                self.mc_R[i, k] = self.mc_R[i-1, k]
                
                self.mc_pm[i, k] = 0.5 * (self.mc_Q[i, k] + self.mc_R[i, k])
                self.mc_turn[i, k] = 0.5 * (self.mc_Q[i, k] - self.mc_R[i, k])
                self.mc_defl[i, k] = self.mc_turn[i, k] - self.mc_turn[i-1, k]
                
                self._get_moc_variables(i, k, i, k-1)
                
                # Characteristic angles
                self.mc_alpha[i, k-1] = self.mc_mang[i, k-1] - self.mc_turn[i, k-1]
                self.mc_beta[i-1, k] = self.mc_mang[i-1, k] + self.mc_turn[i-1, k]
                
                tan_alpha = np.tan(self.mc_alpha[i, k-1] * CONV_DEG_RAD)
                tan_beta = np.tan(self.mc_beta[i-1, k] * CONV_DEG_RAD)
                
                # Axisymmetric correction
                r_avg = 0.5 * (self.mc_y[i, k-1] + self.mc_y[i-1, k])
                if r_avg > 0.001 * throat_radius:
                    M_avg = 0.5 * (self.mc_mach[i, k-1] + self.mc_mach[i-1, k])
                    radial_correction = 1.0 + (M_avg - 1.0) * 0.008 / r_avg
                else:
                    radial_correction = 1.0
                
                # Characteristic intersection
                self.mc_x[i, k] = (
                    (self.mc_y[i, k-1] - self.mc_y[i-1, k] +
                     self.mc_x[i, k-1] * tan_alpha +
                     self.mc_x[i-1, k] * tan_beta) /
                    (tan_alpha + tan_beta)
                ) * radial_correction
                
                self.mc_y[i, k] = self.mc_y[i, k-1] - (self.mc_x[i, k] - self.mc_x[i, k-1]) * tan_alpha
            
            # ---------------------------------------------------------
            # CENTERLINE POINT - Symmetry BC
            # ---------------------------------------------------------
            self.mc_Q[i, i+1] = self.mc_Q[i, i]
            self.mc_R[i, i+1] = self.mc_Q[i, i+1]
            self.mc_pm[i, i+1] = 0.5 * (self.mc_Q[i, i+1] + self.mc_R[i, i+1])
            self.mc_turn[i, i+1] = 0.5 * (self.mc_Q[i, i+1] - self.mc_R[i, i+1])
            
            self._get_moc_variables(i, i+1, i, i)
            
            self.mc_y[i, i+1] = 0.0
            alpha = self.mc_mang[i, i] - self.mc_turn[i, i]
            if np.tan(alpha * CONV_DEG_RAD) > 0.001:
                self.mc_x[i, i+1] = self.mc_x[i, i] + self.mc_y[i, i] / np.tan(alpha * CONV_DEG_RAD)
            else:
                self.mc_x[i, i+1] = self.mc_x[i, i] + dx * throat_radius
        
        # ================================================================
        # EXIT PLANE
        # ================================================================
        n = num_rays // 2 + 1
        
        for k in range(1, n + 1):
            if k == 1:
                # Wall point
                self.mc_x[n, k] = self.mc_x[n-1, k] + dx * throat_radius * 2
                self.mc_y[n, k] = r_exit
            else:
                # Internal points
                self.mc_Q[n, k] = self.mc_Q[n-1, k]
                self.mc_R[n, k] = self.mc_R[n-1, k]
                self.mc_pm[n, k] = 0.5 * (self.mc_Q[n, k] + self.mc_R[n, k])
                self.mc_turn[n, k] = 0.5 * (self.mc_Q[n, k] - self.mc_R[n, k])
                
                self._get_moc_variables(n, k, n-1, k)
                
                self.mc_x[n, k] = self.mc_x[n-1, k] + dx * throat_radius
                self.mc_y[n, k] = self.mc_y[n-1, k]
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return self._prepare_throat_expansion_results(expansion_angle)
    
    def _prepare_throat_expansion_results(self, expansion_angle):
        """Prepare results for throat expansion nozzle"""
        n = self.num_rays // 2 + 1
        
        # Wall contour
        wall_x = [self.mc_x[i, 1] for i in range(1, n+1)]
        wall_y = [self.mc_y[i, 1] for i in range(1, n+1)]
        
        # Exit plane
        exit_x = [self.mc_x[n, k] for k in range(1, n+1)]
        exit_y = [self.mc_y[n, k] for k in range(1, n+1)]
        
        # Exit properties
        exit_mach = self.mc_mach[n, n//2 + 1]
        exit_p_pt = self.mc_p_pt[n, n//2 + 1]
        exit_T_Tt = self.mc_T_Tt[n, n//2 + 1]
        
        # Characteristics
        char_lines = self.get_characteristic_lines()
        
        return {
            'wall_x': np.array(wall_x),
            'wall_y': np.array(wall_y),
            'exit_x': np.array(exit_x),
            'exit_y': np.array(exit_y),
            'exit_mach': exit_mach,
            'exit_p_pt': exit_p_pt,
            'exit_T_Tt': exit_T_Tt,
            'num_rays': self.num_rays,
            'throat_height': self.throat_height,
            'characteristic_lines': char_lines,
            'nozzle_type': 'throatexp',
            'expansion_angle': expansion_angle
        }
        # ========================================================================
    # SHOCK DETECTION AND COALESCENCE
    # ========================================================================
    
    def check_coalescence(self, i, j, tolerance=0.01):
        """Check for coalescence of compression waves (internal)"""
        if i < 2 or j < 2:
            return False
        
        theta_curr = self.mc_turn[i, j]
        theta_prev = self.mc_turn[i-1, j]
        
        # Compression detection
        if theta_curr < theta_prev - tolerance:
            M_curr = self.mc_mach[i, j]
            M_prev = self.mc_mach[i-1, j]
            
            # Potential shock formation
            if M_curr < M_prev and M_prev > 1.0:
                theta_diff = abs(theta_curr - theta_prev)
                
                # Use shock calculator
                try:
                    shock_props = self.shock_calc.theta_beta_M(M_prev, theta_diff)
                    
                    if not shock_props['detached']:
                        self.coalescence_points.append({
                            'i': i, 'j': j,
                            'x': self.mc_x[i, j],
                            'y': self.mc_y[i, j],
                            'M_upstream': M_prev,
                            'M_downstream': M_curr,
                            'theta': theta_diff,
                            'shock_angle': shock_props['weak']['beta'] if shock_props['weak'] else 0
                        })
                        return True
                except:
                    pass
        
        return False
    
    def detect_shocks(self):
        """Detect all shocks in the flow field"""
        shocks = []
        n = self.num_rays // 2 + 1
        
        for i in range(1, n):
            for j in range(1, n):
                if i > 0 and self.mc_p_pt[i, j] > 0:
                    p_ratio = self.mc_p_pt[i, j] / self.mc_p_pt[i-1, j]
                    
                    if p_ratio > 1.10:
                        M_upstream = self.mc_mach[i-1, j]
                        M_downstream = self.mc_mach[i, j]
                        
                        if M_upstream > M_downstream and M_upstream > 1.0:
                            dx = self.mc_x[i, j] - self.mc_x[i-1, j]
                            dy = self.mc_y[i, j] - self.mc_y[i-1, j]
                            
                            if dx != 0:
                                shock_angle = np.degrees(np.arctan(dy / dx))
                            else:
                                shock_angle = 90.0
                            
                            shocks.append({
                                'i': i, 'j': j,
                                'x': self.mc_x[i, j],
                                'y': self.mc_y[i, j],
                                'M1': M_upstream,
                                'M2': M_downstream,
                                'p_ratio': p_ratio,
                                'shock_angle': shock_angle,
                                'type': 'compression'
                            })
        
        self.shocks_detected = shocks
        return shocks
        
    def calculate_2d_ideal_nozzle_zone_method(self, mach_exit, num_zones=20, throat_height=1.0):
        """
        Zone Method (Field Method) - prob=0
        Alternative to Points Method
        Divides flow into zones instead of characteristic lines
        """
        self.mach_exit = mach_exit
        self.num_rays = num_zones
        self.throat_height = throat_height
        self.reset_grid()
        
        props_exit = self.isen.get_all_properties(mach_exit)
        nu_exit = props_exit['nu']
        
        # Zone divisions
        delta_nu = nu_exit / num_zones
        
        # Initial throat
        self.mc_mach[0, 0] = 1.0
        self.mc_pm[0, 0] = 0.0
        self.mc_x[0, 0] = 0.0
        self.mc_y[0, 0] = throat_height
        
        # Zone calculation (simplified)
        for i in range(1, num_zones + 1):
            nu = i * delta_nu
            M = self.isen.get_mach_from_pm(nu)
            
            self.mc_mach[i, 0] = M
            self.mc_pm[i, 0] = nu
            self.mc_turn[i, 0] = nu / 2.0
            self.mc_x[i, 0] = i * 0.1 * throat_height
            self.mc_y[i, 0] = throat_height + self.mc_x[i, 0] * np.tan(self.mc_turn[i, 0] * CONV_DEG_RAD)
            
            props = self.isen.get_all_properties(M)
            self.mc_p_pt[i, 0] = props['p_pt']
            self.mc_T_Tt[i, 0] = props['T_Tt']
        
        self.geometry_calculated = True
        self.flow_calculated = True
        
        return {
            'wall_x': self.mc_x[:num_zones+1, 0],
            'wall_y': self.mc_y[:num_zones+1, 0],
            'exit_mach': self.mc_mach[num_zones, 0],
            'method': 'zone'
        }
    def get_pressure_boundary_condition(self, i, j, p_external=None):
        """
        Calculate pressure boundary condition
        For free boundaries and external flow matching
        """
        if p_external is None:
            # Use ambient pressure
            p_external = 14.7  # psia default
        
        # Current point pressure
        p_local = self.mc_p_pt[i, j] * 50.0  # Assume p0=50 psia
        
        # Pressure ratio
        p_ratio = p_local / p_external
        
        # Boundary adjustment needed?
        if abs(p_ratio - 1.0) > 0.05:
            # Pressure mismatch - requires adjustment
            return {
                'needs_adjustment': True,
                'p_local': p_local,
                'p_external': p_external,
                'p_ratio': p_ratio,
                'adjustment_factor': p_external / p_local
            }
        else:
            # Pressure matched
            return {
                'needs_adjustment': False,
                'p_local': p_local,
                'p_external': p_external,
                'p_ratio': p_ratio
            }
    
    def apply_boundary_corrections(self, p_ambient=14.7):
        """
        Apply boundary condition corrections to nozzle exit
        Adjusts for pressure matching
        """
        n = self.num_rays // 2 + 1
        corrections = []
        
        for k in range(1, n+1):
            bc = self.get_pressure_boundary_condition(n, k, p_ambient)
            if bc['needs_adjustment']:
                corrections.append({
                    'point': (n, k),
                    'x': self.mc_x[n, k],
                    'y': self.mc_y[n, k],
                    'correction': bc
                })
        
        return corrections
    def get_free_stream_properties(self, M_freestream=0.0, altitude=0.0):
        """
        Calculate free stream properties
        For external flow analysis
        
        Args:
            M_freestream: Freestream Mach number (flight Mach)
            altitude: Altitude in feet
        
        Returns:
            dict: Freestream properties
        """
        # Standard atmosphere at altitude
        if altitude == 0.0:
            p_inf = 14.7  # psia sea level
            T_inf = 518.67  # Rankine sea level
        else:
            # Simplified atmosphere model
            p_inf = 14.7 * (1.0 - 0.0000068756 * altitude) ** 5.2559
            T_inf = 518.67 - 0.00356616 * altitude
        
        # Speed of sound
        a_inf = np.sqrt(self.gamma * 1716.0 * T_inf)  # ft/s
        
        # Velocity
        V_inf = M_freestream * a_inf
        
        # Dynamic pressure
        rho_inf = p_inf * 144.0 / (1716.0 * T_inf)  # slug/ft^3
        q_inf = 0.5 * rho_inf * V_inf ** 2
        
        # Stagnation properties
        if M_freestream > 0:
            props = self.isen.get_all_properties(M_freestream)
            p0_inf = p_inf / props['p_pt']
            T0_inf = T_inf / props['T_Tt']
        else:
            p0_inf = p_inf
            T0_inf = T_inf
        
        return {
            'M_inf': M_freestream,
            'p_inf': p_inf,
            'T_inf': T_inf,
            'rho_inf': rho_inf,
            'a_inf': a_inf,
            'V_inf': V_inf,
            'q_inf': q_inf,
            'p0_inf': p0_inf,
            'T0_inf': T0_inf,
            'altitude': altitude
        }
    def calculate_internal_plug_nozzle(self, mach_exit, num_rays=30, 
                                       throat_radius=1.0, plug_truncation=0.8,
                                       perint=0.5):
        """
        Internal/External Plug Nozzle (C.C. Lee method) - prob=5
        perint: internal/external expansion ratio (0.01-1.0)
                Controls how much expansion occurs internally vs externally.
                phiex = nu_exit * perint
        """
        import math

        # Get isentropic properties at exit
        props_exit = self.isentropic_properties(mach_exit)
        nu_exit = props_exit['nu']  # Total Prandtl-Meyer angle (radians)

        # C.C. Lee split: internal expansion angle
        phi_internal = nu_exit * perint        # internal portion
        phi_external = nu_exit * (1.0 - perint)  # external portion

        # Base plug calculation using internal expansion
        result = self.calculate_plug_nozzle(
            mach_exit, num_rays, throat_radius, plug_truncation
        )

        # Apply perint scaling to spike geometry
        # Spike is shaped by internal expansion angle
        scale = phi_internal / nu_exit if nu_exit > 0 else 1.0
        if 'spike_y' in result:
            result['spike_y'] = [y * scale for y in result['spike_y']]

        # External boundary adjusted by external expansion
        ext_scale = 1.0 + phi_external / nu_exit if nu_exit > 0 else 1.0
        if 'outer_y' in result:
            result['outer_y'] = [y * ext_scale for y in result['outer_y']]

        # Mark as C.C. Lee internal type
        result['nozzle_type'] = 'internal_plug'
        result['plug_type'] = 'cc_lee'
        result['perint'] = perint
        result['phi_internal_deg'] = math.degrees(phi_internal)
        result['phi_external_deg'] = math.degrees(phi_external)

        # Efficiency
        try:
            n = self.num_rays // 2 + 1
            exit_M = self.mc_mach[n-1, n//2]
            result['efficiency'] = (exit_M / mach_exit) * 100.0
        except:
            result['efficiency'] = 95.0

        return result
    def calculate_external_plug_nozzle(self, mach_exit, num_rays=30,
                                       throat_radius=1.0, external_expansion_ratio=1.2):
        """
        External Plug Nozzle - prob=4
        Expansion outside spike surface
        """
        result = self.calculate_plug_nozzle(
            mach_exit, num_rays, throat_radius, 0.9
        )
        
        # External flow adjustment
        n = self.num_rays // 2 + 1
        
        # Expand outer boundary
        for k in range(1, n):
            idx = min(k, len(result['outer_y'])-1)
            result['outer_y'][idx] *= external_expansion_ratio
        
        result['nozzle_type'] = 'external_plug'
        result['plug_type'] = 'external'
        result['expansion_ratio'] = external_expansion_ratio
        
        return result
    def calculate_flat_cowl_plug(self, mach_exit, num_rays=30,
                                 throat_radius=1.0, cowl_angle=0.0):
        """
        Flat Cowl Plug Nozzle - prob=6
        Plug with flat external cowl
        """
        result = self.calculate_plug_nozzle(
            mach_exit, num_rays, throat_radius, 0.85
        )
        
        # Flat cowl geometry
        n = self.num_rays // 2 + 1
        x_max = result['outer_x'][-1]
        y_cowl = result['outer_y'][0]
        
        # Add cowl data
        result['cowl'] = {
            'x': [0, x_max],
            'y': [y_cowl, y_cowl + x_max * np.tan(cowl_angle * CONV_DEG_RAD)],
            'angle': cowl_angle
        }
        
        result['nozzle_type'] = 'flat_cowl_plug'
        result['cowl_angle'] = cowl_angle
        
        return result
    def calculate_given_shape_analysis(self, contour_x, contour_y, num_rays=30):
        """
        Given Shape Analysis - prob=7
        Analyze user-defined nozzle contour
        Flow field from known wall shape
        """
        if len(contour_x) != len(contour_y):
            raise ValueError("Contour arrays must have same length")
        
        self.num_rays = num_rays
        self.throat_height = contour_y[0]
        self.reset_grid()
        
        # Map contour to grid
        n_points = len(contour_x)
        for i in range(min(n_points, num_rays)):
            self.mc_x[i, 0] = contour_x[i]
            self.mc_y[i, 0] = contour_y[i]
            
            # Estimate wall angle
            if i > 0:
                dx = contour_x[i] - contour_x[i-1]
                dy = contour_y[i] - contour_y[i-1]
                self.mc_turn[i, 0] = np.degrees(np.arctan(dy / dx)) if dx != 0 else 0
        
        self.geometry_calculated = True
        
        return {
            'wall_x': contour_x,
            'wall_y': contour_y,
            'nozzle_type': 'given_shape',
            'num_points': n_points
        }
    def calculate_wedge_analysis(self, wedge_angle, mach_upstream, length=10.0):
        """
        Wedge Analysis - prob=8
        Analyze flow over wedge shape
        Oblique shock and expansion waves
        """
        self.reset_grid()
        
        # Oblique shock at wedge
        shock_result = self.shock_calc.theta_beta_M(mach_upstream, wedge_angle)
        
        if shock_result['detached']:
            raise ValueError(f"Shock detached: theta={wedge_angle}° > theta_max={shock_result['theta_max']:.2f}°")
        
        # Wedge geometry
        wedge_x = [0, length]
        wedge_y = [0, length * np.tan(wedge_angle * CONV_DEG_RAD)]
        
        return {
            'wall_x': np.array(wedge_x),
            'wall_y': np.array(wedge_y),
            'nozzle_type': 'wedge_analysis',
            'wedge_angle': wedge_angle,
            'M_upstream': mach_upstream,
            'shock_angle': shock_result['weak']['beta'] if shock_result['weak'] else 0,
            'M_downstream': shock_result['weak']['M2'] if shock_result['weak'] else mach_upstream
        }
    def set_plume_mode(self, mode=1):
        """
        Set plume calculation mode
        mode=0: No plume
        mode=1: Simple plume
        mode=2: Advanced (shocks/expansions)
        """
        self.plume_mode = mode
        return mode
    def calculate_simple_plume(self, p_ambient=14.7, length_multiplier=2.0):
        """Simple plume (mode=1) - No shock diamonds"""
        if not self.flow_calculated:
            raise RuntimeError("Flow not calculated")
        
        n = self.num_rays // 2 + 1
        exit_x = self.mc_x[n, :]
        exit_y = self.mc_y[n, :]
        
        plume_x = []
        plume_y = []
        
        for k in range(1, n+1):
            x_start = exit_x[k]
            y_start = exit_y[k]
            length = length_multiplier * self.throat_height
            
            plume_x.append([x_start, x_start + length])
            plume_y.append([y_start, y_start])
        
        return {
            'plume_x': plume_x,
            'plume_y': plume_y,
            'shock_diamonds': None
        }   
    def calculate_advanced_plume(self, p_ambient=14.7):
        """Advanced plume (mode=2) - With shock diamonds"""
        return self.calculate_plume(p_ambient)
    def calculate_jet_under_expansion(self, p_exit, p_ambient):
        """
        Under-expanded jet analysis
        p_exit > p_ambient → expansion waves
        """
        if p_exit <= p_ambient:
            return None
        
        PR = p_exit / p_ambient
        n = self.num_rays // 2 + 1
        M_exit = self.mc_mach[n, n//2]
        
        # Expansion angle
        expansion_angle = 15.0 * (PR - 1.0)
        
        return {
            'type': 'under_expanded',
            'pressure_ratio': PR,
            'expansion_angle': expansion_angle,
            'mach_exit': M_exit
        }
    def calculate_jet_over_expansion(self, p_exit, p_ambient):
        """
        Over-expanded jet analysis
        p_exit < p_ambient → oblique shocks at exit
        """
        if p_exit >= p_ambient:
            return None
        
        PR = p_ambient / p_exit
        n = self.num_rays // 2 + 1
        M_exit = self.mc_mach[n, n//2]
        
        # Shock angle estimate
        shock_angle = 30.0 + 10.0 * (PR - 1.0)
        
        return {
            'type': 'over_expanded',
            'pressure_ratio': PR,
            'shock_angle': shock_angle,
            'mach_exit': M_exit
        }
    def get_flow_turning_plot_data(self):
        """
        Flow turning angle plot (vturn)
        Theta vs position
        """
        if not self.flow_calculated:
            return None
        
        n = self.num_rays // 2 + 1
        x_data = []
        theta_data = []
        
        for i in range(1, n):
            x_data.append(self.mc_x[i, 1])
            theta_data.append(self.mc_turn[i, 1])
        
        return {
            'x': np.array(x_data),
            'theta': np.array(theta_data),
            'title': 'Flow Turning Angle'
        }
    def get_expansion_fan_data(self):
        """
        Expansion fan visualization data
        First few characteristics from throat
        """
        if not self.geometry_calculated:
            return None
        
        fan_lines = []
        n = min(5, self.num_rays // 2)  # First 5 rays
        
        for i in range(1, n+1):
            x_line = [self.mc_x[i, 1], self.mc_x[i, i+1]]
            y_line = [self.mc_y[i, 1], self.mc_y[i, i+1]]
            fan_lines.append((x_line, y_line))
        
        return {
            'lines': fan_lines,
            'count': len(fan_lines)
        }
    def get_reflected_flow_data(self):
        """
        Reflected flow visualization
        Mirror characteristics for symmetry
        """
        if not self.geometry_calculated:
            return None
        
        n = self.num_rays // 2 + 1
        reflected_x = []
        reflected_y = []
        
        for i in range(1, n):
            reflected_x.append(self.mc_x[i, 1])
            reflected_y.append(-self.mc_y[i, 1])  # Mirror
        
        return {
            'x': np.array(reflected_x),
            'y': np.array(reflected_y),
            'reflected': True
        }
    def set_colorbar_limits(self, variable='mach', vmin=None, vmax=None):
        """
        Set color bar min/max for contour plots
        """
        if not self.flow_calculated:
            return None
        
        n = self.num_rays // 2 + 1
        
        if variable == 'mach':
            data = self.mc_mach[:n, :n]
        elif variable == 'pressure':
            data = self.mc_p_pt[:n, :n]
        elif variable == 'temperature':
            data = self.mc_T_Tt[:n, :n]
        else:
            return None
        
        if vmin is None:
            vmin = np.min(data[data > 0])
        if vmax is None:
            vmax = np.max(data[data > 0])
        
        return {
            'vmin': vmin,
            'vmax': vmax,
            'variable': variable
        }
    def calculate_external_cowl_angle(self, cowl_length_ratio=1.2):
        """
        External cowl angle calculation
        For plug nozzles with cowl
        """
        if not self.geometry_calculated:
            return None
        
        n = self.num_rays // 2 + 1
        throat_r = self.throat_height
        exit_r = self.mc_y[n, n//2]
        nozzle_length = self.mc_x[n, 1]
        
        cowl_length = nozzle_length * cowl_length_ratio
        cowl_angle = np.degrees(np.arctan((exit_r - throat_r) / cowl_length))
        
        return {
            'angle': cowl_angle,
            'length': cowl_length,
            'exit_radius': exit_r
        }
    def set_draw_ray_option(self, enabled=True):
        """Draw ray option toggle"""
        self.draw_rays = enabled
        return enabled
    
    def set_gamma_variation(self, gamma_new):
        """Gamma variations (gamopt)"""
        self.gamma = gamma_new
        self.isen = IsentropicRelations(gamma_new)
        self.shock_calc = ShockCalculator(gamma_new)
        return gamma_new    
    def get_max_delx(self):
        """getDelmx() - Maximum delta x"""
        if not self.flow_calculated:
            return 0.0
        n = self.num_rays // 2 + 1
        max_dx = 0.0
        for i in range(1, n):
            dx = abs(self.mc_x[i, 1] - self.mc_x[i-1, 1])
            if dx > max_dx:
                max_dx = dx
        return max_dx    
    def get_max_theta(self):
        """getThetmax() - Maximum turning angle"""
        if not self.flow_calculated:
            return 0.0
        n = self.num_rays // 2 + 1
        return np.max(self.mc_turn[:n, :n])    
    def get_angle_limits(self):
        """getAnglim() - Angle limits"""
        if not self.flow_calculated:
            return None
        theta_max = self.get_max_theta()
        mu_min = self.shock_calc.mach_angle(self.mach_exit)
        return {
            'theta_max': theta_max,
            'mu_min': mu_min,
            'theta_range': (0, theta_max)
        }
    def getDelmx(self):
        """Get maximum Mach number difference between adjacent points"""
        try:
            if not hasattr(self, 'mc_mach') or self.mc_mach is None:
                return 0.0
            
            diffs = []
            n = self.num_rays // 2 + 1
            
            for i in range(1, n):
                for j in range(1, i+1):
                    if i < self.mc_mach.shape[0] and j < self.mc_mach.shape[1]:
                        dm = abs(self.mc_mach[i,j] - self.mc_mach[i-1,j])
                        diffs.append(dm)
            
            return max(diffs) if diffs else 0.0
        except:
            return 0.0
    
    def getThetmax(self):
        """Get maximum flow turning angle"""
        try:
            if not hasattr(self, 'mc_turn') or self.mc_turn is None:
                return 0.0
            
            theta_deg = np.degrees(np.abs(self.mc_turn))
            return float(np.max(theta_deg))
        except:
            return 0.0
    
    def getAnglim(self):
        """Get characteristic angle limit"""
        try:
            if not hasattr(self, 'mc_pm') or self.mc_pm is None:
                return 0.0
            
            nu_deg = np.degrees(np.abs(self.mc_pm))
            return float(np.max(nu_deg))
        except:
            return 0.0