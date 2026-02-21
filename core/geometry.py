"""
Geometry Calculator Module
"""

import numpy as np


class GeometryCalculator:
    """Calculate and manage nozzle geometry"""
    
    def __init__(self):
        self.contour_points = []
        self.characteristic_lines = []
    
    def calculate_throat_area(self, throat_height, width=1.0):
        return throat_height * width
    
    def calculate_exit_area(self, exit_height, width=1.0):
        return exit_height * width
    
    def calculate_area_ratio(self, exit_height, throat_height):
        return exit_height / throat_height
    
    def scale_geometry(self, x_points, y_points, scale_factor):
        x_scaled = np.array(x_points) * scale_factor
        y_scaled = np.array(y_points) * scale_factor
        return x_scaled, y_scaled
    
    def reflect_about_centerline(self, x_points, y_points):
        x_upper = np.array(x_points)
        y_upper = np.array(y_points)
        x_lower = x_upper[::-1]
        y_lower = -y_upper[::-1]
        x_full = np.concatenate([x_upper, x_lower])
        y_full = np.concatenate([y_upper, y_lower])
        return x_full, y_full
    
    def calculate_nozzle_length(self, x_points):
        return np.max(x_points) - np.min(x_points)