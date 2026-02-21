"""
Input Validators
"""

from config import *


class InputValidator:
    """Validate user inputs"""
    
    @staticmethod
    def validate_mach(mach):
        try:
            m = float(mach)
            if m < MIN_MACH or m > MAX_MACH:
                return False, f"Mach must be between {MIN_MACH} and {MAX_MACH}"
            return True, m
        except ValueError:
            return False, "Invalid number"
    
    @staticmethod
    def validate_pressure(pressure):
        try:
            p = float(pressure)
            if p < MIN_PRESSURE or p > MAX_PRESSURE:
                return False, f"Pressure must be between {MIN_PRESSURE} and {MAX_PRESSURE} psia"
            return True, p
        except ValueError:
            return False, "Invalid number"
    
    @staticmethod
    def validate_temperature(temperature):
        try:
            t = float(temperature)
            if t < MIN_TEMPERATURE or t > MAX_TEMPERATURE:
                return False, f"Temperature must be between {MIN_TEMPERATURE} and {MAX_TEMPERATURE} R"
            return True, t
        except ValueError:
            return False, "Invalid number"
    
    @staticmethod
    def validate_num_rays(num_rays):
        try:
            n = int(num_rays)
            if n < MIN_NUM_RAYS or n > MAX_NUM_RAYS:
                return False, f"Number of rays must be between {MIN_NUM_RAYS} and {MAX_NUM_RAYS}"
            if n % 2 != 0:
                return False, "Number of rays must be even"
            return True, n
        except ValueError:
            return False, "Invalid integer"