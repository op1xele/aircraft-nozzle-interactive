"""
Aircraft Nozzle Interactive - Configuration and Constants
Version: 2.0 (Advanced)
Author: Pickle Brothers Engineering Team
"""

import numpy as np

# ============================================================================
# MATHEMATICAL CONSTANTS
# ============================================================================
PI = np.pi
CONV_DEG_RAD = PI / 180.0
CONV_RAD_DEG = 180.0 / PI

# Gas constants
R_GAS = 1716.0  # ft-lbf/(lbm-R) for air
R_AIR = 1716.0  # Same as R_GAS
GC = 32.2       # ft/s^2

# ============================================================================
# DEFAULT FLOW CONDITIONS
# ============================================================================
DEFAULT_GAMMA = 1.4
DEFAULT_MACH_EXIT = 2.5
DEFAULT_TOTAL_PRESSURE = 50.0   # psia (p0)
DEFAULT_TOTAL_TEMP = 3000.0     # Rankine (T0)
DEFAULT_PTO = 50.0              # psia (legacy name)
DEFAULT_TTO = 3000.0            # Rankine (legacy name)
DEFAULT_ALTITUDE = 35000.0      # ft

# ============================================================================
# NOZZLE DESIGN PARAMETERS
# ============================================================================
DEFAULT_NUM_RAYS = 30    # Number of characteristic lines
MIN_NUM_RAYS = 10
MAX_NUM_RAYS = 100

DEFAULT_DELX = 0.1       # Grid spacing
MIN_DELX = 0.01
MAX_DELX = 1.0

DEFAULT_NOZZLE_LENGTH = 1.0
DEFAULT_THROAT_HEIGHT = 1.0

# ============================================================================
# REAL ROCKET ENGINE PRESETS
# ============================================================================
ROCKET_PRESETS = {
    "SpaceX Merlin 1D": {
        "description": "Falcon 9 first stage • RP-1/LOX • Sea level optimized",
        "p0": 1470,      # psia (101 bar)
        "T0": 6332,      # Rankine (3518 K)
        "mach_exit": 2.8,
        "throat_radius": 4.8,  # inches
        "nozzle_type": 6,  # Bell nozzle
    },
    "RS-25 (SSME)": {
        "description": "Space Shuttle Main Engine • LH2/LOX • High performance",
        "p0": 3000,      # psia (206.8 bar)
        "T0": 6170,      # Rankine (3428 K)
        "mach_exit": 4.5,
        "throat_radius": 5.2,
        "nozzle_type": 6,  # Bell nozzle
    },
    "RD-180": {
        "description": "Atlas V • RP-1/LOX • Russian design • Very high pressure",
        "p0": 3722,      # psia (256.5 bar)
        "T0": 6386,      # Rankine (3548 K)
        "mach_exit": 3.2,
        "throat_radius": 4.5,
        "nozzle_type": 6,
    },
    "Saturn V F-1": {
        "description": "Saturn V first stage • RP-1/LOX • Largest rocket engine",
        "p0": 1015,      # psia (70 bar)
        "T0": 6170,      # Rankine
        "mach_exit": 2.7,
        "throat_radius": 11.8,  # inches (very large!)
        "nozzle_type": 6,
    },
    "SpaceX Raptor": {
        "description": "Starship • CH4/LOX • Full-flow staged combustion",
        "p0": 4350,      # psia (300 bar) - highest pressure
        "T0": 6332,      # Rankine
        "mach_exit": 3.5,
        "throat_radius": 5.0,
        "nozzle_type": 6,
    },
    "RL-10": {
        "description": "Centaur upper stage • LH2/LOX • Vacuum optimized",
        "p0": 475,       # psia (32.7 bar)
        "T0": 4824,      # Rankine (2680 K) - hydrogen fuel
        "mach_exit": 4.2,
        "throat_radius": 3.8,
        "nozzle_type": 6,
    },
    "J-2": {
        "description": "Saturn V upper stage • LH2/LOX • Vacuum engine",
        "p0": 785,       # psia (54 bar)
        "T0": 5400,      # Rankine
        "mach_exit": 4.0,
        "throat_radius": 6.2,
        "nozzle_type": 6,
    },
    "BE-4": {
        "description": "Blue Origin Vulcan • CH4/LOX • American methalox",
        "p0": 2000,      # psia (138 bar)
        "T0": 6200,      # Rankine
        "mach_exit": 3.0,
        "throat_radius": 5.5,
        "nozzle_type": 6,
    },
    "Custom": {
        "description": "Custom values (user input)",
        "p0": 50.0,
        "T0": 3000.0,
        "mach_exit": 2.5,
        "throat_radius": 1.0,
        "nozzle_type": 1,
    }
}

# ============================================================================
# PROBLEM TYPES
# ============================================================================
PROBLEM_TYPES = {
    1: "2D Ideal Nozzle",
    2: "Axisymmetric Ideal Nozzle",
    3: "Plug Nozzle (Aerospike)",
    4: "Cone Nozzle",
    5: "Wedge Nozzle (2D)",
    6: "Bell Nozzle",
    7: "Truncated Ideal Nozzle",
    8: "Minimum Length Nozzle",
    9: "Throat Expansion Nozzle"
}

DEFAULT_PROBLEM = 1  # 2D Ideal Nozzle

# ============================================================================
# PLUME OPTIONS
# ============================================================================
PLUME_NONE = 0
PLUME_STATIC = 1
PLUME_SUPERSONIC = 2

DEFAULT_PLUME = 0

# ============================================================================
# GUI SETTINGS
# ============================================================================
WINDOW_WIDTH = 1400
WINDOW_HEIGHT = 900

# Canvas settings
CANVAS_WIDTH = 800
CANVAS_HEIGHT = 600
CANVAS_BG_COLOR = "#FFFFFF"

# Colors
COLOR_NOZZLE = "#2C3E50"
COLOR_CHARACTERISTICS = "#3498DB"
COLOR_BOUNDARY = "#E74C3C"
COLOR_SYMMETRY = "#2ECC71"
COLOR_GRID = "#BDC3C7"

# Plot colors (for contours)
COLORMAP_PRESSURE = "jet"
COLORMAP_TEMPERATURE = "hot"
COLORMAP_MACH = "viridis"

# ============================================================================
# CALCULATION LIMITS
# ============================================================================
MAX_ITERATIONS = 1000
TOLERANCE = 1e-6

# Grid size
MAX_GRID_SIZE = 200

# Mach number limits
MIN_MACH = 1.0
MAX_MACH = 10.0

# Pressure limits
MIN_PRESSURE = 0.1      # psia
MAX_PRESSURE = 10000.0  # psia

# Temperature limits
MIN_TEMPERATURE = 100.0   # Rankine
MAX_TEMPERATURE = 10000.0 # Rankine

# ============================================================================
# FILE SETTINGS
# ============================================================================
PROJECT_FILE_EXTENSION = ".pbproj"
EXPORT_IMAGE_DPI = 300

# ============================================================================
# VISUALIZATION SETTINGS
# ============================================================================
# Display modes
DISPLAY_MESH = "mesh"
DISPLAY_CONTOUR = "contour"
DISPLAY_GEOMETRY = "geometry"
DISPLAY_STREAMLINES = "streamlines"
DISPLAY_MACH = "mach"
DISPLAY_PRESSURE = "pressure"
DISPLAY_TEMP = "temperature"

# Contour levels
NUM_CONTOUR_LEVELS = 20

# Animation settings
ANIMATION_FPS = 30
ANIMATION_DURATION = 5  # seconds

# ============================================================================
# ADVANCED FEATURES
# ============================================================================
ENABLE_3D_VIEW = True
ENABLE_OPTIMIZATION = True
ENABLE_COMPARISON_MODE = True
ENABLE_EXPORT = True

# Optimization settings
OPTIMIZATION_ALGORITHMS = ["Gradient Descent", "Genetic Algorithm", "Particle Swarm"]
DEFAULT_OPTIMIZATION = "Gradient Descent"

# ============================================================================
# UNITS SYSTEM
# ============================================================================
UNITS = {
    'length': 'inches',
    'pressure': 'psia',
    'temperature': 'Rankine',
    'velocity': 'ft/s',
    'mass_flow': 'lbm/s'
}

# ============================================================================
# THEME SETTINGS
# ============================================================================
THEMES = {
    'light': {
        'bg': '#FFFFFF',
        'fg': '#000000',
        'panel_bg': '#F5F5F5',
        'button': '#3498DB',
        'button_hover': '#2980B9'
    },
    'dark': {
        'bg': '#2C3E50',
        'fg': '#ECF0F1',
        'panel_bg': '#34495E',
        'button': '#E74C3C',
        'button_hover': '#C0392B'
    }
}

DEFAULT_THEME = 'light'

# ============================================================================
# VERSION INFO
# ============================================================================
VERSION = "1.0.0"
BUILD_DATE = "2026-02-07"
AUTHOR = "Pickle Brothers Engineering Team"
LICENSE = "Pickle Brothers - All Rights Reserved"
WINDOW_TITLE = "Aircraft Nozzle Interactive - Pickle Brothers"