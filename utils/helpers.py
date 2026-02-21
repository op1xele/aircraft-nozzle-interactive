"""
Helper Functions + Unit Conversion System
Supports Imperial (default) and SI unit systems.
"""


def format_number(value, decimal_places=3):
    try:
        return f"{float(value):.{decimal_places}f}"
    except (ValueError, TypeError):
        return "N/A"


def safe_divide(numerator, denominator, default=0.0):
    try:
        if abs(denominator) < 1e-10:
            return default
        return numerator / denominator
    except (TypeError, ZeroDivisionError):
        return default


# ============================================================================
# UNIT CONVERSION
# ============================================================================

# Active unit system — change this to switch globally
_UNIT_SYSTEM = 'imperial'  # 'imperial' or 'si'


def set_unit_system(system):
    global _UNIT_SYSTEM
    if system not in ('imperial', 'si'):
        raise ValueError("system must be 'imperial' or 'si'")
    _UNIT_SYSTEM = system


def get_unit_system():
    return _UNIT_SYSTEM


# ------------------------------------------------------------------
# Conversion factors: Imperial → SI
# ------------------------------------------------------------------
_CONV = {
    # length
    'in_to_m':    0.0254,
    'm_to_in':    1.0 / 0.0254,
    # pressure
    'psia_to_Pa': 6894.76,
    'Pa_to_psia': 1.0 / 6894.76,
    # temperature (handled as functions below)
    # velocity
    'fps_to_mps': 0.3048,
    'mps_to_fps': 1.0 / 0.3048,
    # mass flow
    'lbms_to_kgs': 0.453592,
    'kgs_to_lbms': 1.0 / 0.453592,
    # force
    'lbf_to_N':   4.44822,
    'N_to_lbf':   1.0 / 4.44822,
    # mass
    'lbm_to_kg':  0.453592,
    'kg_to_lbm':  1.0 / 0.453592,
}


def temp_to_display(T_imperial):
    """Convert temperature from Rankine to display unit."""
    if _UNIT_SYSTEM == 'si':
        return T_imperial * 5.0 / 9.0   # Rankine → Kelvin
    return T_imperial                    # Rankine


def temp_from_display(T_display):
    """Convert temperature from display unit to Rankine."""
    if _UNIT_SYSTEM == 'si':
        return T_display * 9.0 / 5.0    # Kelvin → Rankine
    return T_display


def pressure_to_display(p_imperial):
    """Convert pressure from psia to display unit."""
    if _UNIT_SYSTEM == 'si':
        return p_imperial * _CONV['psia_to_Pa'] / 1000.0  # psia → kPa
    return p_imperial


def pressure_from_display(p_display):
    """Convert pressure from display unit to psia."""
    if _UNIT_SYSTEM == 'si':
        return p_display * 1000.0 * _CONV['Pa_to_psia']   # kPa → psia
    return p_display


def length_to_display(l_imperial):
    """Convert length from inches to display unit."""
    if _UNIT_SYSTEM == 'si':
        return l_imperial * _CONV['in_to_m'] * 100.0      # in → cm
    return l_imperial


def length_from_display(l_display):
    """Convert length from display unit to inches."""
    if _UNIT_SYSTEM == 'si':
        return l_display / 100.0 / _CONV['in_to_m']       # cm → in
    return l_display


def velocity_to_display(v_imperial):
    """Convert velocity from ft/s to display unit."""
    if _UNIT_SYSTEM == 'si':
        return v_imperial * _CONV['fps_to_mps']
    return v_imperial


def force_to_display(f_imperial):
    """Convert force from lbf to display unit."""
    if _UNIT_SYSTEM == 'si':
        return f_imperial * _CONV['lbf_to_N']
    return f_imperial


def massflow_to_display(mdot_imperial):
    """Convert mass flow from lbm/s to display unit."""
    if _UNIT_SYSTEM == 'si':
        return mdot_imperial * _CONV['lbms_to_kgs']
    return mdot_imperial


# ------------------------------------------------------------------
# Unit label helpers
# ------------------------------------------------------------------
UNITS = {
    'imperial': {
        'pressure':    'psia',
        'temperature': 'R',
        'length':      'in',
        'velocity':    'ft/s',
        'force':       'lbf',
        'massflow':    'lbm/s',
        'cstar':       'ft/s',
        'isp':         's',
    },
    'si': {
        'pressure':    'kPa',
        'temperature': 'K',
        'length':      'cm',
        'velocity':    'm/s',
        'force':       'N',
        'massflow':    'kg/s',
        'cstar':       'm/s',
        'isp':         's',
    }
}


def unit_label(quantity):
    """Return unit label for current system."""
    return UNITS.get(_UNIT_SYSTEM, UNITS['imperial']).get(quantity, '')