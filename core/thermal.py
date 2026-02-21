"""
Thermal Analysis Module
Nozzle wall temperature distribution and cooling requirements.

Methods:
- Adiabatic wall temperature
- Heat flux distribution (Bartz correlation)
- Regenerative cooling channel sizing
- Material temperature limits
"""

import numpy as np


# ============================================================================
# MATERIAL DATABASE
# ============================================================================
MATERIALS = {
    'Inconel 718': {
        'density':        8200,    # kg/m³
        'k':              11.4,    # W/(m·K) thermal conductivity
        'E':              200e9,   # Pa Young's modulus
        'yield_strength': 1100e6,  # Pa
        'max_temp_K':     1200,    # K service limit
        'alpha':          13e-6,   # 1/K thermal expansion
        'description':    'Nickel superalloy — most common rocket nozzle material',
    },
    'Copper (C18150)': {
        'density':        8900,
        'k':              350,
        'E':              115e9,
        'yield_strength': 310e6,
        'max_temp_K':     700,
        'alpha':          17e-6,
        'description':    'High conductivity — regeneratively cooled chambers',
    },
    'Titanium (Ti-6Al-4V)': {
        'density':        4430,
        'k':              6.7,
        'E':              114e9,
        'yield_strength': 880e6,
        'max_temp_K':     600,
        'alpha':          8.6e-6,
        'description':    'Lightweight — upper stages, low heat flux',
    },
    'Carbon-Carbon': {
        'density':        1800,
        'k':              40,
        'E':              70e9,
        'yield_strength': 200e6,
        'max_temp_K':     2500,
        'alpha':          1e-6,
        'description':    'Ultra high temp — ablative/reentry nozzles',
    },
    'Stainless 316L': {
        'density':        8000,
        'k':              16,
        'E':              193e9,
        'yield_strength': 290e6,
        'max_temp_K':     870,
        'alpha':          16e-6,
        'description':    'General purpose — low cost, moderate performance',
    },
    'Rhenium': {
        'density':        21020,
        'k':              48,
        'E':              460e9,
        'yield_strength': 1070e6,
        'max_temp_K':     2200,
        'alpha':          6.7e-6,
        'description':    'Extreme temp — thruster nozzles, very expensive',
    },
}


# ============================================================================
# THERMAL ANALYSIS
# ============================================================================

class ThermalAnalysis:
    """
    Nozzle wall thermal analysis using Bartz correlation for heat flux
    and 1D radial conduction for wall temperature.
    """

    def __init__(self, gamma=1.4, T0_R=3000.0, p0_psia=50.0,
                 throat_radius_in=1.0, material='Inconel 718'):
        self.gamma         = gamma
        self.T0_R          = T0_R
        self.T0_K          = T0_R * 5.0/9.0
        self.p0_psia       = p0_psia
        self.p0_Pa         = p0_psia * 6894.76
        self.r_t_in        = throat_radius_in
        self.r_t_m         = throat_radius_in * 0.0254
        self.material      = MATERIALS.get(material, MATERIALS['Inconel 718'])
        self.material_name = material

    # ------------------------------------------------------------------
    # Bartz correlation — heat transfer coefficient
    # ------------------------------------------------------------------
    def bartz_htc(self, mach_dist, wall_x_in, wall_y_in):
        """
        Heat transfer coefficient h [W/(m²·K)] via Bartz (1957).

        h = (0.026 / D_t^0.2) * (mu^0.2 * Cp / Pr^0.6) * (p0/c*)^0.8 * (D_t/R_c)^0.1 * (A_t/A)^0.9 * sigma
        """
        g   = self.gamma
        gm1 = g - 1.0
        gp1 = g + 1.0

        # Gas properties (air-like, approximate)
        Cp_gas  = 1000.0 + 200.0*(self.T0_K/3000.0)  # J/(kg·K), simple fit
        Pr      = 0.85 - 0.1*(g - 1.4)                # Prandtl number
        mu_ref  = 1.8e-5 * (self.T0_K/300.0)**0.7     # Pa·s

        # C* (characteristic velocity)
        cstar = (np.sqrt(g*287*self.T0_K) /
                 (g * np.sqrt((2.0/gp1)**((gp1)/gm1))))

        D_t = 2.0 * self.r_t_m  # throat diameter [m]

        wall_y_m = np.array(wall_y_in) * 0.0254
        throat_y_m = wall_y_m[0]

        htc = np.zeros(len(mach_dist))

        for i, M in enumerate(mach_dist):
            M = max(M, 0.1)
            factor  = 1.0 + 0.5*gm1*M**2
            T_ratio = 1.0/factor
            p_ratio = T_ratio**(g/gm1)

            # Area ratio
            A_ratio = (wall_y_m[i] / throat_y_m)**2 if throat_y_m > 0 else 1.0
            A_ratio = max(A_ratio, 1.0)

            # Adiabatic wall temperature
            r_factor = Pr**0.33  # recovery factor
            T_aw = self.T0_K * (T_ratio + r_factor * 0.5*gm1*M**2 * T_ratio)

            # Film coefficient correction (sigma)
            T_wall_est = min(T_aw * 0.6, self.material['max_temp_K'])
            sigma_num = 0.5*(T_wall_est/self.T0_K)*(1.0 + 0.5*gm1*M**2) + 0.5
            sigma = sigma_num**(-0.68) * (1.0 + 0.5*gm1*M**2)**(-0.12)
            sigma = np.clip(sigma, 0.1, 2.0)

            # Bartz coefficient
            h = (0.026 / D_t**0.2 *
                 (mu_ref**0.2 * Cp_gas / Pr**0.6) *
                 (self.p0_Pa / cstar)**0.8 *
                 (D_t / (1.5*self.r_t_m))**0.1 *
                 (1.0/A_ratio)**0.9 *
                 sigma)

            htc[i] = max(h, 0.0)

        return htc

    # ------------------------------------------------------------------
    # Wall temperature distribution
    # ------------------------------------------------------------------
    def wall_temperature(self, mach_dist, wall_x_in, wall_y_in,
                         wall_thickness_mm=3.0, coolant_temp_K=300.0):
        """
        Compute hot-gas-side and cold-side wall temperatures.

        Returns dict with temperature arrays and heat flux.
        """
        g   = self.gamma
        gm1 = g - 1.0
        Pr  = 0.85

        htc = self.bartz_htc(mach_dist, wall_x_in, wall_y_in)
        t_wall = wall_thickness_mm / 1000.0  # m
        k_wall = self.material['k']

        T_hot  = np.zeros(len(mach_dist))
        T_cold = np.zeros(len(mach_dist))
        q_flux = np.zeros(len(mach_dist))
        T_aw   = np.zeros(len(mach_dist))

        for i, M in enumerate(mach_dist):
            M = max(M, 0.1)
            factor  = 1.0 + 0.5*gm1*M**2
            T_ratio = 1.0/factor
            r_factor = Pr**0.33
            T_aw[i] = self.T0_K * (T_ratio + r_factor*0.5*gm1*M**2*T_ratio)

            # Overall heat transfer: q = (T_aw - T_coolant) / (1/h + t/k)
            R_total = 1.0/htc[i] + t_wall/k_wall if htc[i] > 0 else 1e6
            q_flux[i] = (T_aw[i] - coolant_temp_K) / R_total

            T_hot[i]  = T_aw[i] - q_flux[i]/htc[i] if htc[i] > 0 else T_aw[i]
            T_cold[i] = T_hot[i] - q_flux[i]*t_wall/k_wall

        return {
            'T_adiabatic_wall_K': T_aw,
            'T_hot_wall_K':       T_hot,
            'T_cold_wall_K':      T_cold,
            'heat_flux_MW_m2':    q_flux / 1e6,
            'htc':                htc,
            'max_T_hot_K':        float(np.max(T_hot)),
            'max_q_MW_m2':        float(np.max(q_flux/1e6)),
            'material_limit_K':   self.material['max_temp_K'],
            'material_ok':        bool(np.max(T_hot) < self.material['max_temp_K']),
        }

    # ------------------------------------------------------------------
    # Cooling channel sizing
    # ------------------------------------------------------------------
    def cooling_channel_sizing(self, thermal_results, coolant='water',
                               coolant_flow_kg_s=0.5):
        """
        Estimate regenerative cooling channel dimensions.
        Simple 1D energy balance.
        """
        q_total_MW = float(np.mean(thermal_results['heat_flux_MW_m2']))

        # Coolant properties
        coolants = {
            'water':    {'Cp': 4186, 'rho': 1000, 'mu': 1e-3,   'k': 0.6},
            'RP-1':     {'Cp': 2000, 'rho': 800,  'mu': 2e-3,   'k': 0.14},
            'LH2':      {'Cp': 14310,'rho': 71,   'mu': 1.3e-5, 'k': 0.1},
            'LOX':      {'Cp': 1700, 'rho': 1141, 'mu': 2e-4,   'k': 0.15},
        }
        cp = coolants.get(coolant, coolants['water'])

        # Temperature rise of coolant
        Q_total = q_total_MW * 1e6 * (2*np.pi*self.r_t_m * 0.1)  # approx area
        dT_coolant = Q_total / (coolant_flow_kg_s * cp['Cp'])

        # Channel sizing (assume square channels)
        n_channels = max(20, int(2*np.pi*self.r_t_m*1000 / 5))  # ~5mm pitch
        q_per_channel = Q_total / n_channels
        Re_needed = 1e4  # turbulent
        D_ch = (Re_needed * cp['mu'] / cp['rho'] /
                (coolant_flow_kg_s/n_channels*cp['rho'])) if coolant_flow_kg_s > 0 else 0.003
        D_ch = np.clip(D_ch, 0.001, 0.01)

        return {
            'n_channels':         n_channels,
            'channel_diameter_mm': D_ch*1000,
            'coolant_dT_K':       dT_coolant,
            'total_heat_load_kW': Q_total/1000,
            'coolant':            coolant,
        }

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    def get_material_list(self):
        return list(MATERIALS.keys())

    def check_material_suitability(self, max_T_K):
        results = {}
        for name, mat in MATERIALS.items():
            margin = mat['max_temp_K'] - max_T_K
            results[name] = {
                'max_temp_K': mat['max_temp_K'],
                'margin_K':   margin,
                'suitable':   margin > 0,
                'description': mat['description'],
            }
        return results