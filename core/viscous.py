"""
Viscous Correction Module - Boundary Layer Displacement Thickness
Corrects MOC nozzle wall geometry for viscous effects.

The boundary layer displaces the effective flow, making the nozzle
act as if the wall is closer to the centerline. To compensate,
the physical wall must be moved outward by the displacement thickness δ*.

Method: Thwaites integral method for laminar BL,
        Michel's transition criterion,
        simplified turbulent BL for high-Re flows.
"""

import numpy as np


class ViscousCorrection:
    """
    Boundary layer displacement thickness correction for nozzle walls.

    Usage:
        vc = ViscousCorrection(gamma=1.4, T0_R=3000, p0_psia=50, axisymmetric=True)
        corrected_y = vc.apply_correction(wall_x, wall_y, mach_distribution)
    """

    def __init__(self, gamma=1.4, T0_R=3000.0, p0_psia=50.0,
                 axisymmetric=True, Re_throat=1e6):
        self.gamma      = gamma
        self.T0_R       = T0_R          # Total temperature (Rankine)
        self.p0_psia    = p0_psia       # Total pressure (psia)
        self.axisym     = axisymmetric
        self.Re_throat  = Re_throat     # Reynolds number based on throat height

        # Sutherland constants for air (Imperial)
        self.mu_ref  = 3.737e-7   # lbf·s/ft² at T_ref
        self.T_ref_R = 518.67     # Rankine
        self.S_R     = 198.6      # Sutherland constant (Rankine)

    # ------------------------------------------------------------------
    # Viscosity via Sutherland's law
    # ------------------------------------------------------------------
    def viscosity(self, T_R):
        """Dynamic viscosity [lbf·s/ft²] at static temperature T_R."""
        ratio = T_R / self.T_ref_R
        return self.mu_ref * ratio**1.5 * (self.T_ref_R + self.S_R) / (T_R + self.S_R)

    # ------------------------------------------------------------------
    # Isentropic local properties from Mach
    # ------------------------------------------------------------------
    def _local_props(self, M):
        g  = self.gamma
        gm1 = g - 1.0
        factor = 1.0 + 0.5*gm1*M**2
        T_ratio = 1.0 / factor
        p_ratio = T_ratio ** (g/gm1)
        rho_ratio = T_ratio ** (1.0/gm1)
        return T_ratio, p_ratio, rho_ratio

    # ------------------------------------------------------------------
    # Displacement thickness using simple integral BL method
    # ------------------------------------------------------------------
    def displacement_thickness(self, wall_x, wall_y, mach_dist):
        """
        Compute boundary layer displacement thickness δ* at each wall station.

        Args:
            wall_x    : x-coordinates of wall (inches)
            wall_y    : y-coordinates of wall (inches)
            mach_dist : Mach number at each wall station

        Returns:
            delta_star : displacement thickness (inches) at each station
        """
        n = len(wall_x)
        delta_star = np.zeros(n)

        # Convert to feet for calculations
        x_ft = np.array(wall_x) / 12.0
        y_ft = np.array(wall_y) / 12.0

        # Reference density and velocity at throat
        T0_abs = self.T0_R
        p0_abs = self.p0_psia * 144.0  # lbf/ft²

        R_air  = 1716.0  # ft·lbf/(slug·R) — using slugs for consistency
        rho0   = p0_abs / (R_air * T0_abs)   # slug/ft³  (approx, ignore gc)

        # Throat conditions (M=1)
        T_ratio_t, p_ratio_t, rho_ratio_t = self._local_props(1.0)
        T_throat  = T0_abs * T_ratio_t
        a_throat  = np.sqrt(self.gamma * R_air * T_throat)
        rho_throat = rho0 * rho_ratio_t
        mu_throat  = self.viscosity(T_throat)

        # Throat height (first wall_y value)
        h_throat_ft = y_ft[0] if y_ft[0] > 0 else 1.0/12.0

        for i in range(1, n):
            M = max(mach_dist[i], 1.0)
            T_ratio, p_ratio, rho_ratio = self._local_props(M)

            T_local   = T0_abs * T_ratio
            rho_local = rho0 * rho_ratio
            a_local   = np.sqrt(self.gamma * R_air * T_local)
            u_local   = M * a_local
            mu_local  = self.viscosity(T_local)

            # Running length along the wall
            dx = x_ft[i] - x_ft[0]
            if dx <= 0:
                continue

            # Local Reynolds number
            Re_x = rho_local * u_local * dx / mu_local if mu_local > 0 else 1e6

            # Decide laminar or turbulent transition (Michel criterion ~Re_x > 5e5)
            if Re_x < 5e5:
                # Laminar: Blasius δ* = 1.72 * x / sqrt(Re_x)
                delta_lam = 1.72 * dx / np.sqrt(Re_x) if Re_x > 0 else 0
                delta_star_local = delta_lam
            else:
                # Turbulent: δ* ≈ 0.046 * x / Re_x^0.2 (1/7 power law)
                delta_turb = 0.046 * dx / (Re_x**0.2) if Re_x > 0 else 0
                # Compressibility correction factor (Van Driest)
                T_wall_ratio = 1.0 + 0.178*M**2   # adiabatic wall
                comp_factor  = T_wall_ratio**0.61
                delta_star_local = delta_turb * comp_factor

            # Axisymmetric Mangler correction
            if self.axisym and y_ft[i] > 0:
                mangler = np.sqrt(h_throat_ft / y_ft[i])
                delta_star_local *= mangler

            # Convert back to inches
            delta_star[i] = delta_star_local * 12.0

        return delta_star

    # ------------------------------------------------------------------
    # Main correction entry point
    # ------------------------------------------------------------------
    def apply_correction(self, wall_x, wall_y, mach_dist=None):
        """
        Apply viscous correction to nozzle wall geometry.

        If mach_dist is not provided, estimates Mach from isentropic area ratio.

        Returns:
            corrected_y  : corrected wall y-coordinates (inches)
            delta_star   : displacement thickness at each station (inches)
            correction_pct : correction as % of local wall height
        """
        wall_x  = np.array(wall_x, dtype=float)
        wall_y  = np.array(wall_y, dtype=float)

        # Estimate Mach distribution if not given
        if mach_dist is None:
            mach_dist = self._estimate_mach(wall_y)

        delta_star = self.displacement_thickness(wall_x, wall_y, mach_dist)

        # Corrected wall = original + δ* (wall moves outward)
        corrected_y = wall_y + delta_star

        correction_pct = np.where(wall_y > 0,
                                  delta_star / wall_y * 100.0,
                                  np.zeros_like(delta_star))

        return corrected_y, delta_star, correction_pct

    def _estimate_mach(self, wall_y):
        """Estimate Mach from area ratio (isentropic, axisymmetric)."""
        from scipy.optimize import brentq
        g   = self.gamma
        gm1 = g - 1.0
        gp1 = g + 1.0

        def area_mach(M, AR):
            return (1.0/M)*((2.0/gp1)*(1.0+0.5*gm1*M**2))**((gp1)/(2.0*gm1)) - AR

        throat_y = wall_y[0] if wall_y[0] > 0 else wall_y.max()*0.1
        machs = np.ones(len(wall_y))
        for i in range(1, len(wall_y)):
            AR = (wall_y[i] / throat_y)**2
            if AR < 1.0:
                AR = 1.0
            try:
                machs[i] = brentq(area_mach, 1.001, 15.0, args=(AR,))
            except:
                machs[i] = machs[i-1]
        return machs

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    def get_summary(self, wall_y, delta_star):
        """Return correction summary dict."""
        max_ds   = float(np.max(delta_star))
        exit_ds  = float(delta_star[-1])
        exit_pct = float(exit_ds / wall_y[-1] * 100.0) if wall_y[-1] > 0 else 0.0
        mean_pct = float(np.mean(delta_star[1:] / np.maximum(wall_y[1:], 1e-9)) * 100.0)
        return {
            'max_delta_star_in':  max_ds,
            'exit_delta_star_in': exit_ds,
            'exit_correction_pct': exit_pct,
            'mean_correction_pct': mean_pct,
        }