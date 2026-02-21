"""
Output Panels Module - Enhanced
"""

import tkinter as tk
from tkinter import ttk
from utils import format_number


class OutputPanels:
    """Output panels for displaying results"""

    def __init__(self, parent):
        self.parent = parent
        self.frame = ttk.LabelFrame(parent, text="Output Results", padding=10)
        self.frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self._shock_frame = None
        self._coal_frame = None
        self._plume_frame = None
        self._limits_frame = None

        self.create_geometry_output()
        self.create_flow_output()
        self.create_performance_output()
        self.create_advanced_output()

    # ------------------------------------------------------------------
    def create_geometry_output(self):
        frame = ttk.LabelFrame(self.frame, text="Geometry", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=3)

        from utils.helpers import unit_label
        rows = [
            ("Exit Height:",        "exit_height_label",  unit_label('length')),
            ("Throat Height:",      "throat_height_label",unit_label('length')),
            ("Area Ratio (Ae/At):", "area_ratio_label",   ""),
            ("Nozzle Length:",      "length_label",       unit_label('length')),
            ("Length/Throat:",      "l_d_label",          ""),
        ]
        for i, (text, attr, unit) in enumerate(rows):
            ttk.Label(frame, text=text).grid(row=i, column=0, sticky=tk.W, pady=1)
            lbl = ttk.Label(frame, text="---", foreground='blue', width=10)
            lbl.grid(row=i, column=1, sticky=tk.W, padx=5, pady=1)
            setattr(self, attr, lbl)
            if unit:
                ttk.Label(frame, text=unit).grid(row=i, column=2, sticky=tk.W)

    # ------------------------------------------------------------------
    def create_flow_output(self):
        frame = ttk.LabelFrame(self.frame, text="Flow Properties at Exit", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=3)

        from utils.helpers import unit_label
        rows = [
            ("Exit Mach Number:", "exit_mach_label",    ""),
            ("p_e / p_t:",        "exit_p_ratio_label", ""),
            ("T_e / T_t:",        "exit_T_ratio_label", ""),
            ("Exit Pressure:",    "exit_pressure_label",unit_label('pressure')),
            ("Exit Temperature:", "exit_temp_label",    unit_label('temperature')),
            ("Exit Velocity:",    "exit_velocity_label",unit_label('velocity')),
            ("Mach Angle (μ):",   "mach_angle_label",   "°"),
        ]
        for i, (text, attr, unit) in enumerate(rows):
            ttk.Label(frame, text=text).grid(row=i, column=0, sticky=tk.W, pady=1)
            lbl = ttk.Label(frame, text="---", foreground='blue', width=10)
            lbl.grid(row=i, column=1, sticky=tk.W, padx=5, pady=1)
            setattr(self, attr, lbl)
            if unit:
                ttk.Label(frame, text=unit).grid(row=i, column=2, sticky=tk.W)

    # ------------------------------------------------------------------
    def create_performance_output(self):
        frame = ttk.LabelFrame(self.frame, text="Performance Metrics", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=3)

        from utils.helpers import unit_label
        rows = [
            ("Thrust:",           "thrust_label",      unit_label('force')),
            ("Mass Flow:",        "mass_flow_label",   unit_label('massflow')),
            ("NPR (p_t/p_a):",    "npr_label",         ""),
            ("Isp:",              "isp_label",         "sec"),
            ("C* (C-star):",      "cstar_label",       unit_label('cstar')),
            ("Thrust Coeff Cf:",  "cf_label",          ""),
            ("Efficiency:",       "efficiency_label",  "%"),
            ("Expansion:",        "expansion_label",   ""),
        ]
        for i, (text, attr, unit) in enumerate(rows):
            ttk.Label(frame, text=text).grid(row=i, column=0, sticky=tk.W, pady=1)
            lbl = ttk.Label(frame, text="---", foreground='blue', width=12)
            lbl.grid(row=i, column=1, sticky=tk.W, padx=5, pady=1)
            setattr(self, attr, lbl)
            if unit:
                ttk.Label(frame, text=unit).grid(row=i, column=2, sticky=tk.W)

    # ------------------------------------------------------------------
    def create_advanced_output(self):
        """Shock, coalescence, plume indicators — created once, updated dynamically."""
        self._shock_frame = ttk.LabelFrame(self.frame, text="Shock Statistics", padding=5)
        self._shock_frame.pack(fill=tk.X, padx=5, pady=3)
        self._shock_labels = {}
        for key, text in [("num","Shocks detected:"), ("max_pr","Max pressure ratio:"),
                          ("avg_ang","Avg shock angle:")]:
            ttk.Label(self._shock_frame, text=text).pack(anchor=tk.W)
            lbl = ttk.Label(self._shock_frame, text="---", foreground='blue')
            lbl.pack(anchor=tk.W, padx=10)
            self._shock_labels[key] = lbl

        self._coal_frame = ttk.LabelFrame(self.frame, text="Coalescence", padding=5)
        self._coal_frame.pack(fill=tk.X, padx=5, pady=3)
        self._coal_label = ttk.Label(self._coal_frame, text="---")
        self._coal_label.pack(anchor=tk.W)
        self._coal_warn = ttk.Label(self._coal_frame, text="", foreground='orange',
                                    font=('Arial', 8, 'bold'))
        self._coal_warn.pack(anchor=tk.W)

        self._plume_frame = ttk.LabelFrame(self.frame, text="Plume", padding=5)
        self._plume_frame.pack(fill=tk.X, padx=5, pady=3)
        self._plume_type_lbl = ttk.Label(self._plume_frame, text="---",
                                         font=('Arial', 9, 'bold'))
        self._plume_type_lbl.pack(anchor=tk.W)
        self._plume_pr_lbl = ttk.Label(self._plume_frame, text="")
        self._plume_pr_lbl.pack(anchor=tk.W)

    # ------------------------------------------------------------------
    def update_results(self, results, params=None):
        if results is None:
            self.clear_results()
            return

        import math

        # --- Geometry ---
        throat_h = results.get('throat_height', params.get('throat_height', 1.0) if params else 1.0)
        exit_h = results['wall_y'][-1] if 'wall_y' in results else None

        if exit_h is not None:
            self.exit_height_label.config(text=format_number(exit_h, 3))
        self.throat_height_label.config(text=format_number(throat_h, 3))

        if exit_h is not None:
            ar = exit_h / throat_h
            self.area_ratio_label.config(text=format_number(ar, 3))

        if 'wall_x' in results:
            length = results['wall_x'][-1] - results['wall_x'][0]
            self.length_label.config(text=format_number(length, 3))
            self.l_d_label.config(text=format_number(length / throat_h, 2))

        # --- Flow ---
        if 'exit_mach' in results:
            M = results['exit_mach']
            self.exit_mach_label.config(text=format_number(M, 4))
            if M >= 1.0:
                mu = math.degrees(math.asin(1.0/M))
                self.mach_angle_label.config(text=format_number(mu, 2))

        if 'exit_p_pt' in results:
            self.exit_p_ratio_label.config(text=format_number(results['exit_p_pt'], 5))

        if 'exit_T_Tt' in results:
            self.exit_T_ratio_label.config(text=format_number(results['exit_T_Tt'], 5))

        if 'performance' in results:
            perf = results['performance']
            self.exit_pressure_label.config(text=format_number(perf.get('exit_pressure', 0), 3))
            self.exit_temp_label.config(text=format_number(perf.get('exit_temperature', 0), 1))
            self.exit_velocity_label.config(text=format_number(perf.get('exit_velocity', 0), 1))

        # --- Performance ---
        if 'performance' in results:
            perf = results['performance']
            self.thrust_label.config(text=format_number(perf.get('thrust', 0), 2))
            self.mass_flow_label.config(text=format_number(perf.get('mass_flow', 0), 4))
            self.npr_label.config(text=format_number(perf.get('NPR', 0), 2))
            self.isp_label.config(text=format_number(perf.get('Isp', 0), 1))
            self.cstar_label.config(text=format_number(perf.get('C_star', 0), 1))
            self.cf_label.config(text=format_number(perf.get('C_F', 0), 4))
            self.efficiency_label.config(text=format_number(perf.get('efficiency', 0), 1))
            exp = perf.get('expansion_condition', '---')
            color = 'green' if 'Perfectly' in exp else ('red' if 'Over' in exp else 'blue')
            self.expansion_label.config(text=exp, foreground=color)

        # --- Shock stats ---
        if 'shock_stats' in results:
            s = results['shock_stats']
            self._shock_labels['num'].config(text=str(s.get('num_shocks', 0)))
            self._shock_labels['max_pr'].config(text=format_number(s.get('max_pressure_ratio', 0), 3))
            self._shock_labels['avg_ang'].config(text=format_number(s.get('avg_shock_angle', 0), 1) + "°")
        else:
            for lbl in self._shock_labels.values():
                lbl.config(text="N/A")

        # --- Coalescence ---
        if 'coalescence_data' in results:
            cd = results['coalescence_data']
            total = cd.get('total_count', 0)
            self._coal_label.config(
                text=f"Total: {total}  (C+: {cd.get('plus_count',0)}  C-: {cd.get('minus_count',0)})")
            if total == 0:
                self._coal_warn.config(text="No coalescence (good!)", foreground='green')
            elif total > 10:
                self._coal_warn.config(text="⚠ High — consider more rays", foreground='orange')
            else:
                self._coal_warn.config(text="Minor coalescence", foreground='blue')
        else:
            self._coal_label.config(text="---")
            self._coal_warn.config(text="")

        # --- Plume ---
        if 'plume_data' in results and results['plume_data']:
            p0 = params.get('p0', 50.0) if params else results.get('p0', 50.0)
            p_exit_pt = results.get('exit_p_pt', 0)
            p_exit = p_exit_pt * p0
            p_amb = params.get('p_ambient', 14.7) if params else 14.7
            PR = p_exit / p_amb if p_amb > 0 else 0

            if PR > 1.05:
                ptype, color = "Under-Expanded (p_e > p_a)", 'blue'
            elif PR < 0.95:
                ptype, color = "Over-Expanded (p_e < p_a)", 'red'
            else:
                ptype, color = "Perfectly Expanded", 'green'

            self._plume_type_lbl.config(text=ptype, foreground=color)
            self._plume_pr_lbl.config(text=f"Pressure Ratio: {PR:.3f}")
        else:
            self._plume_type_lbl.config(text="No plume", foreground='gray')
            self._plume_pr_lbl.config(text="")

    # ------------------------------------------------------------------
    def clear_results(self):
        attrs = [
            'exit_height_label', 'throat_height_label', 'area_ratio_label',
            'length_label', 'l_d_label',
            'exit_mach_label', 'exit_p_ratio_label', 'exit_T_ratio_label',
            'exit_pressure_label', 'exit_temp_label', 'exit_velocity_label',
            'mach_angle_label',
            'thrust_label', 'mass_flow_label', 'npr_label', 'isp_label',
            'cstar_label', 'cf_label', 'efficiency_label', 'expansion_label',
        ]
        for attr in attrs:
            if hasattr(self, attr):
                getattr(self, attr).config(text="---", foreground='blue')
        if hasattr(self, '_shock_labels'):
            for lbl in self._shock_labels.values():
                lbl.config(text="---")
        if hasattr(self, '_coal_label'):
            self._coal_label.config(text="---")
            self._coal_warn.config(text="")
        if hasattr(self, '_plume_type_lbl'):
            self._plume_type_lbl.config(text="---", foreground='blue')
            self._plume_pr_lbl.config(text="")