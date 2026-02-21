"""
Input Panels Module
User interface for input parameters
"""

import tkinter as tk
from tkinter import ttk, messagebox
from config import *


class InputValidator:
    """Validate user inputs"""
    
    @staticmethod
    def validate_float(value, name, min_val=None, max_val=None):
        try:
            val = float(value)
            if min_val is not None and val < min_val:
                raise ValueError(f"{name} must be >= {min_val}")
            if max_val is not None and val > max_val:
                raise ValueError(f"{name} must be <= {max_val}")
            return val
        except ValueError as e:
            raise ValueError(f"Invalid {name}: {str(e)}")
    
    @staticmethod
    def validate_int(value, name, min_val=None, max_val=None):
        try:
            val = int(value)
            if min_val is not None and val < min_val:
                raise ValueError(f"{name} must be >= {min_val}")
            if max_val is not None and val > max_val:
                raise ValueError(f"{name} must be <= {max_val}")
            return val
        except ValueError as e:
            raise ValueError(f"Invalid {name}: {str(e)}")


class InputPanels:
    """Input parameter panels"""
    
    def __init__(self, parent, calculate_callback):
        self.parent = parent
        self.calculate_callback = calculate_callback
        self.validator = InputValidator()
        
        # Scrollable container - WIDER for buttons
        self.canvas = tk.Canvas(parent, width=360, bg='white')
        self.scrollbar = ttk.Scrollbar(parent, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)
        
        self.scrollable_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        
        # Adjust canvas width to frame width
        def _configure_canvas(event):
            self.canvas.itemconfig(self.canvas_window, width=event.width)
        self.canvas.bind("<Configure>", _configure_canvas)
        
        self.canvas_window = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.scrollbar.pack(side=tk.LEFT, fill=tk.Y)
        
        # Main frame inside scrollable area
        self.frame = ttk.LabelFrame(self.scrollable_frame, text="Input Parameters", padding=10)
        self.frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Unit system toggle
        unit_bar = ttk.Frame(self.frame)
        unit_bar.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(unit_bar, text="Units:", font=('Arial', 9, 'bold')).pack(side=tk.LEFT)
        self._unit_btn = ttk.Button(unit_bar, text="Imperial (psia / °R / in)",
                                   command=self._toggle_units, width=26)
        self._unit_btn.pack(side=tk.LEFT, padx=5)
        self._unit_system = 'imperial'

        self.create_flow_conditions()
        self.create_design_parameters()
        self.create_view_options()
        self.create_freestream_options()
        self.create_plume_options()
        self.create_viscous_options()
        self.create_grid_options()
        self.create_external_flow_options()
        self.create_given_shape_options()
        self.create_buttons()
        
        # Mouse wheel scrolling - IMPROVED
        def _on_mousewheel(event):
            if event.num == 5 or event.delta < 0:
                self.canvas.yview_scroll(1, "units")
            elif event.num == 4 or event.delta > 0:
                self.canvas.yview_scroll(-1, "units")
        
        # Bind to canvas and all children
        self.canvas.bind_all("<MouseWheel>", _on_mousewheel)
        self.canvas.bind_all("<Button-4>", _on_mousewheel)
        self.canvas.bind_all("<Button-5>", _on_mousewheel)
    
    def _toggle_units(self):
        """Toggle between Imperial and SI units."""
        from utils.helpers import (set_unit_system, unit_label,
                                   pressure_to_display, pressure_from_display,
                                   temp_to_display, temp_from_display,
                                   length_to_display, length_from_display)
        try:
            p0_val  = float(self.p0_var.get())
            T0_val  = float(self.T0_var.get())
            th_val  = float(self.throat_var.get())
        except ValueError:
            return

        if self._unit_system == 'imperial':
            # Switch to SI
            set_unit_system('si')
            self._unit_system = 'si'
            self.p0_var.set(f"{pressure_to_display(p0_val):.2f}")
            self.T0_var.set(f"{temp_to_display(T0_val):.1f}")
            self.throat_var.set(f"{length_to_display(th_val):.3f}")
            self._unit_btn.config(text="SI (kPa / K / cm)")
        else:
            # Switch to Imperial
            set_unit_system('imperial')
            self._unit_system = 'imperial'
            self.p0_var.set(f"{pressure_from_display(float(self.p0_var.get())):.2f}")
            self.T0_var.set(f"{temp_from_display(float(self.T0_var.get())):.1f}")
            self.throat_var.set(f"{length_from_display(float(self.throat_var.get())):.3f}")
            self._unit_btn.config(text="Imperial (psia / °R / in)")

    def create_flow_conditions(self):
        frame = ttk.LabelFrame(self.frame, text="Flow Conditions", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Rocket Engine Preset Selector
        ttk.Label(frame, text="Select Rocket Engine:", font=('Arial', 9, 'bold')).grid(row=0, column=0, sticky=tk.W, pady=(0,5))
        
        self.preset_var = tk.StringVar(value="Custom")
        preset_combo = ttk.Combobox(frame, textvariable=self.preset_var, 
                                    values=list(ROCKET_PRESETS.keys()),
                                    state='readonly', width=20)
        preset_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=(0,5))
        preset_combo.bind('<<ComboboxSelected>>', self.on_preset_change)
        
        self.preset_desc_label = ttk.Label(frame, text="Custom values (user input)", 
                                          font=('Arial', 7), foreground='blue')
        self.preset_desc_label.grid(row=1, column=0, columnspan=2, sticky=tk.W, pady=(0,10))
        
        ttk.Separator(frame, orient='horizontal').grid(row=2, column=0, columnspan=2, sticky='ew', pady=5)
        
        ttk.Label(frame, text="Total Pressure (psia):").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.p0_var = tk.StringVar(value=str(DEFAULT_TOTAL_PRESSURE))
        self.p0_entry = ttk.Entry(frame, textvariable=self.p0_var, width=15)
        self.p0_entry.grid(row=3, column=1, padx=5, pady=2)
        
        ttk.Label(frame, text="Total Temp (°R):").grid(row=4, column=0, sticky=tk.W, pady=2)
        self.T0_var = tk.StringVar(value=str(DEFAULT_TOTAL_TEMP))
        self.T0_entry = ttk.Entry(frame, textvariable=self.T0_var, width=15)
        self.T0_entry.grid(row=4, column=1, padx=5, pady=2)
        
        ttk.Label(frame, text="Gamma (γ):").grid(row=5, column=0, sticky=tk.W, pady=2)
        self.gamma_var = tk.StringVar(value=str(DEFAULT_GAMMA))
        self.gamma_var.trace_add('write', self.on_gamma_change)
        self.gamma_entry = ttk.Entry(frame, textvariable=self.gamma_var, width=15)
        self.gamma_entry.grid(row=5, column=1, padx=5, pady=2)
    
    def on_preset_change(self, event=None):
        """Handle rocket engine preset selection"""
        preset_name = self.preset_var.get()
        preset = ROCKET_PRESETS.get(preset_name, ROCKET_PRESETS["Custom"])
        
        self.preset_desc_label.config(text=preset["description"])
        self.p0_var.set(str(preset["p0"]))
        self.T0_var.set(str(preset["T0"]))
        self.mach_exit_var.set(str(preset["mach_exit"]))
        self.throat_var.set(str(preset["throat_radius"]))
        self.problem_var.set(preset["nozzle_type"])
        
        if preset_name != "Custom":
            self.preset_desc_label.config(foreground='green')
        else:
            self.preset_desc_label.config(foreground='blue')
    
    def create_design_parameters(self):
        frame = ttk.LabelFrame(self.frame, text="Design Parameters", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(frame, text="Exit Mach Number:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.mach_exit_var = tk.StringVar(value=str(DEFAULT_MACH_EXIT))
        self.mach_exit_entry = ttk.Entry(frame, textvariable=self.mach_exit_var, width=15)
        self.mach_exit_entry.grid(row=0, column=1, padx=5, pady=2)
        
        ttk.Label(frame, text="Number of Rays:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.num_rays_var = tk.StringVar(value=str(DEFAULT_NUM_RAYS))
        self.num_rays_entry = ttk.Entry(frame, textvariable=self.num_rays_var, width=15)
        self.num_rays_entry.grid(row=1, column=1, padx=5, pady=2)
        
        ttk.Label(frame, text="Throat Height (in):").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.throat_var = tk.StringVar(value=str(DEFAULT_THROAT_HEIGHT))
        self.throat_entry = ttk.Entry(frame, textvariable=self.throat_var, width=15)
        self.throat_entry.grid(row=2, column=1, padx=5, pady=2)
        
        ttk.Label(frame, text="Problem Type:").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.problem_var = tk.IntVar(value=1)
        problem_combo = ttk.Combobox(frame, textvariable=self.problem_var, 
                                     values=[1, 2, 3, 4, 5, 6, 7, 8, 9],
                                     state='readonly', width=5)
        problem_combo.grid(row=3, column=1, sticky=tk.W, padx=5, pady=2)
        
        ttk.Label(frame, text="1=2D, 2=Axisym, 3=Plug, 4=Cone, 5=Wedge, 6=Bell, 7=Trunc, 8=MinLen, 9=ThroatExp", font=('Arial', 6)).grid(row=3, column=2, sticky=tk.W)

        # PERINT - Internal/External expansion ratio (C.C. Lee plug nozzle, prob=5)
        ttk.Label(frame, text="Int/Ext Expansion Ratio (perint):").grid(row=4, column=0, sticky=tk.W, pady=2)
        self.perint_var = tk.DoubleVar(value=0.5)
        perint_frame = ttk.Frame(frame)
        perint_frame.grid(row=4, column=1, columnspan=2, sticky=tk.W, pady=2)
        perint_slider = ttk.Scale(perint_frame, from_=0.01, to=1.0, orient=tk.HORIZONTAL,
                                  variable=self.perint_var, length=150)
        perint_slider.pack(side=tk.LEFT)
        self.perint_label = ttk.Label(perint_frame, text="0.50", width=5)
        self.perint_label.pack(side=tk.LEFT, padx=3)
        ttk.Label(perint_frame, text="(C.C. Lee plug only)", font=("Arial", 7)).pack(side=tk.LEFT)
        perint_slider.config(command=lambda v: self.perint_label.config(text=f"{float(v):.2f}"))

        # CALCULATION METHOD 
        ttk.Label(frame, text="Calculation Method:").grid(row=5, column=0, sticky=tk.W, pady=2)
        self.zone_method_var = tk.IntVar(value=0)
        method_frame = ttk.Frame(frame)
        method_frame.grid(row=5, column=1, sticky=tk.W, pady=2)
        ttk.Radiobutton(method_frame, text="Points", variable=self.zone_method_var, value=0).pack(side=tk.LEFT)
        ttk.Radiobutton(method_frame, text="Zone", variable=self.zone_method_var, value=1).pack(side=tk.LEFT)
    
    def create_view_options(self):
        frame = ttk.LabelFrame(self.frame, text="View Options", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.view_mode_var = tk.StringVar(value=DISPLAY_GEOMETRY)
        
        ttk.Radiobutton(frame, text="Geometry", variable=self.view_mode_var, 
                       value=DISPLAY_GEOMETRY).pack(anchor=tk.W)
        ttk.Radiobutton(frame, text="Mach Contour", variable=self.view_mode_var, 
                       value=DISPLAY_MACH).pack(anchor=tk.W)
        ttk.Radiobutton(frame, text="Pressure Contour", variable=self.view_mode_var, 
                       value=DISPLAY_PRESSURE).pack(anchor=tk.W)
        ttk.Radiobutton(frame, text="Temperature Contour", variable=self.view_mode_var, 
                       value=DISPLAY_TEMP).pack(anchor=tk.W)
        
        self.mesh_var = tk.IntVar(value=1)
        ttk.Checkbutton(frame, text="Show Characteristic Mesh", 
                       variable=self.mesh_var).pack(anchor=tk.W, pady=5)
        
        self.reflection_var = tk.IntVar(value=0)
        ttk.Checkbutton(frame, text="Show Reflection", 
               variable=self.reflection_var,
               command=self.on_reflection_toggle).pack(anchor=tk.W)
         # COLOR BAR CONTROLS 
        ttk.Separator(frame, orient='horizontal').pack(fill=tk.X, pady=5)
        ttk.Label(frame, text="Contour Plot Controls:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        
        ttk.Label(frame, text="Color Min:").pack(anchor=tk.W)
        self.colorbar_min = tk.DoubleVar(value=0.0)
        ttk.Scale(frame, from_=0.0, to=5.0, variable=self.colorbar_min, 
                 orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5)
        
        ttk.Label(frame, text="Color Max:").pack(anchor=tk.W)
        self.colorbar_max = tk.DoubleVar(value=3.0)
        ttk.Scale(frame, from_=1.0, to=10.0, variable=self.colorbar_max, 
                 orient=tk.HORIZONTAL).pack(fill=tk.X, padx=5)
    def create_freestream_options(self):
        """Free stream conditions input"""
        print("DEBUG: Creating Free Stream panel!")

        frame = ttk.LabelFrame(self.frame, text="Free Stream Conditions", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(frame, text="Altitude (ft):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.altitude_var = tk.DoubleVar(value=0.0)
        ttk.Entry(frame, textvariable=self.altitude_var, width=10).grid(row=0, column=1, sticky=tk.W, padx=5, pady=2)
        
        ttk.Label(frame, text="Free Stream Mach:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.fs_mach_var = tk.DoubleVar(value=0.0)
        ttk.Entry(frame, textvariable=self.fs_mach_var, width=10).grid(row=1, column=1, sticky=tk.W, padx=5, pady=2)
        
        ttk.Button(frame, text="Calculate Free Stream", 
                  command=self.on_freestream_calc).grid(row=2, column=0, columnspan=2, pady=5)
    
    
    def on_freestream_calc(self):
        """Calculate free stream properties"""
        from tkinter import messagebox
        import math
        
        alt = self.altitude_var.get()
        mach = self.fs_mach_var.get()
        
        # Standard atmosphere calculation
        if alt <= 36089:  # Troposphere
            T = 518.67 - 0.00356616 * alt  # Rankine
            p = 14.696 * (T / 518.67) ** 5.2561  # psi
        else:  # Stratosphere
            T = 389.97  # Rankine
            p = 14.696 * 0.2234 * math.exp((36089 - alt) / 20806)
        
        # Convert to more useful units
        T_F = T - 459.67  # Fahrenheit
        
        # Dynamic pressure
        gamma = 1.4
        q = 0.5 * gamma * p * mach**2  # psi
        
        # Total conditions
        p0 = p * (1 + (gamma-1)/2 * mach**2) ** (gamma/(gamma-1))
        T0 = T * (1 + (gamma-1)/2 * mach**2)
        T0_F = T0 - 459.67
        
        result = f"""FREE STREAM PROPERTIES:
        
Altitude: {alt:,.0f} ft
Mach Number: {mach:.2f}

Static Conditions:
  Pressure: {p:.4f} psi
  Temperature: {T_F:.1f} °F

Total Conditions:
  p0: {p0:.4f} psi
  T0: {T0_F:.1f} °F

Dynamic Pressure: {q:.4f} psi"""
        
        messagebox.showinfo("Free Stream Results", result)         
    
    def create_plume_options(self):
        frame = ttk.LabelFrame(self.frame, text="Plume Analysis", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.plume_var = tk.IntVar(value=0)
        ttk.Checkbutton(frame, text="Enable Plume Calculation", 
                       variable=self.plume_var).pack(anchor=tk.W)
        
        # PLUME MODE 
        ttk.Label(frame, text="Plume Mode:").pack(anchor=tk.W, pady=(5,0))
        self.plume_mode_var = tk.IntVar(value=1)
        ttk.Radiobutton(frame, text="0 - No Plume", variable=self.plume_mode_var, value=0).pack(anchor=tk.W, padx=20)
        ttk.Radiobutton(frame, text="1 - Simple", variable=self.plume_mode_var, value=1).pack(anchor=tk.W, padx=20)
        ttk.Radiobutton(frame, text="2 - Advanced (Shocks)", variable=self.plume_mode_var, value=2).pack(anchor=tk.W, padx=20)
        
        ttk.Label(frame, text="Ambient Pressure (psia):").pack(anchor=tk.W, pady=(5,0))
        self.p_ambient_var = tk.StringVar(value="14.7")
        ttk.Entry(frame, textvariable=self.p_ambient_var, width=15).pack(anchor=tk.W, padx=20)
    
    def create_viscous_options(self):
        frame = ttk.LabelFrame(self.frame, text="Viscous Correction (Boundary Layer)", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)

        self.viscous_var = tk.IntVar(value=0)
        ttk.Checkbutton(frame, text="Enable Viscous Correction",
                       variable=self.viscous_var).pack(anchor=tk.W)

        sub = ttk.Frame(frame)
        sub.pack(fill=tk.X, padx=5, pady=3)

        ttk.Label(sub, text="Re (throat):").grid(row=0, column=0, sticky=tk.W)
        self.re_throat_var = tk.StringVar(value="1e6")
        ttk.Entry(sub, textvariable=self.re_throat_var, width=10).grid(row=0, column=1, padx=5)

        self.axisym_bl_var = tk.IntVar(value=1)
        ttk.Checkbutton(sub, text="Axisymmetric (Mangler corr.)",
                       variable=self.axisym_bl_var).grid(row=1, column=0, columnspan=2, sticky=tk.W, pady=2)

        ttk.Label(frame, text="Corrects wall geometry for BL displacement thickness δ*",
                 font=('Arial', 7), foreground='gray').pack(anchor=tk.W)

    def create_grid_options(self):
        frame = ttk.LabelFrame(self.frame, text="Grid Quality Settings", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(frame, text="Grid Quality Preset:", font=('Arial', 8, 'bold')).grid(row=0, column=0, sticky=tk.W, pady=2)
        
        self.grid_quality_var = tk.StringVar(value="Medium")
        quality_frame = ttk.Frame(frame)
        quality_frame.grid(row=1, column=0, columnspan=2, sticky=tk.W, padx=5)
        
        qualities = [
            ("Coarse (Fast)", "Coarse", 15),
            ("Medium (Balanced)", "Medium", 30),
            ("Fine (Accurate)", "Fine", 50),
            ("Very Fine (Slow)", "VeryFine", 80)
        ]
        
        for i, (text, value, rays) in enumerate(qualities):
            rb = ttk.Radiobutton(quality_frame, text=text, 
                                variable=self.grid_quality_var, 
                                value=value,
                                command=lambda r=rays: self.on_quality_change(r))
            rb.pack(anchor=tk.W, pady=1)
        
        ttk.Label(frame, text="Manual Grid Control:", font=('Arial', 8, 'bold')).grid(row=2, column=0, sticky=tk.W, pady=(10, 2))
        ttk.Label(frame, text="Number of Rays:", font=('Arial', 8)).grid(row=3, column=0, sticky=tk.W, pady=2)
        
        slider_frame = ttk.Frame(frame)
        slider_frame.grid(row=4, column=0, columnspan=2, sticky=tk.EW, padx=5)
        
        self.rays_slider = tk.Scale(slider_frame, from_=10, to=100, 
                                    orient=tk.HORIZONTAL,
                                    variable=self.num_rays_var,
                                    command=self.on_slider_change,
                                    length=200)
        self.rays_slider.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        self.rays_label = ttk.Label(slider_frame, text="30 rays", width=10)
        self.rays_label.pack(side=tk.LEFT, padx=5)
        
        ttk.Label(frame, text="Grid Statistics:", font=('Arial', 8, 'bold')).grid(row=5, column=0, sticky=tk.W, pady=(10, 2))
        
        self.grid_stats_label = ttk.Label(frame, text="Points: ~465 | Size: 31x31", 
                                         font=('Arial', 7), foreground='blue')
        self.grid_stats_label.grid(row=6, column=0, columnspan=2, sticky=tk.W, padx=5)
        
        self.comp_time_label = ttk.Label(frame, text="Est. time: <1 sec", 
                                        font=('Arial', 7), foreground='green')
        self.comp_time_label.grid(row=7, column=0, columnspan=2, sticky=tk.W, padx=5)
        
        ttk.Label(frame, text="Expected Quality:", font=('Arial', 8, 'bold')).grid(row=8, column=0, sticky=tk.W, pady=(10, 2))
        
        self.quality_indicator = ttk.Label(frame, text="● Good", 
                                          font=('Arial', 8, 'bold'), foreground='green')
        self.quality_indicator.grid(row=9, column=0, sticky=tk.W, padx=5)
    
    def on_quality_change(self, rays_value):
        self.num_rays_var.set(str(rays_value))
        self.update_grid_stats()
    
    def on_slider_change(self, value):
        rays = int(float(value))
        self.rays_label.config(text=f"{rays} rays")
        self.update_grid_stats()
    
    def update_grid_stats(self):
        try:
            rays = int(self.num_rays_var.get())
        except:
            rays = 30
        
        n = rays // 2 + 1
        total_points = (n * (n + 1)) // 2
        grid_size = n
        
        self.grid_stats_label.config(text=f"Points: ~{total_points} | Size: {grid_size}x{grid_size}")
        
        time_estimate = (rays ** 1.5) / 1000.0
        if time_estimate < 1.0:
            time_str = f"Est. time: <1 sec"
            time_color = 'green'
        elif time_estimate < 5.0:
            time_str = f"Est. time: ~{time_estimate:.1f} sec"
            time_color = 'orange'
        else:
            time_str = f"Est. time: ~{time_estimate:.1f} sec"
            time_color = 'red'
        
        self.comp_time_label.config(text=time_str, foreground=time_color)
        
        if rays < 20:
            quality_text = "● Low (Coarse mesh)"
            quality_color = 'red'
        elif rays < 40:
            quality_text = "● Good (Standard)"
            quality_color = 'green'
        elif rays < 70:
            quality_text = "● Excellent (Fine mesh)"
            quality_color = 'blue'
        else:
            quality_text = "● Very High (May be slow)"
            quality_color = 'purple'
        
        self.quality_indicator.config(text=quality_text, foreground=quality_color)
    
    def create_external_flow_options(self):
        frame = ttk.LabelFrame(self.frame, text="External Flow Analysis", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.external_flow_var = tk.IntVar(value=0)
        
        ttk.Checkbutton(frame, text="Enable External Flow Calculations", 
                       variable=self.external_flow_var,
                       onvalue=1, offvalue=0,
                       command=self.on_external_flow_toggle).pack(anchor=tk.W, pady=2)
        
        self.external_options_frame = ttk.Frame(frame)
        self.external_options_frame.pack(fill=tk.X, padx=20, pady=5)
        
        ttk.Label(self.external_options_frame, text="Ambient Pressure (psia):", 
                 font=('Arial', 8)).grid(row=0, column=0, sticky=tk.W, pady=2)
        self.ext_p_ambient_var = tk.StringVar(value="14.7")
        ttk.Entry(self.external_options_frame, textvariable=self.ext_p_ambient_var, 
                 width=10).grid(row=0, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(self.external_options_frame, text="Num. Streamlines:", 
                 font=('Arial', 8)).grid(row=1, column=0, sticky=tk.W, pady=2)
        self.ext_streamlines_var = tk.StringVar(value="10")
        ttk.Entry(self.external_options_frame, textvariable=self.ext_streamlines_var, 
                 width=10).grid(row=1, column=1, sticky=tk.W, padx=5)
        
        self.cowl_var = tk.IntVar(value=0)
        ttk.Checkbutton(self.external_options_frame, text="Show Cowl Design", 
                       variable=self.cowl_var).grid(row=2, column=0, columnspan=2, 
                                                   sticky=tk.W, pady=2)
        
        self.afterbody_var = tk.IntVar(value=0)
        ttk.Checkbutton(self.external_options_frame, text="Calculate Afterbody Pressure", 
                       variable=self.afterbody_var).grid(row=3, column=0, columnspan=2, 
                                                         sticky=tk.W, pady=2)
        
        self._toggle_external_options_state('disabled')
    def create_given_shape_options(self):
        """Given shape analysis input"""
        frame = ttk.LabelFrame(self.frame, text="Given Shape Analysis", padding=5)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(frame, text="Upload your nozzle contour:", font=('Arial', 9)).pack(anchor=tk.W)
        ttk.Label(frame, text="Format: CSV with x,y columns", font=('Arial', 7), foreground='gray').pack(anchor=tk.W)
        
        button_frame = ttk.Frame(frame)
        button_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(button_frame, text="Browse File...", 
                  command=self.browse_shape_file).pack(side=tk.LEFT, padx=5)
        
        self.shape_file_label = ttk.Label(button_frame, text="No file selected", 
                                         foreground='gray', font=('Arial', 8))
        self.shape_file_label.pack(side=tk.LEFT, padx=5)
        
        self.shape_file_path = None
        
        ttk.Button(frame, text="Analyze Given Shape", 
                  command=self.analyze_given_shape,
                  state='disabled').pack(fill=tk.X, pady=5)
        
        self.analyze_shape_btn = frame.winfo_children()[-1]
    
    def browse_shape_file(self):
        """Browse for shape file"""
        from tkinter import filedialog
        filename = filedialog.askopenfilename(
            title="Select Nozzle Shape File",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            self.shape_file_path = filename
            import os
            self.shape_file_label.config(text=os.path.basename(filename), foreground='blue')
            self.analyze_shape_btn.config(state='normal')
    
    def analyze_given_shape(self):
        """Trigger given shape analysis"""
        if self.shape_file_path and hasattr(self, 'shape_analysis_callback'):
            self.shape_analysis_callback(self.shape_file_path)
    
    def on_external_flow_toggle(self):
        if self.external_flow_var.get() == 1:
            self._toggle_external_options_state('normal')
        else:
            self._toggle_external_options_state('disabled')
    
    def _toggle_external_options_state(self, state):
        for child in self.external_options_frame.winfo_children():
            try:
                child.configure(state=state)
            except:
                pass
    
    def create_buttons(self):
        frame = ttk.Frame(self.frame)
        frame.pack(fill=tk.X, padx=5, pady=10)
        
        self.calculate_btn = ttk.Button(frame, text="Calculate", 
                                       command=self.on_calculate,
                                       style='Accent.TButton')
        self.calculate_btn.pack(fill=tk.X, pady=5)
        
        self.results_btn = ttk.Button(frame, text="Show Results", 
                                     command=self.on_show_results)
        self.results_btn.pack(fill=tk.X, pady=5)
    
    def on_calculate(self):
        try:
            params = self.get_parameters()
            self.calculate_callback(params)
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
    
    def on_show_results(self):
        if hasattr(self, 'results_window_callback'):
            self.results_window_callback()
    
    def get_parameters(self):
        """CRITICAL: Returns parameters with correct key names"""
        return {
            'p0': self.validator.validate_float(self.p0_var.get(), "Total Pressure", 0.1, 10000),
            'T0': self.validator.validate_float(self.T0_var.get(), "Total Temperature", 100, 10000),
            'gamma': self.validator.validate_float(self.gamma_var.get(), "Gamma", 1.1, 1.67),
            'mach_exit': self.validator.validate_float(self.mach_exit_var.get(), "Exit Mach", 1.01, 10),
            'num_rays': self.validator.validate_int(self.num_rays_var.get(), "Number of Rays", 5, 100),
            'throat_height': self.validator.validate_float(self.throat_var.get(), "Throat Height", 0.01, 100),
            'problem_type': self.problem_var.get(),
            'view_mode': self.view_mode_var.get(),
            'show_mesh': bool(self.mesh_var.get()),
            'show_reflection': bool(self.reflection_var.get()),
            'plume_enabled': bool(self.plume_var.get()),
            'p_ambient': self.validator.validate_float(self.p_ambient_var.get(), "Ambient Pressure", 0.01, 100),
            'plume_mode': self.plume_mode_var.get(),
            'enable_external_flow': bool(self.external_flow_var.get()),
            'ext_p_ambient': float(self.ext_p_ambient_var.get()),
            'ext_streamlines': int(self.ext_streamlines_var.get()),
            'show_cowl': bool(self.cowl_var.get()),
            'calc_afterbody': bool(self.afterbody_var.get()),
            'colorbar_min': self.colorbar_min.get(),  
            'colorbar_max': self.colorbar_max.get(),
            'use_zone_method': bool(self.zone_method_var.get()),
            'perint': self.perint_var.get(),
            'enable_viscous': bool(self.viscous_var.get()),
            're_throat': float(self.re_throat_var.get()),
            'axisym_bl': bool(self.axisym_bl_var.get())
        }
    
    def reset_inputs(self):
        """Reset all inputs to defaults"""
        self.preset_var.set("Custom")
        self.p0_var.set(str(DEFAULT_TOTAL_PRESSURE))
        self.T0_var.set(str(DEFAULT_TOTAL_TEMP))
        self.gamma_var.set(str(DEFAULT_GAMMA))
        self.mach_exit_var.set(str(DEFAULT_MACH_EXIT))
        self.num_rays_var.set(str(DEFAULT_NUM_RAYS))
        self.throat_var.set(str(DEFAULT_THROAT_HEIGHT))
        self.problem_var.set(DEFAULT_PROBLEM)
        self.view_mode_var.set(DISPLAY_GEOMETRY)
        self.mesh_var.set(1)
        self.reflection_var.set(0)
        self.plume_var.set(0)
        self.p_ambient_var.set("14.7")
        self.external_flow_var.set(0)
        self.on_external_flow_toggle()
        self.update_grid_stats()
    def on_reflection_toggle(self):
        """Reflection checkbox toggle"""
        if hasattr(self, 'canvas_ref'):
            self.canvas_ref.show_reflection = bool(self.reflection_var.get())
            if self.canvas_ref.current_data:
                self.canvas_ref.plot_nozzle_geometry(self.canvas_ref.current_data, self.canvas_ref.show_mesh)
    def on_gamma_change(self, *args):
        """Gamma value changed"""
        #  recalculation 
        try:
            gamma = float(self.gamma_var.get())
            if 1.1 <= gamma <= 1.67:
                self.gamma_entry.config(foreground='blue')
            else:
                self.gamma_entry.config(foreground='red')
        except:
            self.gamma_entry.config(foreground='red')