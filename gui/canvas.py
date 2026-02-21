"""
Canvas Module - Matplotlib visualization
"""

import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
from config import *


class NozzleCanvas:
    """Matplotlib canvas for nozzle visualization"""
    
    def __init__(self, parent):
        self.parent = parent
        self.fig = Figure(figsize=(14, 10), dpi=100, facecolor='white')
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=False, side=tk.TOP)
        self.toolbar = NavigationToolbar2Tk(self.canvas, parent)
        self.toolbar.update()
        self.display_mode = DISPLAY_GEOMETRY
        self.show_mesh = True
        self.show_reflection = False
        self.current_data = None
        self.colorbar = None
        self.setup_plot()
    
    def setup_plot(self):
        self.ax.clear()
        self.ax.set_xlabel('Axial Distance (x)')
        self.ax.set_ylabel('Radial Distance (y)')
        self.ax.set_title('Aircraft Nozzle Interactive')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_aspect('equal')
        self.canvas.draw()
    
    def plot_nozzle_geometry(self, results, show_characteristics=True):
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.colorbar = None
        self.current_data = results
        
        nozzle_type = results.get('nozzle_type', 'standard')
        
        # Plot nozzle walls based on type
        if nozzle_type == 'plug':
            spike_x = results['spike_x']
            spike_y = results['spike_y']
            outer_x = results['outer_x']
            outer_y = results['outer_y']
            
            self.ax.plot(spike_x, spike_y, 'r-', linewidth=2.5, label='Spike Surface')
            self.ax.plot(outer_x, outer_y, 'b-', linewidth=2, label='Outer Boundary')
            self.ax.plot([spike_x[0], outer_x[0]], [spike_y[0], outer_y[0]], 
                       'k--', linewidth=1, alpha=0.5)
            
            if self.show_reflection:
                self.ax.plot(spike_x, -spike_y, 'r-', linewidth=2.5)
                self.ax.plot(outer_x, -outer_y, 'b-', linewidth=2)
                self.ax.plot([spike_x[0], outer_x[0]], [-spike_y[0], -outer_y[0]], 
                           'k--', linewidth=1, alpha=0.5)
        
        elif nozzle_type == 'cone':
            wall_x = results['wall_x']
            wall_y = results['wall_y']
            cone_angle = results.get('cone_half_angle', 15.0)
            self.ax.plot(wall_x, wall_y, 'b-', linewidth=2, 
                       label=f'Cone Wall: Straight {cone_angle}° | 3D axisymmetric | Simple geometry')
            self.ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, label='Centerline')
            
            if self.show_reflection:
                self.ax.plot(wall_x, -wall_y, 'b-', linewidth=2)
        
        elif nozzle_type == 'wedge':
            wall_x = results['wall_x']
            wall_y = results['wall_y']
            wedge_angle = results.get('wedge_angle', 12.0)
            self.ax.plot(wall_x, wall_y, 'b-', linewidth=2, 
                       label=f'Wedge Wall: Straight {wedge_angle}° | 2D planar | Scramjet use')
            self.ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, label='Centerline')
            
            if self.show_reflection:
                self.ax.plot(wall_x, -wall_y, 'b-', linewidth=2)
        
        elif nozzle_type == 'bell':
            wall_x = results['wall_x']
            wall_y = results['wall_y']
            theta_exit = results.get('theta_exit', 10.0)
            self.ax.plot(wall_x, wall_y, 'b-', linewidth=2, 
                       label=f'Bell Wall: Contoured parabolic | Exit angle {theta_exit}° | Most common rocket')
            self.ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, label='Centerline')
            
            if self.show_reflection:
                self.ax.plot(wall_x, -wall_y, 'b-', linewidth=2)
        
        elif nozzle_type == 'truncated':
            wall_x = results['wall_x']
            wall_y = results['wall_y']
            trunc_ratio = results.get('truncation_ratio', 0.75)
            self.ax.plot(wall_x, wall_y, 'b-', linewidth=2, 
                       label=f'Truncated Wall: Cut at {trunc_ratio*100:.0f}% | ~97-98% efficient | Practical length')
            self.ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, label='Centerline')
            
            if self.show_reflection:
                self.ax.plot(wall_x, -wall_y, 'b-', linewidth=2)
        
        elif nozzle_type == 'minlength':
            wall_x = results['wall_x']
            wall_y = results['wall_y']
            self.ax.plot(wall_x, wall_y, 'b-', linewidth=2, 
                       label=f'MinLength Wall: Shortest possible | Rao optimization | Max weight savings')
            self.ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, label='Centerline')
            
            if self.show_reflection:
                self.ax.plot(wall_x, -wall_y, 'b-', linewidth=2)
        
        elif nozzle_type == 'throatexp':
            wall_x = results['wall_x']
            wall_y = results['wall_y']
            exp_angle = results.get('expansion_angle', 45.0)
            self.ax.plot(wall_x, wall_y, 'b-', linewidth=2, 
                       label=f'ThroatExp Wall: Rapid expansion {exp_angle}° | Compact design | Space-constrained')
            self.ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, label='Centerline')
            
            if self.show_reflection:
                self.ax.plot(wall_x, -wall_y, 'b-', linewidth=2)
        
        else:
            wall_x = results['wall_x']
            wall_y = results['wall_y']
            self.ax.plot(wall_x, wall_y, 'b-', linewidth=2, label='Nozzle Wall')
            self.ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, label='Centerline')
            
            if self.show_reflection:
                self.ax.plot(wall_x, -wall_y, 'b-', linewidth=2)
        
        # Characteristic lines
        if show_characteristics and self.show_mesh:
            if 'characteristic_lines' in results:
                char_data = results['characteristic_lines']
                
                for x_line, y_line in char_data['plus_characteristics']:
                    self.ax.plot(x_line, y_line, 'b-', linewidth=0.5, alpha=0.4)
                
                for x_line, y_line in char_data['minus_characteristics']:
                    self.ax.plot(x_line, y_line, 'r-', linewidth=0.5, alpha=0.4)
                
                self.ax.plot([], [], 'b-', linewidth=1, alpha=0.6, label='C+ Characteristics')
                self.ax.plot([], [], 'r-', linewidth=1, alpha=0.6, label='C- Characteristics')
        
        # Plume visualization
        if 'plume_data' in results:
            plume_data = results['plume_data']
            if plume_data and 'plume_x' in plume_data and plume_data['plume_x']:
                for i, (px, py) in enumerate(zip(plume_data['plume_x'], plume_data['plume_y'])):
                    if i == 0:
                        self.ax.plot(px, py, 'orange', linewidth=1.5, 
                                   alpha=0.7, linestyle='--', label='Jet Plume')
                    else:
                        self.ax.plot(px, py, 'orange', linewidth=1, 
                                   alpha=0.6, linestyle='--')

        if self.show_reflection and 'plume_data' in results:
            plume_data = results['plume_data']
            if plume_data and 'plume_x' in plume_data:
                for px, py in zip(plume_data['plume_x'], plume_data['plume_y']):
                    self.ax.plot(px, [-y for y in py], 'orange', linewidth=1, alpha=0.4, linestyle='--')
        
        # Shock diamonds
        if 'plume_data' in results and 'shock_diamonds' in results['plume_data']:
            shock_data = results['plume_data']['shock_diamonds']
            
            if shock_data and shock_data['type'] != 'perfectly_expanded':
                first_diamond = True
                for diamond in shock_data.get('diamonds', []):
                    x_start = diamond['x_start']
                    x_end = diamond['x_end']
                    x_center = diamond['x_center']
                    radius = diamond['radius']
                    
                    label = 'Shock Diamonds' if first_diamond else None
                    first_diamond = False
                    
                    self.ax.plot([x_start, x_center], [0, radius], 
                               'r-', linewidth=1.5, alpha=0.8, label=label)
                    self.ax.plot([x_start, x_center], [0, -radius], 
                               'r-', linewidth=1.5, alpha=0.8)
                    self.ax.plot([x_center, x_end], [radius, 0], 
                               'r-', linewidth=1.5, alpha=0.8)
                    self.ax.plot([x_center, x_end], [-radius, 0], 
                               'r-', linewidth=1.5, alpha=0.8)
        
        # External flow streamlines
        if 'external_flow' in results:
            ext_data = results['external_flow']
            if 'streamlines' in ext_data:
                for i, streamline in enumerate(ext_data['streamlines']):
                    if i == 0:
                        self.ax.plot(streamline['x'], streamline['y'], 
                                   'cyan', linewidth=1, alpha=0.5, 
                                   linestyle=':', label='External Flow')
                    else:
                        self.ax.plot(streamline['x'], streamline['y'], 
                                   'cyan', linewidth=1, alpha=0.5, linestyle=':')
        
        # Cowl
        if 'cowl' in results:
            cowl_data = results['cowl']
            self.ax.plot(cowl_data['x'], cowl_data['y'], 
                       'purple', linewidth=2, alpha=0.6, label='Cowl')

        # Shock waves
        if 'shocks' in results and len(results['shocks']) > 0:
            for i, shock in enumerate(results['shocks'][:3]):
                x, y = shock['x'], shock['y']
                self.ax.plot(x, y, 'rx', markersize=8, markeredgewidth=2, 
                           label='Shock' if i == 0 else '')
                if i < 2:
                    self.ax.text(x+0.2, y+0.2, f"M1={shock['M1']:.2f}\nM2={shock['M2']:.2f}", 
                               fontsize=7, color='red',
                               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

        if self.show_reflection and 'shocks' in results and len(results['shocks']) > 0:
            for shock in results['shocks'][:3]:
                x, y = shock['x'], shock['y']
                self.ax.plot(x, -y, 'rx', markersize=8, markeredgewidth=2, alpha=0.5)

        # Expansion fan
        expansion_data = results.get('expansion_fan_data', None)
        if expansion_data and 'lines' in expansion_data:
            for idx, (x_line, y_line) in enumerate(expansion_data['lines']):
                label = 'Expansion Fan' if idx == 0 else ''
                self.ax.plot(x_line, y_line, 'g--', linewidth=1.5, alpha=0.7, label=label)

        if self.show_reflection and expansion_data and 'lines' in expansion_data:
            for x_line, y_line in expansion_data['lines']:
                self.ax.plot(x_line, [-y for y in y_line], 'g--', linewidth=1.5, alpha=0.5)

        # Labels and title
        self.ax.set_xlabel('Axial Distance (x)', fontsize=10)
        self.ax.set_ylabel('Radial Distance (y)', fontsize=10)
        
        if nozzle_type == 'plug':
            title = 'Plug Nozzle (Aerospike)'
            desc = 'Altitude compensating • Spike surface + free boundary\nHigh efficiency at all altitudes'
        elif nozzle_type == 'cone':
            cone_angle = results.get('cone_half_angle', 15.0)
            title = f'Cone Nozzle (α={cone_angle}°)'
            desc = f'Conical expansion • Simple geometry\nEfficiency ~95-98% • Easy manufacturing'
        elif nozzle_type == 'wedge':
            wedge_angle = results.get('wedge_angle', 12.0)
            title = f'Wedge Nozzle (θ={wedge_angle}°)'
            desc = f'2D planar expansion at {wedge_angle}° wedge angle\nBlue line: Straight wedge wall\nRed/Blue mesh: Flow characteristics\nUsed in scramjet & hypersonic vehicles'
        elif nozzle_type == 'bell':
            theta_exit = results.get('theta_exit', 10.0)
            title = f'Bell Nozzle (θ_exit={theta_exit}°)'
            desc = f'Contoured wall (parabolic) • Most common rocket nozzle\nBlue curve: Optimized bell contour\n~1-2% more efficient than cone • Shorter & lighter\nUsed in: SpaceX Merlin, Saturn V, Space Shuttle'
        elif nozzle_type == 'truncated':
            trunc_ratio = results.get('truncation_ratio', 0.75)
            title = f'Truncated Ideal Nozzle ({trunc_ratio*100:.0f}% length)'
            desc = f'Ideal nozzle cut short at {trunc_ratio*100:.0f}% of full length\nBlue curve: Truncated ideal contour\nPractical compromise: ~97-98% efficiency\nShorter = lighter • Common in real applications'
        elif nozzle_type == 'minlength':
            title = 'Minimum Length Nozzle'
            desc = f'Shortest possible nozzle (Rao method)\nBlue curve: Optimized minimum length contour\nRapid expansion then straightening\n~95-97% efficiency • Maximum weight savings'
        elif nozzle_type == 'throatexp':
            exp_angle = results.get('expansion_angle', 45.0)
            title = f'Throat Expansion Nozzle (θ={exp_angle}°)'
            desc = f'Rapid expansion immediately at throat ({exp_angle}°)\nBlue curve: Throat expansion profile\nCompact design for space-constrained applications\nSpecial geometric requirements'
        else:
            title = 'Nozzle Geometry'
            desc = 'Method of Characteristics design'
        
        self.ax.set_title(title, fontsize=12, fontweight='bold')
        self.fig.text(0.15, 0.96, desc, fontsize=9, 
                     ha='left', va='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
        
        self.ax.grid(True, alpha=0.3)
        self.ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9, frameon=True)
        self.ax.set_aspect('equal')
        
        # Auto-scale
        if nozzle_type == 'plug':
            all_x = np.concatenate([results['spike_x'], results['outer_x']])
            all_y = np.concatenate([results['spike_y'], results['outer_y']])
            x_min, x_max = np.min(all_x), np.max(all_x)
            y_min, y_max = np.min(all_y), np.max(all_y)
        else:
            wall_x_data = results['wall_x']
            wall_y_data = results['wall_y']
            x_min, x_max = np.min(wall_x_data), np.max(wall_x_data)
            y_min, y_max = -np.max(wall_y_data), np.max(wall_y_data)
        
        padding = 0.1
        x_range = x_max - x_min
        y_range = y_max - y_min
        
        self.ax.set_xlim([x_min - padding * x_range, x_max + padding * x_range])
        self.ax.set_ylim([y_min - padding * y_range, y_max + padding * y_range])
        
        self.canvas.draw()
    
    def plot_contour(self, X, Y, Z, variable='mach'):
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        
        if variable == 'mach':
            cmap = COLORMAP_MACH
            label = 'Mach Number'
        elif variable == 'pressure':
            cmap = COLORMAP_PRESSURE
            label = 'Pressure Ratio (p/p_t)'
        elif variable == 'temperature':
            cmap = COLORMAP_TEMPERATURE
            label = 'Temperature Ratio (T/T_t)'
        else:
            cmap = 'viridis'
            label = variable
        
        levels = NUM_CONTOUR_LEVELS
        vmin = getattr(self, 'colorbar_min', None)
        vmax = getattr(self, 'colorbar_max', None)
        
        contour = self.ax.contourf(X, Y, Z, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
        
        if self.colorbar is not None:
            try:
                self.colorbar.remove()
            except:
                pass
            self.colorbar = None
        
        cbar = self.fig.colorbar(contour, ax=self.ax)
        self.colorbar = cbar
        cbar.set_label(label, fontsize=10)
        
        if self.current_data is not None:
            nozzle_type = self.current_data.get('nozzle_type', 'standard')
            
            if nozzle_type == 'plug':
                spike_x = self.current_data['spike_x']
                spike_y = self.current_data['spike_y']
                outer_x = self.current_data['outer_x']
                outer_y = self.current_data['outer_y']
                
                self.ax.plot(spike_x, spike_y, 'k-', linewidth=2)
                self.ax.plot(outer_x, outer_y, 'k-', linewidth=2)
                
                if self.show_reflection:
                    self.ax.plot(spike_x, -spike_y, 'k-', linewidth=2)
                    self.ax.plot(outer_x, -outer_y, 'k-', linewidth=2)
            else:
                wall_x = self.current_data['wall_x']
                wall_y = self.current_data['wall_y']
                self.ax.plot(wall_x, wall_y, 'k-', linewidth=2)
                
                if self.show_reflection:
                    self.ax.plot(wall_x, -wall_y, 'k-', linewidth=2)
        
        self.ax.set_xlabel('Axial Distance (x)', fontsize=10)
        self.ax.set_ylabel('Radial Distance (y)', fontsize=10)
        self.ax.set_title(f'{label} Contour', fontsize=12, fontweight='bold')
        self.ax.set_aspect('equal')
        
        x_min, x_max = np.min(X), np.max(X)
        y_max = np.max(np.abs(Y))
        padding = 0.02
        x_range = x_max - x_min
        
        self.ax.set_xlim([x_min - padding * x_range, x_max + padding * x_range])
        self.ax.set_ylim([-y_max * (1 + padding), y_max * (1 + padding)])
        
        self.canvas.draw()
    
    def toggle_reflection(self):
        self.show_reflection = not self.show_reflection
        if self.current_data is not None:
            self.plot_nozzle_geometry(self.current_data, self.show_mesh)
    
    def toggle_mesh(self):
        self.show_mesh = not self.show_mesh
        if self.current_data is not None:
            self.plot_nozzle_geometry(self.current_data, self.show_mesh)
    
    def clear(self):
        self.ax.clear()
        self.setup_plot()
    
    def save_figure(self, filename):
        self.fig.savefig(filename, dpi=EXPORT_IMAGE_DPI, bbox_inches='tight')