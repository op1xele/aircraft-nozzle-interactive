"""
PDF Report Generator
Pickle Brothers - Aircraft Nozzle Interactive
Usage: from core.pdf_report import generate_pdf
"""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, HRFlowable, Image)
from reportlab.lib.enums import TA_CENTER, TA_LEFT
import io
import datetime


def generate_pdf(filepath, params, results, project_name="Nozzle Design",
                 engineer="", revision="A"):
    """
    Generate a professional PDF report.
    
    Args:
        filepath   : output .pdf path
        params     : dict from input_panels.get_parameters()
        results    : dict from MOC calculation (current_results)
        project_name, engineer, revision: metadata
    """
    doc = SimpleDocTemplate(
        filepath, pagesize=A4,
        leftMargin=2*cm, rightMargin=2*cm,
        topMargin=2*cm, bottomMargin=2*cm
    )

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle('Title', parent=styles['Title'],
                                 fontSize=18, textColor=colors.HexColor('#1A237E'),
                                 spaceAfter=6)
    h1 = ParagraphStyle('H1', parent=styles['Heading1'],
                        fontSize=13, textColor=colors.HexColor('#1A237E'),
                        spaceBefore=12, spaceAfter=4)
    h2 = ParagraphStyle('H2', parent=styles['Heading2'],
                        fontSize=11, textColor=colors.HexColor('#37474F'),
                        spaceBefore=8, spaceAfter=2)
    normal = styles['Normal']
    small  = ParagraphStyle('Small', parent=normal, fontSize=8,
                            textColor=colors.grey)

    story = []

    # ---- Header ----
    story.append(Paragraph("Pickle Brothers", small))
    story.append(Paragraph("Aircraft Nozzle Interactive", small))
    story.append(Spacer(1, 0.3*cm))
    story.append(HRFlowable(width="100%", thickness=2,
                            color=colors.HexColor('#1A237E')))
    story.append(Spacer(1, 0.3*cm))
    story.append(Paragraph(project_name, title_style))

    date_str = datetime.datetime.now().strftime("%Y-%m-%d  %H:%M")
    meta_data = [
        ['Date:', date_str,    'Engineer:', engineer or '---'],
        ['Revision:', revision, 'Software:', 'Aircraft Nozzle Interactive v2.0'],
    ]
    meta_table = Table(meta_data, colWidths=[3*cm, 5*cm, 3*cm, 6*cm])
    meta_table.setStyle(TableStyle([
        ('FONTSIZE',    (0,0), (-1,-1), 9),
        ('FONTNAME',    (0,0), (0,-1),  'Helvetica-Bold'),
        ('FONTNAME',    (2,0), (2,-1),  'Helvetica-Bold'),
        ('TEXTCOLOR',   (0,0), (0,-1),  colors.HexColor('#1A237E')),
        ('TEXTCOLOR',   (2,0), (2,-1),  colors.HexColor('#1A237E')),
        ('BOTTOMPADDING',(0,0),(-1,-1), 3),
    ]))
    story.append(meta_table)
    story.append(HRFlowable(width="100%", thickness=1, color=colors.lightgrey))
    story.append(Spacer(1, 0.5*cm))

    # ---- Design Inputs ----
    story.append(Paragraph("1. Design Inputs", h1))

    nozzle_types = {1:'2D Ideal', 2:'Axisymmetric', 3:'Plug (Aerospike)',
                    4:'Cone', 5:'Wedge', 6:'Bell', 7:'Truncated Ideal',
                    8:'Minimum Length', 9:'Throat Expansion'}
    prob = params.get('problem_type', 1)

    input_data = [
        ['Parameter', 'Value', 'Unit'],
        ['Total Pressure (p₀)',   f"{params.get('p0', '---'):.2f}",   'psia'],
        ['Total Temperature (T₀)',f"{params.get('T0', '---'):.1f}",   '°R'],
        ['Gamma (γ)',              f"{params.get('gamma', 1.4):.4f}",  '—'],
        ['Exit Mach Number',      f"{params.get('mach_exit', '---'):.3f}", '—'],
        ['Throat Height',         f"{params.get('throat_height', '---'):.3f}", 'in'],
        ['Number of Rays',        str(params.get('num_rays', 30)),     '—'],
        ['Nozzle Type',           nozzle_types.get(prob, str(prob)),   '—'],
        ['Ambient Pressure',      f"{params.get('p_ambient', 14.7):.3f}", 'psia'],
    ]
    story.append(_make_table(input_data))
    story.append(Spacer(1, 0.4*cm))

    # ---- Geometry ----
    if results:
        story.append(Paragraph("2. Nozzle Geometry", h1))
        wall_y = results.get('wall_y', [])
        wall_x = results.get('wall_x', [])
        throat_h = params.get('throat_height', 1.0)
        exit_h   = wall_y[-1] if wall_y else 0
        length   = (wall_x[-1] - wall_x[0]) if wall_x else 0
        area_ratio = exit_h / throat_h if throat_h > 0 else 0

        geom_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Throat Height',      f"{throat_h:.4f}",    'in'],
            ['Exit Height',        f"{exit_h:.4f}",      'in'],
            ['Area Ratio (Ae/At)', f"{area_ratio:.4f}",  '—'],
            ['Nozzle Length',      f"{length:.4f}",      'in'],
            ['L/D Ratio',          f"{length/throat_h:.3f}" if throat_h else '---', '—'],
            ['Exit Mach',          f"{results.get('exit_mach', 0):.4f}", '—'],
            ['Exit p/pt',          f"{results.get('exit_p_pt', 0):.5f}",  '—'],
            ['Exit T/Tt',          f"{results.get('exit_T_Tt', 0):.5f}",  '—'],
        ]
        story.append(_make_table(geom_data))
        story.append(Spacer(1, 0.4*cm))

        # ---- Performance ----
        perf = results.get('performance', {})
        if perf:
            story.append(Paragraph("3. Performance Metrics", h1))
            perf_data = [
                ['Parameter', 'Value', 'Unit'],
                ['Thrust',           f"{perf.get('thrust', 0):.2f}",       'lbf'],
                ['Mass Flow Rate',   f"{perf.get('mass_flow', 0):.4f}",    'lbm/s'],
                ['Specific Impulse', f"{perf.get('Isp', 0):.2f}",          's'],
                ['C* (C-star)',       f"{perf.get('C_star', 0):.1f}",       'ft/s'],
                ['Thrust Coeff Cf',  f"{perf.get('C_F', 0):.5f}",          '—'],
                ['Nozzle Eff.',      f"{perf.get('efficiency', 0):.2f}",    '%'],
                ['NPR (p₀/pₐ)',      f"{perf.get('NPR', 0):.2f}",          '—'],
                ['Exit Pressure',    f"{perf.get('exit_pressure', 0):.3f}", 'psia'],
                ['Exit Temperature', f"{perf.get('exit_temperature', 0):.1f}", '°R'],
                ['Exit Velocity',    f"{perf.get('exit_velocity', 0):.1f}", 'ft/s'],
                ['Expansion',        perf.get('expansion_condition', '---'),'—'],
            ]
            story.append(_make_table(perf_data))
            story.append(Spacer(1, 0.4*cm))

        # ---- Viscous correction ----
        visc = results.get('viscous_correction', {})
        if visc:
            story.append(Paragraph("4. Viscous Correction", h1))
            visc_data = [
                ['Parameter', 'Value', 'Unit'],
                ['Max δ* (displacement thickness)', f"{visc.get('max_delta_star_in', 0):.5f}", 'in'],
                ['Exit δ*',    f"{visc.get('exit_delta_star_in', 0):.5f}",  'in'],
                ['Exit correction', f"{visc.get('exit_correction_pct', 0):.3f}", '%'],
                ['Mean correction', f"{visc.get('mean_correction_pct', 0):.3f}", '%'],
            ]
            story.append(_make_table(visc_data))
            story.append(Spacer(1, 0.4*cm))

    # ---- Footer ----
    story.append(Spacer(1, 1*cm))
    story.append(HRFlowable(width="100%", thickness=1, color=colors.lightgrey))
    story.append(Paragraph(
        f"Pickle Brothers Engineering — Aircraft Nozzle Interactive v2.0 — {date_str}",
        small))

    doc.build(story)
    return filepath


def _make_table(data):
    """Helper: create styled table from data list."""
    col_widths = [6*cm, 5*cm, 4*cm]
    t = Table(data, colWidths=col_widths)
    t.setStyle(TableStyle([
        # Header row
        ('BACKGROUND',  (0,0), (-1,0),  colors.HexColor('#1A237E')),
        ('TEXTCOLOR',   (0,0), (-1,0),  colors.white),
        ('FONTNAME',    (0,0), (-1,0),  'Helvetica-Bold'),
        ('FONTSIZE',    (0,0), (-1,0),  10),
        # Data rows
        ('FONTSIZE',    (0,1), (-1,-1), 9),
        ('ROWBACKGROUNDS', (0,1), (-1,-1),
         [colors.white, colors.HexColor('#F5F5F5')]),
        ('FONTNAME',    (0,1), (0,-1),  'Helvetica-Bold'),
        ('TEXTCOLOR',   (0,1), (0,-1),  colors.HexColor('#37474F')),
        # Grid
        ('GRID',        (0,0), (-1,-1), 0.5, colors.HexColor('#B0BEC5')),
        ('TOPPADDING',  (0,0), (-1,-1), 4),
        ('BOTTOMPADDING',(0,0),(-1,-1), 4),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
    ]))
    return t