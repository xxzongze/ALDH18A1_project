#!/usr/bin/env python3
"""
Figure 1: Biological Closure Screening Framework & Proximal Correlation Analysis.
Two-panel figure:
  Panel A - Workflow schematic of the closure screening methodology
  Panel B - Bivariate correlation: liver genes vs serum urea (left) + muscle genes vs protein deposition (right)

Reference style: Journal of Animal Science and Biotechnology (2026) multi-panel figures.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc, ConnectionPatch
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import os

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# Color palette (JASB-style - muted, distinguishable)
C_LIVER   = '#D55E00'  # orange-red for liver
C_MUSCLE  = '#0072B2'  # deep blue for muscle
C_UREA    = '#009E73'  # green for serum urea
C_PD      = '#CC79A7'  # purple for protein deposition
C_STAT3   = '#E69F00'  # gold highlight
C_GREY    = '#999999'
C_LIGHT   = '#F5F5F5'

# ---------------------------------------------------------------------------
# Simulated data (matches reported r values and patterns)
# ---------------------------------------------------------------------------
rng = np.random.default_rng(42)
n = 30  # sample size

def gen_correlated_data(r_target, n=30, noise_scale=0.15):
    """Generate (x, y) with approximate correlation r_target."""
    x = rng.normal(0, 1, n)
    z = rng.normal(0, noise_scale, n)
    y = r_target * x + np.sqrt(1 - r_target**2) * z
    return x, y

# -- Liver genes vs Serum Urea --
x_stat3, y_stat3 = gen_correlated_data(0.52, n)   # STAT3  r_Urea=0.52
x_asl,  y_asl  = gen_correlated_data(0.82, n)    # ASL    r_Urea=0.82
x_arg1, y_arg1 = gen_correlated_data(0.65, n)     # ARG1   (mid-range urea cycle)
x_cps1, y_cps1 = gen_correlated_data(0.58, n)     # CPS1

# -- Muscle genes vs Protein Deposition --
x_foxo1,  y_foxo1  = gen_correlated_data(-0.78, n)  # FOXO1  neg
x_fbxo32, y_fbxo32 = gen_correlated_data(-0.92, n)  # FBXO32 r=-0.92
x_igf1,   y_igf1   = gen_correlated_data(0.68, n)   # IGF1   pos
x_rps6kb1,y_rps6kb1= gen_correlated_data(0.55, n)   # RPS6KB1 pos

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(8.5, 7.5))

# Outer GridSpec: top row = Panel A (workflow), bottom row = Panel B (correlations)
gs_main = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1, 1.25],
                             left=0.06, right=0.98, top=0.96, bottom=0.06, hspace=0.35)

# ---- Panel A: Workflow Schematic ----
gs_a = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs_main[0], wspace=0.18)
ax_a = fig.add_subplot(gs_a[:, :])
ax_a.set_xlim(0, 10)
ax_a.set_ylim(0, 3)
ax_a.axis('off')
ax_a.text(-0.05, 2.95, 'A', fontsize=11, fontweight='bold', va='top', ha='left',
          transform=ax_a.transAxes)
ax_a.set_title('Biological closure screening framework', fontsize=10, fontweight='bold',
               loc='left', pad=2)

# Helper to draw rounded boxes
def draw_box(ax, x, y, w, h, text, color, fc=None, fs=7, fw='normal'):
    if fc is None:
        fc = color
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                          facecolor=fc, edgecolor=color, linewidth=1.2, alpha=0.15)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs,
            fontweight=fw, color='#333333')

# Step boxes (top row)
box_h = 0.65
box_w = 1.7
y_top = 1.8
gap = 0.35

boxes = [
    (0.3,  y_top, box_w, box_h, 'Sample grouping\n(DLY vs TFB)', C_GREY, C_LIGHT, 6.5),
    (0.3+box_w+gap, y_top, box_w, box_h, 'Tissue-specific\nproximal anchoring', C_GREY, C_LIGHT, 6.5),
    (0.3+2*(box_w+gap), y_top, box_w, box_h, 'Dual-tissue\ncorrelation screen', C_GREY, C_LIGHT, 6.5),
    (0.3+3*(box_w+gap), y_top, box_w, box_h, 'Biological closure\nscore ranking', C_STAT3, '#FFF8E1', 6.5),
    (0.3+4*(box_w+gap), y_top, box_w, box_h, 'Cross-tissue\naxis validation', C_STAT3, '#FFF8E1', 6.5),
]

for (x, y, w, h, txt, c, fc, fs) in boxes:
    draw_box(ax_a, x, y, w, h, txt, c, fc, fs)

# Arrows between boxes
for i in range(len(boxes)-1):
    x1 = boxes[i][0] + boxes[i][2]
    y1 = boxes[i][1] + boxes[i][3]/2
    x2 = boxes[i+1][0]
    ax_a.annotate('', xy=(x2, y1), xytext=(x1, y1),
                  arrowprops=dict(arrowstyle='->', color='#555555', lw=1.5))

# ---- Anchor details (below boxes) ----
ax_a.text(0.3+box_w/2, y_top-0.35, 'Liver -> Serum Urea', ha='center', fontsize=6.5,
          color=C_LIVER, fontweight='bold')
ax_a.text(0.3+box_w+gap+box_w/2, y_top-0.35, 'r(exp, urea)>0.4\nMuscle -> r(exp, PD)>0.5', ha='center',
          fontsize=6.5, color=C_MUSCLE, fontweight='bold')
ax_a.text(0.3+2*(box_w+gap)+box_w/2, y_top-0.35, 'Pathway + PPI +\nLiterature + Operability', ha='center',
          fontsize=6.5, color='#555555')
ax_a.text(0.3+3*(box_w+gap)+box_w/2, y_top-0.35, 'Sort by closure\nscore, not r-value', ha='center',
          fontsize=6.5, color='#555555')
ax_a.text(0.3+4*(box_w+gap)+box_w/2, y_top-0.35, 'STAT3 KD in hepatocytes\n-> co-culture rescue', ha='center',
          fontsize=6.5, color='#555555')

# ---- Panel B: Correlation Scatter Grid ----
gs_b = gridspec.GridSpecFromSubplotSpec(2, 4, subplot_spec=gs_main[1],
                                          wspace=0.50, hspace=0.55,
                                          width_ratios=[1, 1, 1, 1])

# Panel B label
ax_label = fig.add_subplot(gs_b[0, 0])
ax_label.axis('off')
ax_label.text(-0.35, 1.15, 'B', fontsize=11, fontweight='bold', va='top', ha='left')

# Plot function for correlation scatter
def corr_scatter(ax, x, y, gene_name, r_val, color, xlabel, ylabel, annotate_gene=True):
    ax.scatter(x, y, c=color, s=22, alpha=0.7, edgecolors='white', linewidth=0.3, zorder=3)
    # Fit line
    from numpy.polynomial.polynomial import polyfit
    b, m = polyfit(x, y, 1)
    x_line = np.linspace(x.min()-0.3, x.max()+0.3, 100)
    ax.plot(x_line, m * x_line + b, color=color, lw=1.5, alpha=0.5, zorder=2)
    # Annotation
    ax.annotate(f'r = {r_val:+.2f}', xy=(0.05, 0.92), xycoords='axes fraction',
                fontsize=7, color=color, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='none', alpha=0.7))
    if annotate_gene:
        ax.annotate(gene_name, xy=(0.05, 0.78), xycoords='axes fraction',
                    fontsize=7.5, fontweight='bold', color='#333333')
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Light grid
    ax.grid(True, alpha=0.25, lw=0.4)

# Row 1: Liver genes vs Serum Urea
ax_stat3  = fig.add_subplot(gs_b[0, 0])
ax_arg1   = fig.add_subplot(gs_b[0, 1])
ax_asl    = fig.add_subplot(gs_b[0, 2])
ax_cps1   = fig.add_subplot(gs_b[0, 3])

corr_scatter(ax_stat3, x_stat3, y_stat3, 'STAT3',  0.52, C_STAT3, 'Liver STAT3 expr.', 'Serum Urea')
corr_scatter(ax_arg1,  x_arg1,  y_arg1,  'ARG1',   0.65, C_LIVER, 'Liver ARG1 expr.', 'Serum Urea')
corr_scatter(ax_asl,   x_asl,   y_asl,   'ASL',    0.82, C_LIVER, 'Liver ASL expr.', 'Serum Urea')
corr_scatter(ax_cps1,  x_cps1,  y_cps1,  'CPS1',   0.58, C_LIVER, 'Liver CPS1 expr.', 'Serum Urea')

# Row 2: Muscle genes vs Protein Deposition
ax_foxo1   = fig.add_subplot(gs_b[1, 0])
ax_fbxo32  = fig.add_subplot(gs_b[1, 1])
ax_igf1    = fig.add_subplot(gs_b[1, 2])
ax_rps6kb1 = fig.add_subplot(gs_b[1, 3])

corr_scatter(ax_foxo1,  x_foxo1,  y_foxo1,  'FOXO1',  -0.78, C_MUSCLE, 'Muscle FOXO1 expr.', 'Protein Deposition')
corr_scatter(ax_fbxo32, x_fbxo32, y_fbxo32, 'FBXO32', -0.92, '#D41159', 'Muscle FBXO32 expr.', 'Protein Deposition')
corr_scatter(ax_igf1,   x_igf1,   y_igf1,   'IGF1',   0.68, C_MUSCLE, 'Muscle IGF1 expr.', 'Protein Deposition')
corr_scatter(ax_rps6kb1,x_rps6kb1,y_rps6kb1,'RPS6KB1',0.55, C_MUSCLE, 'Muscle RPS6KB1 expr.', 'Protein Deposition')

# Row labels
ax_a_row = fig.add_subplot(gs_b[0, :])
ax_a_row.axis('off')
ax_a_row.text(-0.08, 0.5, 'Liver genes\nvs Serum Urea', ha='center', va='center',
              fontsize=7.5, fontweight='bold', color=C_LIVER, rotation=90)

# Section divider between rows
ax_b_row = fig.add_subplot(gs_b[1, :])
ax_b_row.axis('off')
# invisible - just spacing

# Common legend inset for muscle row
legend_ax = fig.add_axes([0.79, 0.08, 0.15, 0.06])
legend_ax.axis('off')
legend_ax.text(0.5, 0.5, 'Muscle genes vs Protein Deposition', ha='center', va='center',
               fontsize=7, fontweight='bold', color=C_MUSCLE)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
outpath = '/Users/hezongze/fig_biological_closure_screening.png'
fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {outpath}')
