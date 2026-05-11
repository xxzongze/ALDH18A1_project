#!/usr/bin/env python3
"""
Figure 3: Cross-tissue STAT3-centered Nitrogen Allocation Axis.
Mechanism diagram showing the IL6→STAT3→Urea Cycle→Serum Urea→FoxO→UPS→Protein Deposition
biological closure chain, spanning liver and muscle compartments.

Reference style: JASB Graphical Abstract / mechanism figures.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc, Ellipse, Rectangle
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

fig, ax = plt.subplots(1, 1, figsize=(10, 5.5))
ax.set_xlim(0, 12)
ax.set_ylim(0, 6)
ax.axis('off')

# ---- Color scheme ----
C_LIVER_BG    = '#FFF3E0'
C_MUSCLE_BG   = '#E3F2FD'
C_LIVER_BORDER = '#D55E00'
C_MUSCLE_BORDER= '#0072B2'
C_ARROW       = '#555555'
C_HIGHLIGHT   = '#D41159'
C_GOLD        = '#E69F00'
C_UREA_SIGNAL = '#009E73'
C_NODE_EDGE   = '#444444'

# ---- Compartment backgrounds ----
# Liver compartment
liver_box = FancyBboxPatch((0.3, 0.5), 5.0, 5.0, boxstyle="round,pad=0.3",
                            facecolor=C_LIVER_BG, edgecolor=C_LIVER_BORDER,
                            linewidth=1.5, alpha=0.6, zorder=0)
ax.add_patch(liver_box)
ax.text(2.8, 5.25, 'LIVER', ha='center', fontsize=9, fontweight='bold',
        color=C_LIVER_BORDER, zorder=5)
ax.text(2.8, 4.95, 'Proximal anchor: r(exp, Serum Urea)', ha='center',
        fontsize=6, color='#888888', style='italic', zorder=5)

# Muscle compartment
muscle_box = FancyBboxPatch((6.8, 0.5), 4.8, 5.0, boxstyle="round,pad=0.3",
                             facecolor=C_MUSCLE_BG, edgecolor=C_MUSCLE_BORDER,
                             linewidth=1.5, alpha=0.6, zorder=0)
ax.add_patch(muscle_box)
ax.text(9.2, 5.25, 'SKELETAL MUSCLE', ha='center', fontsize=9, fontweight='bold',
        color=C_MUSCLE_BORDER, zorder=5)
ax.text(9.2, 4.95, 'Proximal anchor: r(exp, Protein Deposition)', ha='center',
        fontsize=6, color='#888888', style='italic', zorder=5)

# ---- Helper functions ----
def draw_node(ax, x, y, w, h, label, color, fontsize=7, fontweight='bold',
              textcolor='white', radius=0.12, zorder=10):
    """Draw a rounded gene/protein node."""
    box = FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle=f"round,pad={radius}",
                          facecolor=color, edgecolor=C_NODE_EDGE, linewidth=0.8, zorder=zorder)
    ax.add_patch(box)
    ax.text(x, y, label, ha='center', va='center', fontsize=fontsize,
            fontweight=fontweight, color=textcolor, zorder=zorder+1)

def draw_arrow(ax, x1, y1, x2, y2, color=C_ARROW, lw=1.5, zorder=5, style='simple'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                               connectionstyle='arc3,rad=0'), zorder=zorder)

def draw_label(ax, x, y, text, fontsize=6, color='#666666', ha='center', va='center',
               rotation=0, style='italic'):
    ax.text(x, y, text, ha=ha, va=va, fontsize=fontsize, color=color,
            rotation=rotation, style='italic' if style=='italic' else 'normal')

# ---- Node positions ----
# Liver section nodes
nodes_liver = {
    'IL6':    (1.4, 3.6, '#E69F00'),      # JAK/STAT upstream
    'STAT3':  (3.0, 3.6, '#CC79A7'),      # Master regulator
    'CPS1':   (1.4, 2.1, '#D55E00'),      # Urea cycle 1
    'ARG1':   (2.2, 1.35, '#D55E00'),     # Urea cycle 2
    'ASS1':   (3.0, 1.35, '#D55E00'),     # Urea cycle 3
    'ASL':    (3.8, 1.35, '#D55E00'),     # Urea cycle 4
}

# Draw liver nodes
for name, (x, y, color) in nodes_liver.items():
    w = 0.65 if name not in ('STAT3',) else 0.7
    draw_node(ax, x, y, w, 0.42, name, color, fontsize=6.5, fontweight='bold')

# Muscle section nodes
nodes_muscle = {
    'FOXO1':   (7.8, 3.2, '#0072B2'),
    'FOXO3':   (9.0, 3.2, '#0072B2'),
    'IGF1':    (8.4, 4.2, '#009E73'),
    'RPS6KB1': (9.8, 4.2, '#56B4E9'),
    'FBXO32':  (7.8, 1.8, '#D41159'),
    'TRIM63':  (9.4, 1.8, '#D41159'),
}

for name, (x, y, color) in nodes_muscle.items():
    w = 0.75 if name not in ('RPS6KB1',) else 0.82
    draw_node(ax, x, y, w, 0.42, name, color, fontsize=6.5)

# ---- Arrows: Liver signaling ----
# IL6 → STAT3
draw_arrow(ax, 1.9, 3.6, 2.5, 3.6, C_ARROW, lw=2.0)
draw_label(ax, 2.2, 3.85, 'JAK/STAT\nactivation', fontsize=5.5, color=C_GOLD)

# STAT3 → Urea Cycle genes (branching)
# STAT3 to ARG1
ax.annotate('', xy=(2.45, 1.55), xytext=(2.95, 3.2),
            arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=1.2,
                           connectionstyle='arc3,rad=-0.3'), zorder=5)
draw_label(ax, 2.8, 2.35, '4 binding\nsites', fontsize=5.5, color='#888888')

# STAT3 to CPS1
ax.annotate('', xy=(1.55, 2.3), xytext=(2.6, 3.35),
            arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=1.2,
                           connectionstyle='arc3,rad=0.25'), zorder=5)

# STAT3 to ASS1
draw_arrow(ax, 3.0, 3.2, 3.0, 1.75, C_ARROW, lw=1.2)

# STAT3 to ASL
ax.annotate('', xy=(3.8, 1.55), xytext=(3.3, 3.2),
            arrowprops=dict(arrowstyle='->', color=C_ARROW, lw=1.0,
                           connectionstyle='arc3,rad=0.2'), zorder=5)

# Urea Cycle cascade: CPS1 → ARG1 → ASS1 → ASL
draw_arrow(ax, 1.85, 2.1, 2.1, 1.6, '#D55E00', lw=1.0)
draw_arrow(ax, 2.6, 1.45, 2.7, 1.45, '#D55E00', lw=1.0)
draw_arrow(ax, 3.35, 1.45, 3.45, 1.45, '#D55E00', lw=1.0)

# ---- Urea output / Serum Urea node (between tissues) ----
urea_x, urea_y = 5.6, 2.6
ellipse = Ellipse((urea_x, urea_y), 1.0, 0.6, facecolor=C_UREA_SIGNAL, edgecolor='white',
                  linewidth=2, alpha=0.9, zorder=8)
ax.add_patch(ellipse)
ax.text(urea_x, urea_y, 'Serum\nUrea ↑', ha='center', va='center', fontsize=7,
        fontweight='bold', color='white', zorder=9)

# Liver → Urea
all_liver_y = [nodes_liver[n][1] for n in ['ARG1', 'ASS1', 'ASL']]
avg_liver_y = np.mean(all_liver_y)
ax.annotate('', xy=(5.0, urea_y), xytext=(4.15, avg_liver_y),
            arrowprops=dict(arrowstyle='->', color=C_UREA_SIGNAL, lw=2.5,
                           connectionstyle='arc3,rad=-0.15'), zorder=6)
draw_label(ax, 4.65, 1.95, 'Nitrogen shunt\nr_Urea=0.52~0.82', fontsize=5.5, color=C_UREA_SIGNAL)

# Urea → Muscle FoxO
ax.annotate('', xy=(7.1, 2.9), xytext=(6.15, urea_y+0.1),
            arrowprops=dict(arrowstyle='->', color=C_UREA_SIGNAL, lw=2.5,
                           connectionstyle='arc3,rad=-0.1'), zorder=6)
draw_label(ax, 6.6, 3.15, 'Systemic\ncatabolic signal?', fontsize=5.5, color=C_UREA_SIGNAL)

# ---- Arrows: Muscle signaling ----
# FoxO → FBXO32/TRIM63
draw_arrow(ax, 7.8, 2.9, 7.8, 2.2, '#0072B2', lw=2.0)
draw_label(ax, 7.45, 2.55, 'transcriptional\nactivation', fontsize=5.5, color='#0072B2')

# FOXO3 → E3 ligases
ax.annotate('', xy=(9.0, 2.2), xytext=(9.0, 2.9),
            arrowprops=dict(arrowstyle='->', color='#0072B2', lw=1.5,
                           connectionstyle='arc3,rad=0.15'), zorder=5)

# IGF1 → RPS6KB1 (mTOR pathway - parallel)
draw_arrow(ax, 8.75, 4.2, 9.45, 4.2, '#009E73', lw=1.8)
draw_label(ax, 9.1, 4.4, 'mTORC1', fontsize=5.5, color='#009E73')

# ---- Protein Deposition endpoint ----
pd_box = FancyBboxPatch((7.6, 0.75), 2.5, 0.55, boxstyle="round,pad=0.1",
                         facecolor='#F5F5F5', edgecolor=C_HIGHLIGHT, linewidth=1.5, zorder=10)
ax.add_patch(pd_box)
ax.text(8.85, 1.03, 'Protein Deposition ↓', ha='center', va='center', fontsize=8,
        fontweight='bold', color=C_HIGHLIGHT, zorder=11)

# FBXO32/TRIM63 → Protein Deposition
draw_arrow(ax, 8.0, 1.55, 8.85, 1.25, C_HIGHLIGHT, lw=2.2)
draw_label(ax, 8.1, 1.65, 'Ubiquitin-\nproteasome', fontsize=5.5, color=C_HIGHLIGHT)
# TRIM63 → PD
ax.annotate('', xy=(8.7, 1.2), xytext=(9.3, 1.55),
            arrowprops=dict(arrowstyle='->', color=C_HIGHLIGHT, lw=1.5,
                           connectionstyle='arc3,rad=-0.1'), zorder=5)

# IGF1/RPS6KB1 → PD (anabolic arm, dashed)
ax.annotate('', xy=(8.5, 1.2), xytext=(10.0, 3.95),
            arrowprops=dict(arrowstyle='->', color='#009E73', lw=1.3, ls='--',
                           connectionstyle='arc3,rad=-0.5'), zorder=5)
draw_label(ax, 10.5, 2.5, 'Anabolic\narm', fontsize=5.5, color='#009E73')

# ---- Correlation callouts ----
# r values at key edges
corr_annotations = [
    (4.5, 3.8, 'r_Urea=0.52', C_GOLD),
    (4.35, 2.1, 'r_Urea=0.82', '#D55E00'),
    (7.0, 3.8, 'r_PD = -0.92', C_HIGHLIGHT),
    (7.1, 0.95, 'r_PD = -0.92', C_HIGHLIGHT),
]
for x, y, txt, color in corr_annotations:
    ax.annotate(txt, xy=(x, y), fontsize=6.5, fontweight='bold', color=color,
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                         edgecolor=color, linewidth=0.6, alpha=0.85))

# ---- Pathway labels (top of compartments) ----
pathway_labels = {
    (1.4, 4.3):  ('JAK/STAT', C_GOLD),
    (3.0, 4.3):  ('Transcriptional\nHub', '#CC79A7'),
    (1.4, 2.8):  ('Urea Cycle\n(NH3 -> Urea)', '#D55E00'),
    (7.8, 3.9):  ('FoxO/\nProteolysis', '#0072B2'),
    (8.4, 4.9):  ('mTOR/\nAnabolic', '#009E73'),
    (7.8, 2.5):  ('UPS/E3 Ligase', C_HIGHLIGHT),
}
for (x, y), (txt, color) in pathway_labels.items():
    ax.annotate(txt, xy=(x, y), fontsize=6, color=color, fontweight='bold',
                ha='center', va='center',
                arrowprops=dict(arrowstyle='->', color=color, lw=0.5, ls='--'))

# ---- Bottom legend / scoring formula ----
formula_box = FancyBboxPatch((0.5, 0.08), 11.0, 0.45, boxstyle="round,pad=0.08",
                              facecolor='white', edgecolor='#CCCCCC', linewidth=0.8, zorder=15)
ax.add_patch(formula_box)
formula_text = (
    'Closure Score = w1 x Pathway Coherence + w2 x PPI Connectivity + '
    'w3 x Literature Support + w4 x Experimental Operability    |    '
    'Ranked by biological logic, not r-value magnitude'
)
ax.text(6.0, 0.3, formula_text, ha='center', va='center', fontsize=6.5,
        color='#555555', zorder=16)

# ---- Save ----
outpath = '/Users/hezongze/fig_cross_tissue_axis.png'
fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {outpath}')
