#!/usr/bin/env python3
"""
Fig 6: Working model — ALDH18A1/Glu metabolic branch point determines nitrogen partitioning
between anabolic (Pro -> mTORC1 -> protein deposition) and catabolic (Urea -> N excretion) fates.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc, Rectangle, Polygon, Ellipse
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 6.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

fig, ax = plt.subplots(1, 1, figsize=(11, 8.5))
ax.set_xlim(0, 14)
ax.set_ylim(0, 10.5)
ax.axis('off')

# ---------------------------------------------------------------------------
# Color scheme
# ---------------------------------------------------------------------------
C_LIVER    = '#D55E00'
C_BLOOD    = '#CC79A7'
C_MUSCLE   = '#0072B2'
C_ANABOLIC = '#009E73'
C_CATABOLIC = '#E69F00'
C_ENZYME   = '#D41159'
C_DATA     = '#56B4E9'
C_ARROW    = '#555555'
C_BG_LIVER  = '#FFF8F0'
C_BG_BLOOD  = '#FFF0F7'
C_BG_MUSCLE = '#F0F5FF'

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def draw_box(ax, x, y, w, h, text, color, bg, fs=6.5, fw='bold', za=5, alpha=1.0):
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                          facecolor=bg, edgecolor=color, linewidth=1.2, zorder=za, alpha=alpha)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs,
            fontweight=fw, color='#333333', zorder=za+1)

def draw_arrow(ax, x1, y1, x2, y2, color=C_ARROW, lw=1.3, za=3, style='->'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw), zorder=za)

def draw_metab_icon(ax, x, y, text, color, size=0.28):
    circle = Ellipse((x, y), size*2, size, facecolor=color, edgecolor='white',
                     linewidth=0.8, alpha=0.85, zorder=10)
    ax.add_patch(circle)
    ax.text(x, y, text, ha='center', va='center', fontsize=5.5, fontweight='bold',
            color='white', zorder=11)

def data_callout(ax, x, y, text, color=C_DATA):
    ax.annotate(text, xy=(x, y), fontsize=5.5, color=color, fontweight='bold',
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                          edgecolor=color, alpha=0.9, lw=0.7))

# ---------------------------------------------------------------------------
# Compartment backgrounds
# ---------------------------------------------------------------------------
# Liver compartment
rect_liver = FancyBboxPatch((0.3, 4.8), 5.8, 5.2, boxstyle="round,pad=0.3",
                            facecolor=C_BG_LIVER, edgecolor=C_LIVER, linewidth=2, alpha=0.5, zorder=1)
ax.add_patch(rect_liver)
ax.text(3.2, 9.7, 'LIVER', fontsize=10, fontweight='bold', color=C_LIVER, ha='center', zorder=2)

# Blood compartment
rect_blood = FancyBboxPatch((0.3, 2.4), 13.4, 2.2, boxstyle="round,pad=0.2",
                            facecolor=C_BG_BLOOD, edgecolor=C_BLOOD, linewidth=1.5, alpha=0.4, zorder=1)
ax.add_patch(rect_blood)
ax.text(7.0, 4.3, 'CIRCULATION', fontsize=9, fontweight='bold', color=C_BLOOD, ha='center', zorder=2)

# Muscle compartment
rect_muscle = FancyBboxPatch((6.4, 0.3), 7.3, 5.2, boxstyle="round,pad=0.3",
                             facecolor=C_BG_MUSCLE, edgecolor=C_MUSCLE, linewidth=2, alpha=0.5, zorder=1)
ax.add_patch(rect_muscle)
ax.text(10.05, 5.2, 'SKELETAL MUSCLE', fontsize=10, fontweight='bold', color=C_MUSCLE, ha='center', zorder=2)

# Section labels
ax.text(0.5, 10.2, 'a', fontsize=13, fontweight='bold', color=C_LIVER, ha='center')
ax.text(0.5, 4.5, 'b', fontsize=13, fontweight='bold', color=C_BLOOD, ha='center')
ax.text(6.6, 5.4, 'c', fontsize=13, fontweight='bold', color=C_MUSCLE, ha='center')

# ===========================================================================
# LIVER PANEL (a)
# ===========================================================================

# TCA cycle / alpha-KG source
draw_box(ax, 1.0, 8.3, 2.0, 0.7, 'TCA Cycle\n(alpha-KG)', C_LIVER, 'white', 6, 'normal', 4)

# Glu node
draw_metab_icon(ax, 2.0, 7.5, 'Glu', C_LIVER)

# Arrow alpha-KG -> Glu
draw_arrow(ax, 2.0, 8.3, 2.0, 7.8, C_LIVER, 1.5, 3)

# ALDH18A1 enzyme box (the BRANCH POINT - highlighted)
draw_box(ax, 0.7, 6.3, 2.6, 0.8, 'ALDH18A1\n(P5CS)', C_ENZYME, '#FFF0F0', 7, 'bold', 8)

# Arrow Glu -> ALDH18A1
draw_arrow(ax, 2.0, 7.2, 2.0, 6.9, C_ENZYME, 2.0, 7)

# Branch arrows from ALDH18A1
# Left branch: Pro (anabolic)
draw_arrow(ax, 1.3, 6.3, 0.8, 5.4, C_ANABOLIC, 1.8, 5)
# Right branch: Orn -> Urea (catabolic)
draw_arrow(ax, 2.7, 6.3, 4.2, 6.3, C_CATABOLIC, 1.8, 5)

# Pro synthesis
draw_box(ax, 0.2, 4.9, 1.2, 0.6, 'P5C', C_ANABOLIC, '#F0FFF5', 5.5, 'normal', 5)
draw_metab_icon(ax, 0.8, 4.4, 'Pro', C_ANABOLIC)

# Orn -> Urea axis
draw_box(ax, 3.5, 5.9, 1.0, 0.55, 'Orn', C_CATABOLIC, '#FFF8E1', 5.5, 'normal', 5)
draw_box(ax, 4.7, 5.9, 1.0, 0.55, 'Cit', C_CATABOLIC, '#FFF8E1', 5.5, 'normal', 5)
draw_box(ax, 3.5, 5.2, 2.2, 0.5, 'Urea Cycle', C_CATABOLIC, '#FFF0E5', 5.5, 'bold', 5)
draw_metab_icon(ax, 4.6, 4.6, 'Urea', C_CATABOLIC)
draw_arrow(ax, 3.5, 6.17, 4.7, 6.17, C_CATABOLIC, 1, 5)
draw_arrow(ax, 4.6, 5.2, 4.6, 5.0, C_CATABOLIC, 1.2, 5)

# Pro export to blood
draw_arrow(ax, 0.8, 4.2, 0.8, 3.8, C_ANABOLIC, 1.5, 6)
# Urea export to blood
draw_arrow(ax, 4.6, 4.4, 4.6, 3.8, C_CATABOLIC, 1.2, 6)

# PYCR1
draw_box(ax, 0.2, 3.5, 1.2, 0.6, 'PYCR1', C_LIVER, 'white', 5.5, 'normal', 5)
draw_arrow(ax, 0.8, 4.9, 0.8, 4.1, C_LIVER, 0.8, 4)

# Data callouts in liver
data_callout(ax, 3.5, 8.8, 'DLY: ALDH18A1 mRNA high\nFC=+0.85 (75kg)')
data_callout(ax, 5.5, 7.2, 'DLY: CPS1 low\nTFB: CPS1 high')

# ===========================================================================
# BLOOD PANEL (b)
# ===========================================================================

# Serum metabolites
metab_positions = {
    'Pro':  (1.5, 3.3, C_ANABOLIC),
    'Glu':  (4.0, 3.3, C_LIVER),
    'Orn':  (6.5, 3.3, C_CATABOLIC),
    'Cit':  (9.0, 3.3, C_CATABOLIC),
    'Urea': (11.5, 3.3, C_CATABOLIC),
}

for name, (mx, my, mc) in metab_positions.items():
    draw_metab_icon(ax, mx, my, name, mc, 0.30)

# Data callouts
data_callout(ax, 1.5, 2.9, 'DLY > TFB\n(Phase 1 test)', C_DATA)
data_callout(ax, 4.0, 2.9, 'DLY > TFB\np=0.005 (45kg)', C_DATA)
data_callout(ax, 6.5, 2.9, 'TFB > DLY\np=0.009 (75kg)', C_DATA)
data_callout(ax, 9.0, 2.9, 'TFB > DLY\np<0.001 (45kg)', C_DATA)
data_callout(ax, 11.5, 2.9, 'TFB > DLY\np<0.001', C_DATA)

# Pro enters muscle
draw_arrow(ax, 1.8, 3.3, 7.0, 2.5, C_ANABOLIC, 1.8, 6)
ax.text(4.0, 3.05, 'Pro transport', fontsize=5.5, color=C_ANABOLIC, fontstyle='italic', ha='center')

# Urea excretion arrow
draw_arrow(ax, 11.5, 3.0, 11.5, 2.2, C_CATABOLIC, 1.2, 6)
ax.text(11.5, 2.65, 'N excretion', fontsize=5, color=C_CATABOLIC, fontstyle='italic', ha='center')

# ===========================================================================
# MUSCLE PANEL (c)
# ===========================================================================

# Pro enters muscle cell
draw_box(ax, 7.2, 4.3, 1.8, 0.6, 'Pro uptake\n(SLC36A1/PAT1)', C_MUSCLE, 'white', 5.5, 'normal', 6)

# EPRS sensor (KEY MECHANISM)
draw_box(ax, 8.0, 3.3, 2.5, 0.7, 'EPRS senses\nProlyl-tRNA charge', C_ENZYME, '#FFF0F0', 6, 'bold', 8)
draw_arrow(ax, 8.1, 4.3, 9.25, 4.0, C_MUSCLE, 1.2, 5)
data_callout(ax, 10.8, 3.65, 'Non-classical\nmTORC1 activation\n(AKT-independent)', C_ENZYME)

# mTORC1 translocation
draw_box(ax, 7.8, 2.3, 3.0, 0.7, 'mTORC1 lysosomal\ntranslocation (Rheb)', C_MUSCLE, 'white', 6, 'bold', 7)
draw_arrow(ax, 9.25, 3.3, 9.3, 3.0, C_MUSCLE, 1.5, 6)

# S6K1 and 4E-BP1
draw_box(ax, 7.0, 1.3, 2.0, 0.65, 'S6K1-P\n(rpS6 -> translation)', C_ANABOLIC, '#F0FFF5', 5.5, 'bold', 7)
draw_box(ax, 9.5, 1.3, 2.0, 0.65, '4E-BP1-P\n(eIF4E release)', C_ANABOLIC, '#F0FFF5', 5.5, 'bold', 7)

draw_arrow(ax, 8.5, 2.3, 8.0, 1.95, C_ANABOLIC, 1.2, 6)
draw_arrow(ax, 10.0, 2.3, 10.5, 1.95, C_ANABOLIC, 1.2, 6)

# Final output
draw_box(ax, 7.8, 0.4, 3.0, 0.65, 'Global protein translation\n-> MUSCLE PROTEIN DEPOSITION', C_ANABOLIC, '#E8F8F5',
         7, 'bold', 8)
draw_arrow(ax, 8.0, 1.3, 8.5, 1.05, C_ANABOLIC, 1.2, 6)
draw_arrow(ax, 10.5, 1.3, 10.0, 1.05, C_ANABOLIC, 1.2, 6)

# Data callouts in muscle
data_callout(ax, 7.0, 2.8, 'DLY: p-4E-BP1 high\n(data confirmed)', C_DATA)
data_callout(ax, 12.5, 1.9, 'p-AKT: NS\n(DLY = TFB)\n-> non-classical pathway', C_ENZYME)

# ===========================================================================
# Bottom summary box
# ===========================================================================
summary_box = FancyBboxPatch((0.5, 0.05), 13.0, 0.25, boxstyle="round,pad=0.05",
                             facecolor='#333333', edgecolor='none', alpha=0.08, zorder=1)
ax.add_patch(summary_box)
ax.text(7.0, 0.17, 'ALDH18A1 controls Glu partitioning: high ALDH18A1 -> Pro -> mTORC1 -> protein synthesis  |  low ALDH18A1 -> Orn -> Urea -> N loss',
        ha='center', fontsize=7, fontweight='bold', color='#333333', fontstyle='italic')

# ===========================================================================
# Legend
# ===========================================================================
legend_items = [
    (0.5,  10.5, C_ENZYME, 'Key enzyme / regulatory node'),
    (3.0,  10.5, C_ANABOLIC, 'Anabolic axis (N retention)'),
    (5.5,  10.5, C_CATABOLIC, 'Catabolic axis (N loss)'),
    (8.0,  10.5, C_DATA, 'Data-supported evidence'),
]
for x, y, c, label in legend_items:
    rect = FancyBboxPatch((x, y-0.06), 0.25, 0.12, boxstyle="round,pad=0.02",
                          facecolor=c, edgecolor='none', alpha=0.8)
    ax.add_patch(rect)
    ax.text(x + 0.35, y, label, fontsize=6, va='center', color='#333333')

# ---------------------------------------------------------------------------
outpath = '/Users/hezongze/ALDH18A1_project/03_main_figures_ALDH18A1/fig6_working_model.png'
# Data source annotation
fig.text(0.5, 0.008, 'Working model integrating: ALDH18A1/P5CS biochemistry (Perez-Arellano et al. 2010), Pro-EPRS-mTORC1 non-classical pathway (Kim et al. 2015; Arif et al. 2017), mTORC1 lysosomal translocation (Sancak et al. 2008). Data callouts indicate nodes supported by measured data in this study.', ha='center', fontsize=5, color='#555555', fontstyle='italic')
fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {outpath}')
