#!/usr/bin/env python3
"""
Figure: Experimental design overview for ALDH18A1/Glu metabolic branch point validation.
Three-phase experiment with decision gates.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.set_xlim(0, 14)
ax.set_ylim(0, 9)
ax.axis('off')

# Color scheme
C_PHASE1 = '#0072B2'
C_PHASE2 = '#009E73'
C_PHASE3 = '#D55E00'
C_GATE   = '#D41159'
C_BG1    = '#EBF5FB'
C_BG2    = '#E8F8F5'
C_BG3    = '#FEF5E7'

def draw_box(ax, x, y, w, h, color, bg, text, fs=7, fw='bold', z=5):
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                          facecolor=bg, edgecolor=color, linewidth=1.5, zorder=z)
    ax.add_patch(rect)

def arrow(ax, x1, y1, x2, y2, color='#555555', lw=1.5, z=3):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw), zorder=z)

def gate(ax, x, y, text, color=C_GATE):
    diamond = plt.Polygon([[x, y+0.35], [x+0.5, y], [x, y-0.35], [x-0.5, y]],
                          facecolor='white', edgecolor=color, linewidth=2, zorder=8)
    ax.add_patch(diamond)
    ax.text(x, y, text, ha='center', va='center', fontsize=5.5, fontweight='bold', color=color, zorder=9)

# ==============================================================
# Phase 1: Molecular verification
# ==============================================================
ax.text(0.3, 8.7, 'PHASE 1', fontsize=12, fontweight='bold', color=C_PHASE1)
ax.text(2.3, 8.7, 'Molecular Branch Point Confirmation (1-2 months)', fontsize=9, color=C_PHASE1)

# Exp 1.1
draw_box(ax, 0.3, 7.2, 3.8, 1.2, C_PHASE1, C_BG1, '')
ax.text(2.2, 8.1, 'Exp 1.1: Serum Pro Retrospective', fontsize=7.5, fontweight='bold', color=C_PHASE1)
ax.text(2.2, 7.8, 'UPLC-MS/MS targeted AA panel\n48 serum samples (DLY/TFB x 4 stages)\nPro, Hyp, Glu, Orn, Cit, Arg, BCAA', fontsize=6, ha='center', va='top')

# Exp 1.2
draw_box(ax, 4.5, 7.2, 3.8, 1.2, C_PHASE1, C_BG1, '')
ax.text(6.4, 8.1, 'Exp 1.2: Liver ALDH18A1 Validation', fontsize=7.5, fontweight='bold', color=C_PHASE1)
ax.text(6.4, 7.8, 'WB: ALDH18A1, PYCR1, OAT, CPS1\nIHC: zonal distribution\nLiver AA metabolome (Glu, Pro, P5C, Orn, Cit)', fontsize=6, ha='center', va='top')

# Exp 1.3
draw_box(ax, 8.7, 7.2, 3.8, 1.2, C_PHASE1, C_BG1, '')
ax.text(10.6, 8.1, 'Exp 1.3: Muscle mTORC1 Status', fontsize=7.5, fontweight='bold', color=C_PHASE1)
ax.text(10.6, 7.8, 'WB: p-S6K1, p-4E-BP1, p-AKT\nEPRS protein (prolyl-tRNA synthetase)\nMuscle fiber CSA', fontsize=6, ha='center', va='top')

# Gate 1
arrow(ax, 6.4, 7.1, 6.4, 6.4, C_PHASE1, 1.5)
gate(ax, 6.4, 6.1, 'Pro > TFB?\nALDH18A1\nprotein OK?', C_GATE)
ax.text(6.4, 5.65, 'YES → Phase 2', fontsize=6, color=C_PHASE2, ha='center', fontweight='bold')
ax.text(3.5, 6.2, 'NO → Revise\nhypothesis', fontsize=5, color='#999999', ha='center')

# ==============================================================
# Phase 2: In vitro mechanism
# ==============================================================
ax.text(0.3, 5.2, 'PHASE 2', fontsize=12, fontweight='bold', color=C_PHASE2)
ax.text(2.3, 5.2, 'Hepatocyte Mechanism & Muscle Crosstalk (3-4 months)', fontsize=9, color=C_PHASE2)

# Exp 2.1
draw_box(ax, 0.3, 3.7, 3.8, 1.3, C_PHASE2, C_BG2, '')
ax.text(2.2, 4.7, 'Exp 2.1: ALDH18A1 KD in Hepatocytes', fontsize=7.5, fontweight='bold', color=C_PHASE2)
ax.text(2.2, 4.35, 'Porcine primary hepatocytes\nsiRNA-ALDH18A1 x 3 + NC\nMeasure: Glu/Pro/Orn/Urea in media\n^13C-Glu tracing: ^13C-Pro vs ^13C-Urea\nRescue: +Pro / +B6 (PLP cofactor)', fontsize=6, ha='center', va='top')

# Exp 2.2
draw_box(ax, 4.5, 3.7, 3.8, 1.3, C_PHASE2, C_BG2, '')
ax.text(6.4, 4.7, 'Exp 2.2: Conditioned Medium Crosstalk', fontsize=7.5, fontweight='bold', color=C_PHASE2)
ax.text(6.4, 4.35, 'KD-hepatocyte CM on C2C12 myotubes\nSUnSET (protein synthesis)\nProteasome activity (degradation)\np-S6K1, p-4E-BP1, FOXO1 localization\nRescue: KD-CM + Pro (0.2 mM)', fontsize=6, ha='center', va='top')

# Exp 2.3: OE (optional)
draw_box(ax, 8.7, 3.7, 3.8, 1.3, C_PHASE2, C_BG2, '')
ax.text(10.6, 4.7, 'Exp 2.3: ALDH18A1 OE (Optional)', fontsize=7.5, fontweight='bold', color=C_PHASE2)
ax.text(10.6, 4.35, 'ALDH18A1 overexpression plasmid\nTFB-genotype hepatocyte line\nVerify: Pro secretion increase\nUrea reduction\nReversal of TFB metabolic pattern', fontsize=6, ha='center', va='top')

# Gate 2
arrow(ax, 6.4, 3.6, 6.4, 2.9, C_PHASE2, 1.5)
gate(ax, 6.4, 2.6, 'Pro/Urea\nshift + rescue\nconfirmed?', C_GATE)
ax.text(6.4, 2.15, 'YES → Phase 3', fontsize=6, color=C_PHASE3, ha='center', fontweight='bold')

# ==============================================================
# Phase 3: In vivo validation
# ==============================================================
ax.text(0.3, 1.6, 'PHASE 3', fontsize=12, fontweight='bold', color=C_PHASE3)
ax.text(2.3, 1.6, 'In Vivo Causal Validation (4-6 months)', fontsize=9, color=C_PHASE3)

# Exp 3.1
draw_box(ax, 0.3, 0.3, 6.0, 1.1, C_PHASE3, C_BG3, '')
ax.text(3.3, 1.15, 'Exp 3.1: Dietary Pro+Glu Supplementation', fontsize=8, fontweight='bold', color=C_PHASE3)
ax.text(3.3, 0.85, '2x2 factorial: DLY/TFB x Basal/+Pro(0.5%)+Glu(1%), n=10/group, 25-75kg, ~60d', fontsize=6.5, ha='center')
ax.text(3.3, 0.6, 'N balance (total collection) + Serum AA timecourse + Carcass traits + mTORC1 WB', fontsize=6, ha='center')

# Exp 3.2 (optional)
draw_box(ax, 6.7, 0.3, 5.8, 1.1, C_PHASE3, C_BG3, '')
ax.text(9.6, 1.15, 'Exp 3.2: AAV8-ALDH18A1 (Optional)', fontsize=8, fontweight='bold', color=C_PHASE3)
ax.text(9.6, 0.85, 'Liver-specific ALDH18A1 OE via AAV8-TBG in TFB pigs (n=5)\nDirect causal test: does ALDH18A1 OE rescue TFB protein deposition?', fontsize=6.5, ha='center')

# ==============================================================
# Timeline bar at bottom
# ==============================================================
timeline_y = 0.05
ax.text(0.3, timeline_y, 'Month: 1    2    3    4    5    6    7    8    9    10',
        fontsize=6, color='#888888', fontfamily='monospace')

# Phase bars
phases_timeline = [
    (0.3, timeline_y-0.03, 1.8, C_PHASE1, 'Phase 1'),
    (2.2, timeline_y-0.03, 2.5, C_PHASE2, 'Phase 2'),
    (4.8, timeline_y-0.03, 4.5, C_PHASE3, 'Phase 3'),
]
for x, y, w, c, label in phases_timeline:
    rect = FancyBboxPatch((x, y), w, 0.04, boxstyle="round,pad=0.01",
                          facecolor=c, edgecolor='none', alpha=0.6)
    ax.add_patch(rect)

# ==============================================================
# Central hypothesis callout
# ==============================================================
hypo_box = FancyBboxPatch((0.3, 4.9), 12.2, 0.55, boxstyle="round,pad=0.1",
                           facecolor='white', edgecolor='#999999', linewidth=1, alpha=0.9, zorder=1)
ax.add_patch(hypo_box)

# Mini pathway diagram inside
ax.text(6.4, 5.35, 'Glu', fontsize=7, fontweight='bold', color='#333333', ha='center')
ax.text(6.4, 5.15, 'ALDH18A1 (P5CS)', fontsize=7, fontweight='bold', color=C_GATE, ha='center')
ax.text(5.0, 4.97, 'Pro → mTORC1 → Protein Deposition', fontsize=6, color=C_PHASE2, ha='left')
ax.text(8.0, 4.97, 'Orn → Cit → Urea → N loss', fontsize=6, color=C_PHASE3, ha='left')
ax.annotate('', xy=(5.5, 5.15), xytext=(5.9, 5.15),
            arrowprops=dict(arrowstyle='->', color=C_PHASE2, lw=1), zorder=2)
ax.annotate('', xy=(7.3, 5.15), xytext=(6.9, 5.15),
            arrowprops=dict(arrowstyle='->', color=C_PHASE3, lw=1), zorder=2)

fig.savefig('/Users/hezongze/ALDH18A1_project/03_main_figures_ALDH18A1/fig_experiment_design.png', dpi=300, facecolor='white')
plt.close(fig)
print('Saved: /Users/hezongze/ALDH18A1_project/03_main_figures_ALDH18A1/fig_experiment_design.png')
