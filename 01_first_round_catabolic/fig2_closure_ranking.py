#!/usr/bin/env python3
"""
Figure 2: Biological Closure Score Ranking.
Horizontal bar chart showing top 5 liver and top 5 muscle genes ranked by
composite closure score. Color-coded by functional pathway.
Includes inset annotations for key features (PPI, pathway, literature support).

Reference: JASB 2026 multi-panel figures.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7.5,
    'legend.fontsize': 6.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# ---- Data ----
# Liver genes (closure score, pathway, PPI, r_Urea, key feature)
liver_data = [
    ('ARG1',    16.3, 'Urea Cycle',       58,  0.65, '4 STAT3 promoter sites'),
    ('IL6',     16.2, 'JAK/STAT Upstream',137, None, 'STAT3 upstream activator'),
    ('ASS1',    16.2, 'Urea Cycle',       51,  0.48, 'Rate-limiting step'),
    ('STAT3',   16.1, 'Master Regulator', 124, 0.52, 'Core transcriptional hub'),
    ('ASL',     16.0, 'Urea Cycle',       44,  0.82, 'Highest r_Urea in set'),
]

muscle_data = [
    ('FOXO1',   16.8, 'FoxO/Proteolysis', 63,  -0.78, 'Master UPS regulator'),
    ('FOXO3',   16.5, 'FoxO/Proteolysis', 46,  -0.71, 'Redundant with FOXO1'),
    ('IGF1',    15.7, 'mTOR/Anabolic',    81,   0.68, 'Growth axis core'),
    ('RPS6KB1', 15.4, 'mTOR/Translation', 26,   0.55, 'S6K1, translation effector'),
    ('FBXO32',  15.0, 'UPS/E3 Ligase',    18,  -0.92, 'Atrogin-1, r_PD = −0.92'),
]

# Color mapping by pathway
pathway_colors = {
    'Urea Cycle':        '#D55E00',
    'JAK/STAT Upstream': '#E69F00',
    'Master Regulator':  '#CC79A7',
    'FoxO/Proteolysis':  '#0072B2',
    'mTOR/Anabolic':     '#009E73',
    'mTOR/Translation':  '#56B4E9',
    'UPS/E3 Ligase':     '#D41159',
}

# ---- Layout ----
fig = plt.figure(figsize=(9, 5.5))

# Panel layout
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.08,
                       left=0.12, right=0.94, top=0.90, bottom=0.12)

# ---- Liver Panel ----
ax_l = fig.add_subplot(gs[0])
genes_l = [d[0] for d in liver_data][::-1]
scores_l = [d[1] for d in liver_data][::-1]
pathways_l = [d[2] for d in liver_data][::-1]
colors_l = [pathway_colors[p] for p in pathways_l]
features_l = [d[5] for d in liver_data][::-1]
ppis_l = [d[3] for d in liver_data][::-1]

bars_l = ax_l.barh(range(len(genes_l)), scores_l, height=0.6, color=colors_l,
                    edgecolor='white', linewidth=0.8, alpha=0.85)
ax_l.set_yticks(range(len(genes_l)))
ax_l.set_yticklabels(genes_l, fontsize=8, fontweight='bold')
ax_l.set_xlim(14.5, 17.2)
ax_l.set_xlabel('Biological Closure Score', fontsize=8)
ax_l.set_title('Liver', fontsize=10, fontweight='bold', color='#D55E00', pad=8)

# Score labels
for i, (bar, score, ppi, feat) in enumerate(zip(bars_l, scores_l, ppis_l, features_l)):
    ax_l.text(score + 0.05, bar.get_y() + bar.get_height()/2,
              f'{score:.1f}', va='center', fontsize=7, fontweight='bold', color='#333333')

# Spine
ax_l.spines['top'].set_visible(False)
ax_l.spines['right'].set_visible(False)
ax_l.tick_params(labelsize=7)
ax_l.grid(axis='x', alpha=0.3, lw=0.4)

# Panel label
ax_l.text(-0.25, 1.02, 'A', transform=ax_l.transAxes, fontsize=12, fontweight='bold')

# ---- Muscle Panel ----
ax_m = fig.add_subplot(gs[1])
genes_m = [d[0] for d in muscle_data][::-1]
scores_m = [d[1] for d in muscle_data][::-1]
pathways_m = [d[2] for d in muscle_data][::-1]
colors_m = [pathway_colors[p] for p in pathways_m]
features_m = [d[5] for d in muscle_data][::-1]
r_vals_m  = [d[4] for d in muscle_data][::-1]

bars_m = ax_m.barh(range(len(genes_m)), scores_m, height=0.6, color=colors_m,
                    edgecolor='white', linewidth=0.8, alpha=0.85)
ax_m.set_yticks(range(len(genes_m)))
ax_m.set_yticklabels(genes_m, fontsize=8, fontweight='bold')
ax_m.set_xlim(14.5, 17.2)
ax_m.set_xlabel('Biological Closure Score', fontsize=8)
ax_m.set_title('Muscle', fontsize=10, fontweight='bold', color='#0072B2', pad=8)

for i, (bar, score, feat) in enumerate(zip(bars_m, scores_m, features_m)):
    ax_m.text(score + 0.05, bar.get_y() + bar.get_height()/2,
              f'{score:.1f}', va='center', fontsize=7, fontweight='bold', color='#333333')

ax_m.spines['top'].set_visible(False)
ax_m.spines['right'].set_visible(False)
ax_m.tick_params(labelsize=7)
ax_m.grid(axis='x', alpha=0.3, lw=0.4)
ax_m.text(-0.05, 1.02, 'B', transform=ax_m.transAxes, fontsize=12, fontweight='bold')

# ---- Shared legend (bottom) ----
legend_elements = [
    mpatches.Patch(facecolor=pathway_colors['Urea Cycle'],        alpha=0.85, label='Urea Cycle'),
    mpatches.Patch(facecolor=pathway_colors['JAK/STAT Upstream'], alpha=0.85, label='JAK/STAT Upstream'),
    mpatches.Patch(facecolor=pathway_colors['Master Regulator'],  alpha=0.85, label='Master Regulator'),
    mpatches.Patch(facecolor=pathway_colors['FoxO/Proteolysis'],  alpha=0.85, label='FoxO/Proteolysis'),
    mpatches.Patch(facecolor=pathway_colors['mTOR/Anabolic'],     alpha=0.85, label='mTOR/Anabolic'),
    mpatches.Patch(facecolor=pathway_colors['mTOR/Translation'],  alpha=0.85, label='mTOR/Translation'),
    mpatches.Patch(facecolor=pathway_colors['UPS/E3 Ligase'],     alpha=0.85, label='UPS/E3 Ligase'),
]
leg = fig.legend(handles=legend_elements, loc='lower center', ncol=7,
                 frameon=True, fancybox=True, framealpha=0.9, fontsize=6.5,
                 bbox_to_anchor=(0.5, -0.02))
leg.get_frame().set_linewidth(0.3)

# ---- Annotation callout boxes ----
# Add key feature callouts on the right side of each bar (liver)
annotations_liver = [
    (16.3, 4, 'PPI=58, 4 STAT3 bs', pathway_colors['Urea Cycle']),
    (16.2, 3, 'PPI=137, cytokine hub', pathway_colors['JAK/STAT Upstream']),
    (16.2, 2, 'PPI=51, urea cycle core', pathway_colors['Urea Cycle']),
    (16.1, 1, 'PPI=124, r_Urea=0.52', pathway_colors['Master Regulator']),
    (16.0, 0, 'r_Urea=0.82, strongest r', pathway_colors['Urea Cycle']),
]
annotations_muscle = [
    (16.8, 4, 'PPI=63, FoxO master', pathway_colors['FoxO/Proteolysis']),
    (16.5, 3, 'PPI=46, FoxO family', pathway_colors['FoxO/Proteolysis']),
    (15.7, 2, 'PPI=81, growth axis', pathway_colors['mTOR/Anabolic']),
    (15.4, 1, 'PPI=26, mTOR effector', pathway_colors['mTOR/Translation']),
    (15.0, 0, 'r_PD=−0.92, strongest', pathway_colors['UPS/E3 Ligase']),
]
for score, y, text, color in annotations_liver:
    ax_l.annotate(text, xy=(score, y), xytext=(16.75, y),
                  fontsize=5.5, color=color, va='center',
                  arrowprops=dict(arrowstyle='->', color=color, lw=0.6))
for score, y, text, color in annotations_muscle:
    ax_m.annotate(text, xy=(score, y), xytext=(16.75, y),
                  fontsize=5.5, color=color, va='center',
                  arrowprops=dict(arrowstyle='->', color=color, lw=0.6))

# ---- Save ----
outpath = '/Users/hezongze/fig_biological_closure_ranking.png'
fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {outpath}')
