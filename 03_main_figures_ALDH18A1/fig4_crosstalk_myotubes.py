#!/usr/bin/env python3
"""
Fig 4: Conditioned medium from ALDH18A1-KD hepatocytes reduces C2C12 myotube protein synthesis;
Pro rescue restores the phenotype.
Panel A: SUnSET puromycin incorporation (protein synthesis rate).
Panel B: p-S6K1/S6K1 ratio.
Panel C: Myotube diameter (MyHC+ fibers).
Panel D: Representative immunofluorescence images (schematic).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.titlesize': 8,
    'axes.labelsize': 7,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

rng = np.random.default_rng(456)
n = 3  # 3 independent CM preparations

C_NCCM   = '#888888'
C_KDCM   = '#D55E00'
C_KDPRO  = '#009E73'
C_NCUREA = '#CC79A7'

groups  = ['NC-CM', 'KD-CM', 'KD-CM\n+Pro', 'NC-CM\n+Urea']
colors  = [C_NCCM, C_KDCM, C_KDPRO, C_NCUREA]
g_names = ['NC-CM', 'KD-CM', 'KD-CM+Pro', 'NC-CM+Urea']

def sig_str(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'ns'

# ---------------------------------------------------------------------------
# Simulated data
# ---------------------------------------------------------------------------
# Panel A: SUnSET (normalized to NC-CM = 1.0)
sunset = {
    'NC-CM':      np.array([1.00, 0.96, 1.04]),
    'KD-CM':      np.array([0.62, 0.58, 0.68]),
    'KD-CM+Pro':  np.array([0.85, 0.80, 0.90]),
    'NC-CM+Urea': np.array([0.98, 0.94, 1.02]),
}

# Panel B: p-S6K1/S6K1
ps6k1 = {
    'NC-CM':      np.array([1.00, 0.95, 1.05]),
    'KD-CM':      np.array([0.55, 0.50, 0.62]),
    'KD-CM+Pro':  np.array([0.80, 0.75, 0.88]),
    'NC-CM+Urea': np.array([0.96, 0.92, 1.02]),
}

# Panel C: Myotube diameter (um)
diameter = {
    'NC-CM':      np.array([18.5, 19.2, 18.0, 19.8, 18.8, 17.5, 19.0, 18.2, 19.5, 20.0,
                             17.8, 18.6, 19.3, 18.1, 19.0, 17.2, 18.5, 20.2, 18.9, 19.1]),
    'KD-CM':      np.array([15.2, 14.8, 16.0, 15.5, 14.2, 15.8, 14.5, 16.2, 15.0, 14.8,
                             15.5, 15.2, 16.0, 14.5, 15.8, 16.5, 14.2, 15.0, 15.5, 14.0]),
    'KD-CM+Pro':  np.array([17.5, 17.0, 18.2, 16.8, 17.8, 17.2, 18.0, 16.5, 17.5, 18.2,
                             16.8, 17.8, 17.0, 18.5, 17.2, 17.5, 16.8, 18.0, 17.8, 17.0]),
    'NC-CM+Urea': np.array([18.8, 18.2, 19.5, 18.0, 19.2, 17.8, 18.5, 19.0, 18.2, 18.8,
                             19.5, 18.0, 18.5, 19.2, 18.8, 17.5, 19.0, 18.2, 19.2, 18.5]),
}

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(9.5, 7.5))
gs_outer = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[1, 1],
                              left=0.07, right=0.97, top=0.95, bottom=0.07,
                              hspace=0.45, wspace=0.55)

# ---- Panel A: SUnSET ----
ax_a = fig.add_subplot(gs_outer[0, 0])
ax_a.text(-0.15, 1.05, 'A', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_a.transAxes)
ax_a.set_title('Protein synthesis (SUnSET)', fontsize=8, fontweight='bold', pad=4)

for i, g in enumerate(g_names):
    d = sunset[g]
    m, s = np.mean(d), np.std(d, ddof=1) / np.sqrt(n)
    ax_a.bar(i, m, 0.5, color=colors[i], edgecolor='white', linewidth=0.5)
    ax_a.errorbar(i, m, yerr=s, color='#333333', capsize=3, lw=0.8)
    for v in d:
        ax_a.scatter(i + rng.uniform(-0.06, 0.06), v, color=colors[i], s=10, alpha=0.5, zorder=5)

# Key comparisons
y_m = 1.1
ax_a.plot([0, 0, 1, 1], [y_m, y_m+0.04, y_m+0.04, y_m], color='#555555', lw=0.6)
ax_a.text(0.5, y_m+0.05, '***', ha='center', fontsize=8, fontweight='bold', color='#D41159')
ax_a.plot([1, 1, 2, 2], [1.0, 1.04, 1.04, 1.0], color='#555555', lw=0.6)
ax_a.text(1.5, 1.05, '*', ha='center', fontsize=8, fontweight='bold', color='#D41159')

ax_a.set_xticks(range(4))
ax_a.set_xticklabels(groups, fontsize=6)
ax_a.set_ylabel('Puromycin incorporation\n(relative to NC-CM)', fontsize=7)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.grid(axis='y', alpha=0.25, lw=0.4)
ax_a.set_ylim(0, 1.22)

# ---- Panel B: p-S6K1 ----
ax_b = fig.add_subplot(gs_outer[0, 1])
ax_b.text(-0.15, 1.05, 'B', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_b.transAxes)
ax_b.set_title('mTORC1 activity (p-S6K1/S6K1)', fontsize=8, fontweight='bold', pad=4)

for i, g in enumerate(g_names):
    d = ps6k1[g]
    m, s = np.mean(d), np.std(d, ddof=1) / np.sqrt(n)
    ax_b.bar(i, m, 0.5, color=colors[i], edgecolor='white', linewidth=0.5)
    ax_b.errorbar(i, m, yerr=s, color='#333333', capsize=3, lw=0.8)
    for v in d:
        ax_b.scatter(i + rng.uniform(-0.06, 0.06), v, color=colors[i], s=10, alpha=0.5, zorder=5)

y_m2 = 1.12
ax_b.plot([0, 0, 1, 1], [y_m2, y_m2+0.04, y_m2+0.04, y_m2], color='#555555', lw=0.6)
ax_b.text(0.5, y_m2+0.05, '***', ha='center', fontsize=8, fontweight='bold', color='#D41159')
ax_b.plot([1, 1, 2, 2], [1.02, 1.06, 1.06, 1.02], color='#555555', lw=0.6)
ax_b.text(1.5, 1.07, '*', ha='center', fontsize=8, fontweight='bold', color='#D41159')

ax_b.set_xticks(range(4))
ax_b.set_xticklabels(groups, fontsize=6)
ax_b.set_ylabel('p-S6K1 / S6K1\n(relative to NC-CM)', fontsize=7)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.grid(axis='y', alpha=0.25, lw=0.4)
ax_b.set_ylim(0, 1.22)

# ---- Panel C: Myotube diameter ----
ax_c = fig.add_subplot(gs_outer[0, 2])
ax_c.text(-0.15, 1.05, 'C', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_c.transAxes)
ax_c.set_title('Myotube diameter', fontsize=8, fontweight='bold', pad=4)

positions = [0, 1, 2, 3]
bp_data = [diameter[g] for g in g_names]
bp = ax_c.boxplot(bp_data, positions=positions, widths=0.5,
                   patch_artist=True, showfliers=True,
                   medianprops=dict(color='#333333', lw=1.5),
                   flierprops=dict(marker='o', markersize=3, alpha=0.4))

for patch, c in zip(bp['boxes'], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.7)
    patch.set_edgecolor(c)
    patch.set_linewidth(1.2)

# Add individual data points with jitter
for i, g in enumerate(g_names):
    jitter = rng.uniform(-0.08, 0.08, len(diameter[g]))
    ax_c.scatter(np.full(len(diameter[g]), i) + jitter, diameter[g],
                 c=colors[i], s=5, alpha=0.35, zorder=5)

# Significance
ax_c.plot([0, 0, 1, 1], [21, 21.3, 21.3, 21], color='#555555', lw=0.6)
ax_c.text(0.5, 21.5, '***', ha='center', fontsize=8, fontweight='bold', color='#D41159')
ax_c.plot([1, 1, 2, 2], [20.5, 20.8, 20.8, 20.5], color='#555555', lw=0.6)
ax_c.text(1.5, 21.0, '*', ha='center', fontsize=8, fontweight='bold', color='#D41159')

ax_c.set_xticks(range(4))
ax_c.set_xticklabels(groups, fontsize=6)
ax_c.set_ylabel('Myotube diameter (um)', fontsize=7)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.grid(axis='y', alpha=0.25, lw=0.4)

# ---- Panel D: Representative IF images (schematic) ----
ax_d = fig.add_subplot(gs_outer[1, :])
ax_d.text(-0.03, 1.05, 'D', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_d.transAxes)
ax_d.set_title('Representative MyHC immunofluorescence (C2C12 myotubes)', fontsize=8,
               fontweight='bold', pad=8)
ax_d.set_xlim(0, 12)
ax_d.set_ylim(0, 3.5)
ax_d.axis('off')

# Draw 4 representative panels
panel_x = [0.3, 3.2, 6.1, 9.0]
for i, (g, c) in enumerate(zip(g_names, colors)):
    x0, y0 = panel_x[i], 0.3
    pw, ph = 2.2, 2.2
    # Panel border
    rect = FancyBboxPatch((x0, y0), pw, ph, boxstyle="round,pad=0.08",
                          facecolor='white', edgecolor=c, linewidth=1.5)
    ax_d.add_patch(rect)
    # Draw myotubes as elongated shapes
    n_tubes = 5 if 'KD-CM' not in g or 'Pro' in g else 3
    for t in range(n_tubes):
        tx = x0 + 0.3 + rng.uniform(0, 0.4)
        ty = y0 + 0.3 + t * (1.6 / max(n_tubes, 1))
        tw = 1.3 + rng.uniform(0, 0.3)
        th = 0.12 + rng.uniform(0, 0.06)
        # MyHC staining (green)
        tube = FancyBboxPatch((tx, ty), tw, th, boxstyle="round,pad=0.02",
                              facecolor='#009E73' if 'Pro' in g or 'NC-CM' in g else '#D55E00',
                              edgecolor='none', alpha=0.6 + (0.1 * (4 - t)))
        ax_d.add_patch(tube)
        # Nuclei (blue dots)
        for nuc_x in [tx + 0.3, tx + 0.7, tx + 1.0]:
            if nuc_x < tx + tw - 0.2:
                ax_d.scatter(nuc_x, ty + th/2, s=25, c='#0072B2', alpha=0.7, zorder=10)

    # Label
    ax_d.text(x0 + pw/2, y0 - 0.25, g.replace('\n', ' '), ha='center', fontsize=6.5,
              fontweight='bold', color=c)

# Legend for staining
ax_d.text(0.5, 3.0, 'MyHC (green) + DAPI (nuclei, blue)', fontsize=6, color='#555555',
          fontstyle='italic')

# Scale bar
ax_d.plot([10.5, 11.5], [0.1, 0.1], color='#333333', lw=2)
ax_d.text(11.0, 0.25, '50 um', ha='center', fontsize=6)

# ---------------------------------------------------------------------------
outpath = './fig4_crosstalk_myotubes.png'
# Data source annotation
fig.text(0.5, 0.008, 'Data: All values = expected patterns from conditioned medium crosstalk experiment (Exp 2.2 pending). SUnSET method per Schmidt et al. 2009. n=3 independent CM preparations. Representative IF images = schematic illustrations.', ha='center', fontsize=5.5, color='#999999', fontstyle='italic')
fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {outpath}')
