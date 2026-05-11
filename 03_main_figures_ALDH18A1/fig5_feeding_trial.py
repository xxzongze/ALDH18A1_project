#!/usr/bin/env python3
"""
Fig 5: Dietary Pro+Glu supplementation increases nitrogen retention and protein deposition —
a 2x2 factorial trial (Breed x Diet).
Panel A: Average daily gain (ADG).
Panel B: Nitrogen retention (g/d).
Panel C: Serum Pro timecourse (day 0, 28, 42, 56).
Panel D: Serum Urea timecourse.
Panel E: Interaction plot — Breed x Diet on ADG.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.titlesize': 8,
    'axes.labelsize': 7,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

rng = np.random.default_rng(2025)
n_per_group = 8

C_DLY_CON  = '#D55E00'
C_DLY_SUPP = '#E69F00'
C_TFB_CON  = '#0072B2'
C_TFB_SUPP = '#56B4E9'

group_names = ['DLY\nCON', 'DLY\nSUPP', 'TFB\nCON', 'TFB\nSUPP']
group_colors = [C_DLY_CON, C_DLY_SUPP, C_TFB_CON, C_TFB_SUPP]
g_keys = ['DLY-CON', 'DLY-SUPP', 'TFB-CON', 'TFB-SUPP']

days = [0, 28, 42, 56]

def sig_str(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'ns'

# ---------------------------------------------------------------------------
# Simulated data
# ---------------------------------------------------------------------------

# -- ADG (g/d) --
adg = {
    'DLY-CON':  rng.normal(720, 45, n_per_group),
    'DLY-SUPP': rng.normal(745, 48, n_per_group),
    'TFB-CON':  rng.normal(580, 50, n_per_group),
    'TFB-SUPP': rng.normal(655, 52, n_per_group),
}

# -- N retention (g/d) at end of trial --
n_ret = {
    'DLY-CON':  rng.normal(28.5, 2.5, n_per_group),
    'DLY-SUPP': rng.normal(29.8, 2.6, n_per_group),
    'TFB-CON':  rng.normal(22.0, 2.8, n_per_group),
    'TFB-SUPP': rng.normal(26.5, 2.7, n_per_group),
}

# -- Serum Pro (uM) timecourse --
pro_tc = {
    'DLY-CON':  [rng.normal(185, 12, n_per_group) for _ in range(4)],
    'DLY-SUPP': [rng.normal(185, 12, n_per_group),
                 rng.normal(220, 15, n_per_group),
                 rng.normal(238, 16, n_per_group),
                 rng.normal(245, 17, n_per_group)],
    'TFB-CON':  [rng.normal(155, 11, n_per_group) for _ in range(4)],
    'TFB-SUPP': [rng.normal(155, 11, n_per_group),
                 rng.normal(200, 14, n_per_group),
                 rng.normal(218, 15, n_per_group),
                 rng.normal(225, 16, n_per_group)],
}

# -- Serum Urea (mM) timecourse --
urea_tc = {
    'DLY-CON':  [rng.normal(3.0, 0.3, n_per_group) for _ in range(4)],
    'DLY-SUPP': [rng.normal(3.0, 0.3, n_per_group),
                 rng.normal(2.7, 0.28, n_per_group),
                 rng.normal(2.5, 0.25, n_per_group),
                 rng.normal(2.4, 0.25, n_per_group)],
    'TFB-CON':  [rng.normal(4.0, 0.35, n_per_group) for _ in range(4)],
    'TFB-SUPP': [rng.normal(4.0, 0.35, n_per_group),
                 rng.normal(3.2, 0.30, n_per_group),
                 rng.normal(2.9, 0.28, n_per_group),
                 rng.normal(2.7, 0.28, n_per_group)],
}

# P-values for annotation
p_breed_adg = 0.0001
p_diet_adg  = 0.015
p_inter_adg = 0.042
p_breed_n   = 0.0005
p_diet_n    = 0.008
p_inter_n   = 0.035

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(10, 8.5))
gs_outer = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 0.9],
                              left=0.07, right=0.97, top=0.95, bottom=0.06,
                              hspace=0.45, wspace=0.50)

# ---- Panel A: ADG ----
ax_a = fig.add_subplot(gs_outer[0, 0])
ax_a.text(-0.15, 1.05, 'A', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_a.transAxes)
ax_a.set_title('Average Daily Gain', fontsize=9, fontweight='bold', pad=4)

for i, g in enumerate(g_keys):
    d = adg[g]
    m, s = np.mean(d), np.std(d, ddof=1) / np.sqrt(n_per_group)
    ax_a.bar(i, m, 0.55, color=group_colors[i], edgecolor='white', linewidth=0.5)
    ax_a.errorbar(i, m, yerr=s, color='#333333', capsize=3, lw=0.8)
    for v in d:
        ax_a.scatter(i + rng.uniform(-0.07, 0.07), v, color=group_colors[i], s=8, alpha=0.4, zorder=5)

# Annotate p-values
ax_a.text(0.5, 820, f'Breed: {sig_str(p_breed_adg)}\nDiet: {sig_str(p_diet_adg)}\nBxD: {sig_str(p_inter_adg)}',
          fontsize=7, ha='center', va='top', color='#333333',
          bbox=dict(boxstyle='round,pad=0.35', facecolor='#FFF8E1', edgecolor='#E69F00', alpha=0.8, lw=0.6))

ax_a.set_xticks(range(4))
ax_a.set_xticklabels(group_names, fontsize=6.5)
ax_a.set_ylabel('ADG (g/d)', fontsize=7.5)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.grid(axis='y', alpha=0.25, lw=0.4)

# ---- Panel B: N retention ----
ax_b = fig.add_subplot(gs_outer[0, 1])
ax_b.text(-0.15, 1.05, 'B', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_b.transAxes)
ax_b.set_title('Nitrogen Retention', fontsize=9, fontweight='bold', pad=4)

for i, g in enumerate(g_keys):
    d = n_ret[g]
    m, s = np.mean(d), np.std(d, ddof=1) / np.sqrt(n_per_group)
    ax_b.bar(i, m, 0.55, color=group_colors[i], edgecolor='white', linewidth=0.5)
    ax_b.errorbar(i, m, yerr=s, color='#333333', capsize=3, lw=0.8)
    for v in d:
        ax_b.scatter(i + rng.uniform(-0.07, 0.07), v, color=group_colors[i], s=8, alpha=0.4, zorder=5)

ax_b.text(0.5, 34, f'Breed: {sig_str(p_breed_n)}\nDiet: {sig_str(p_diet_n)}\nBxD: {sig_str(p_inter_n)}',
          fontsize=7, ha='center', va='top', color='#333333',
          bbox=dict(boxstyle='round,pad=0.35', facecolor='#FFF8E1', edgecolor='#E69F00', alpha=0.8, lw=0.6))

ax_b.set_xticks(range(4))
ax_b.set_xticklabels(group_names, fontsize=6.5)
ax_b.set_ylabel('N retention (g/d)', fontsize=7.5)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.grid(axis='y', alpha=0.25, lw=0.4)

# ---- Panel C: Serum Pro timecourse ----
ax_c = fig.add_subplot(gs_outer[1, 0])
ax_c.text(-0.15, 1.05, 'C', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_c.transAxes)
ax_c.set_title('Serum Proline Timecourse', fontsize=9, fontweight='bold', pad=4)

for g, c, ls in zip(g_keys, group_colors, ['-', '-', '--', '--']):
    means = [np.mean(d) for d in pro_tc[g]]
    sems  = [np.std(d, ddof=1)/np.sqrt(n_per_group) for d in pro_tc[g]]
    marker = 'o' if 'DLY' in g else 's'
    ax_c.errorbar(days, means, yerr=sems, color=c, lw=1.8, marker=marker, ms=5,
                  capsize=3, capthick=0.8, ls=ls, label=g)

ax_c.set_xlabel('Day', fontsize=7)
ax_c.set_ylabel('Serum Pro (uM)', fontsize=7)
ax_c.legend(frameon=False, fontsize=6, ncol=2)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.grid(alpha=0.25, lw=0.4)
ax_c.set_xticks(days)

# ---- Panel D: Serum Urea timecourse ----
ax_d = fig.add_subplot(gs_outer[1, 1])
ax_d.text(-0.15, 1.05, 'D', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_d.transAxes)
ax_d.set_title('Serum Urea Timecourse', fontsize=9, fontweight='bold', pad=4)

for g, c, ls in zip(g_keys, group_colors, ['-', '-', '--', '--']):
    means = [np.mean(d) for d in urea_tc[g]]
    sems  = [np.std(d, ddof=1)/np.sqrt(n_per_group) for d in urea_tc[g]]
    marker = 'o' if 'DLY' in g else 's'
    ax_d.errorbar(days, means, yerr=sems, color=c, lw=1.8, marker=marker, ms=5,
                  capsize=3, capthick=0.8, ls=ls, label=g)

ax_d.set_xlabel('Day', fontsize=7)
ax_d.set_ylabel('Serum Urea (mM)', fontsize=7)
ax_d.legend(frameon=False, fontsize=6, ncol=2)
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)
ax_d.grid(alpha=0.25, lw=0.4)
ax_d.set_xticks(days)

# ---- Panel E: Interaction plot (Breed x Diet on ADG) ----
ax_e = fig.add_subplot(gs_outer[2, :])
ax_e.text(-0.03, 1.05, 'E', fontsize=11, fontweight='bold', va='bottom', ha='left',
          transform=ax_e.transAxes)
ax_e.set_title('Breed x Diet Interaction on ADG', fontsize=9, fontweight='bold', pad=4)

# Interaction plot
breed_labels = ['DLY', 'TFB']
x_inter = [0, 1]
# DLY
dly_con_mean = np.mean(adg['DLY-CON'])
dly_supp_mean = np.mean(adg['DLY-SUPP'])
dly_con_sem  = np.std(adg['DLY-CON'], ddof=1) / np.sqrt(n_per_group)
dly_supp_sem  = np.std(adg['DLY-SUPP'], ddof=1) / np.sqrt(n_per_group)
# TFB
tfb_con_mean = np.mean(adg['TFB-CON'])
tfb_supp_mean = np.mean(adg['TFB-SUPP'])
tfb_con_sem  = np.std(adg['TFB-CON'], ddof=1) / np.sqrt(n_per_group)
tfb_supp_sem  = np.std(adg['TFB-SUPP'], ddof=1) / np.sqrt(n_per_group)

# DLY line
ax_e.errorbar([0, 1], [dly_con_mean, dly_supp_mean],
              yerr=[dly_con_sem, dly_supp_sem],
              color=C_DLY_CON, lw=2.5, marker='o', ms=9, capsize=4, capthick=1.2,
              label='DLY')
# TFB line
ax_e.errorbar([0, 1], [tfb_con_mean, tfb_supp_mean],
              yerr=[tfb_con_sem, tfb_supp_sem],
              color=C_TFB_CON, lw=2.5, marker='s', ms=9, capsize=4, capthick=1.2,
              label='TFB')

# Fill to emphasize divergence
ax_e.fill_between([0, 1], [dly_con_mean - dly_con_sem, dly_supp_mean - dly_supp_sem],
                  [dly_con_mean + dly_con_sem, dly_supp_mean + dly_supp_sem],
                  color=C_DLY_CON, alpha=0.1)
ax_e.fill_between([0, 1], [tfb_con_mean - tfb_con_sem, tfb_supp_mean - tfb_supp_sem],
                  [tfb_con_mean + tfb_con_sem, tfb_supp_mean + tfb_supp_sem],
                  color=C_TFB_CON, alpha=0.1)

ax_e.set_xticks([0, 1])
ax_e.set_xticklabels(['CON (Basal)', 'SUPP (+0.5%Pro+1%Glu)'], fontsize=7.5)
ax_e.set_ylabel('ADG (g/d)', fontsize=7.5)
ax_e.legend(frameon=False, fontsize=7.5)
ax_e.spines['top'].set_visible(False)
ax_e.spines['right'].set_visible(False)
ax_e.grid(axis='y', alpha=0.25, lw=0.4)

# Interactive p-value annotation
ax_e.annotate(f'P_interaction = 0.042', xy=(0.5, 0.15), xycoords='axes fraction',
              fontsize=8, fontweight='bold', color='#D41159', ha='center',
              bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#D41159', lw=1))

# ---------------------------------------------------------------------------
outpath = '/Users/hezongze/ALDH18A1_project/03_main_figures_ALDH18A1/fig5_feeding_trial.png'
# Data source annotation
fig.text(0.5, 0.008, 'Data: All values = expected patterns from 2x2 factorial feeding trial (Exp 3.1 pending). n=8/group, 25-75 kg, ~60 days. Statistical model: lmer(y ~ Breed*Diet + (1|Block)). Interaction p-values are hypothetical and for illustration only.', ha='center', fontsize=5.5, color='#999999', fontstyle='italic')
fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {outpath}')
