#!/usr/bin/env python3
"""
Fig 1: Liver ALDH18A1 branch-point gene mRNA expression and serum metabolic fingerprint.
Panel A: Liver mRNA log2FC (DLY/TFB) for ALDH18A1 and urea cycle genes across 4 stages.
         REAL RNA-seq data from DLY_TFB_肝脏靶点完整注释表.xlsx.
Panel B: Serum metabolite timecourse — Glu, Orn, Cit, Urea across 4 growth stages.
         REAL UPLC-MS/MS data from serum index 0514.xlsx (n=8/group/stage).

DATA: All values in Panel A are measured RNA-seq log2FC. All values in Panel B
are measured serum concentrations. Serum Pro is pending Exp 1.1.

Reference style: Journal of Animal Science and Biotechnology.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.titlesize': 8.5,
    'axes.labelsize': 7.5,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 6.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

C_DLY  = '#D55E00'
C_TFB  = '#0072B2'
C_STAR = '#D41159'
C_PEND = '#999999'

stages = ['15', '45', '75', '105']
x_stage = np.arange(4)

def sig_str(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'ns'

# ===========================================================================
# Panel A: Liver mRNA log2FC — REAL RNA-seq data
# ===========================================================================
# From DLY_TFB_肝脏靶点完整注释表.xlsx. log2FC(DLY/TFB): + = DLY higher.
# * p<0.05, dagger p<0.1 (as marked in source Excel).

mrna_data = {
    'ALDH18A1': {
        'log2FC': np.array([1.61, -1.88, 0.85, 0.76]),
        'sig':    ['*', '*', '*', '†'],
        'color':  C_DLY,
        'label':  'DLY-high (branch-point enzyme)',
    },
    'OAT': {
        'log2FC': np.array([-2.37, -4.73, 0.47, -1.40]),
        'sig':    ['*', '*', 'ns', '*'],
        'color':  C_TFB,
        'label':  'TFB-high (Orn→Pro shunt)',
    },
    'ARG1': {
        'log2FC': np.array([-1.78, -3.64, -0.57, -0.80]),
        'sig':    ['*', '*', 'ns', '*'],
        'color':  '#009E73',
        'label':  'TFB-high (urea cycle)',
    },
    'ASL': {
        'log2FC': np.array([-1.19, -1.55, -0.45, -1.10]),
        'sig':    ['*', '*', '*', '*'],
        'color':  '#E69F00',
        'label':  'TFB-high (urea cycle)',
    },
    'ASS1': {
        'log2FC': np.array([-1.30, -4.61, 0.45, -0.93]),
        'sig':    ['*', '*', 'ns', '*'],
        'color':  '#CC79A7',
        'label':  'TFB-high (urea cycle)',
    },
}

# ===========================================================================
# Panel B: Serum metabolites — REAL UPLC-MS/MS measured data
# ===========================================================================
serum_data = {
    'Glu': {
        'DLY':     np.array([0.32, 0.34, 0.32, 0.29]),
        'TFB':     np.array([0.31, 0.24, 0.25, 0.26]),
        'DLY_sem': np.array([0.02, 0.02, 0.02, 0.02]),
        'TFB_sem': np.array([0.03, 0.02, 0.02, 0.03]),
        'p_vals':  [0.8457, 0.0066, 0.0160, 0.4387],
        'ylabel':  'Serum Glu (mM)',
    },
    'Orn': {
        'DLY':     np.array([0.11, 0.09, 0.07, 0.09]),
        'TFB':     np.array([0.12, 0.10, 0.10, 0.10]),
        'DLY_sem': np.array([0.01, 0.00, 0.01, 0.01]),
        'TFB_sem': np.array([0.02, 0.01, 0.01, 0.01]),
        'p_vals':  [0.6758, 0.3765, 0.0329, 0.7307],
        'ylabel':  'Serum Orn (mM)',
    },
    'Cit': {
        'DLY':     np.array([0.12, 0.09, 0.08, 0.08]),
        'TFB':     np.array([0.13, 0.13, 0.09, 0.11]),
        'DLY_sem': np.array([0.01, 0.01, 0.01, 0.01]),
        'TFB_sem': np.array([0.01, 0.01, 0.01, 0.01]),
        'p_vals':  [0.2027, 0.0008, 0.2759, 0.0565],
        'ylabel':  'Serum Cit (mM)',
    },
    'Urea': {
        'DLY':     np.array([0.81, 2.30, 3.24, 3.93]),
        'TFB':     np.array([3.16, 5.02, 3.31, 4.34]),
        'DLY_sem': np.array([0.11, 0.20, 0.18, 0.25]),
        'TFB_sem': np.array([0.23, 0.24, 0.27, 0.52]),
        'p_vals':  [0.0000, 0.0000, 0.8220, 0.4927],
        'ylabel':  'Serum Urea (mM)',
    },
}
metab_names = ['Glu', 'Orn', 'Cit', 'Urea']

# ===========================================================================
# Layout
# ===========================================================================
fig = plt.figure(figsize=(11, 6.8))

gs_outer = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[0.48, 0.52],
                              left=0.07, right=0.97, top=0.96, bottom=0.08, hspace=0.35)

# ---- Panel A: Liver mRNA ----
gs_a = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs_outer[0], wspace=0.50)
ax_a_label = fig.add_subplot(gs_outer[0])
ax_a_label.axis('off')
ax_a_label.text(-0.04, 1.02, 'A', fontsize=11, fontweight='bold', va='bottom', ha='left',
                transform=ax_a_label.transAxes)
ax_a_label.text(0.5, 1.02, 'Liver mRNA log2FC (DLY/TFB) — RNA-seq measured',
                fontsize=9, fontweight='bold', ha='center', transform=ax_a_label.transAxes)

for i, (gene, d) in enumerate(mrna_data.items()):
    ax = fig.add_subplot(gs_a[i])
    fc = d['log2FC']
    sig_flags = d['sig']
    c = d['color']

    colors_bar = [c if v > 0 else c for v in fc]
    bars = ax.bar(x_stage, fc, width=0.55, color=c, alpha=0.75, edgecolor='white', linewidth=0.5)

    # Zero line
    ax.axhline(y=0, color='#333333', lw=0.8, alpha=0.5)

    # Significance markers
    for j, (v, s) in enumerate(zip(fc, sig_flags)):
        if s != 'ns':
            y_offset = 0.15 if v >= 0 else -0.15
            ax.text(j, v + y_offset, s, ha='center', fontsize=6, fontweight='bold',
                    color=C_STAR)

    ax.set_xticks(x_stage)
    ax.set_xticklabels(stages, fontsize=6.5)
    ax.set_xlabel('Body weight (kg)', fontsize=6.5)
    ax.set_ylabel('log2FC (DLY/TFB)', fontsize=6.5)
    ax.set_title(gene, fontsize=8, fontweight='bold', color=c, pad=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.25, lw=0.4)

    # Set symmetric y-limits around zero for consistent visual comparison
    y_abs_max = max(abs(fc)) * 1.3
    ax.set_ylim(-y_abs_max, y_abs_max)

# ---- Panel B: Serum metabolites ----
gs_b = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=gs_outer[1], wspace=0.50)

ax_b_label = fig.add_subplot(gs_outer[1])
ax_b_label.axis('off')
ax_b_label.text(-0.04, 1.02, 'B', fontsize=11, fontweight='bold', va='bottom', ha='left',
                transform=ax_b_label.transAxes)
ax_b_label.text(0.5, 1.02, 'Serum metabolite timecourse — UPLC-MS/MS measured (n=8/group)',
                fontsize=9, fontweight='bold', ha='center', transform=ax_b_label.transAxes)

for i, name in enumerate(metab_names):
    ax = fig.add_subplot(gs_b[i])
    d = serum_data[name]

    ax.errorbar(x_stage, d['DLY'], yerr=d['DLY_sem'], color=C_DLY, lw=1.8, marker='o', ms=5,
                capsize=3, capthick=0.8, label='DLY')
    ax.errorbar(x_stage, d['TFB'], yerr=d['TFB_sem'], color=C_TFB, lw=1.8, marker='s', ms=5,
                capsize=3, capthick=0.8, label='TFB')

    for j, p in enumerate(d['p_vals']):
        if p < 0.05:
            y_range = max(d['DLY']) - min(min(d['DLY']), min(d['TFB']))
            y_pos = max(d['DLY'][j], d['TFB'][j]) + max(d['DLY_sem'][j], d['TFB_sem'][j]) + y_range * 0.08
            ax.text(j, y_pos, sig_str(p), ha='center', fontsize=6, fontweight='bold',
                    color=C_STAR)

    ax.set_xticks(x_stage)
    ax.set_xticklabels(stages, fontsize=6.5)
    ax.set_xlabel('Body weight (kg)', fontsize=6.5)
    ax.set_ylabel(d['ylabel'], fontsize=6.5)
    ax.set_title(name, fontsize=8, fontweight='bold', color='#333333', pad=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.25, lw=0.4)
    if i == 0:
        ax.legend(frameon=False, fontsize=6, loc='upper left')

# ---- 5th panel: Serum Pro — NOT MEASURED ----
ax_pro = fig.add_subplot(gs_b[4])
ax_pro.set_title('Pro', fontsize=8, fontweight='bold', color=C_PEND, pad=3)
ax_pro.set_xticks(x_stage)
ax_pro.set_xticklabels(stages, fontsize=6.5)
ax_pro.set_xlabel('Body weight (kg)', fontsize=6.5)
ax_pro.set_ylabel('Serum Pro (uM)', fontsize=6.5, color=C_PEND)
ax_pro.set_ylim(0, 300)
ax_pro.set_xlim(-0.4, 3.4)
ax_pro.spines['top'].set_visible(False)
ax_pro.spines['right'].set_visible(False)
ax_pro.grid(alpha=0.25, lw=0.4)
ax_pro.text(1.5, 150, 'Pending measurement\n(Exp 1.1)', ha='center', va='center',
            fontsize=7, fontweight='bold', color=C_PEND, fontstyle='italic',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#F5F5F5', edgecolor=C_PEND,
                      lw=0.8, alpha=0.85))
for y_line in [240, 170]:
    ax_pro.axhline(y=y_line, color=C_PEND, lw=0.5, ls='--', alpha=0.4)
ax_pro.text(3.2, 243, 'DLY\n(expected)', fontsize=5, color=C_DLY, ha='center', alpha=0.5)
ax_pro.text(3.2, 167, 'TFB\n(expected)', fontsize=5, color=C_TFB, ha='center', alpha=0.5)

# ===========================================================================
fig.text(0.5, 0.012,
         'Panel A: Liver mRNA log2FC(DLY/TFB) from RNA-seq (DLY_TFB_肝脏靶点完整注释表.xlsx). '
         '* p<0.05, dagger p<0.1. ALDH18A1 is DLY-high; OAT/ARG1/ASL/ASS1 are TFB-high. '
         'PYCR1 and CPS1: not measured in this RNA-seq dataset. '
         'Panel B: Serum Glu/Orn/Cit/Urea measured by UPLC-MS/MS (n=8/group/stage), '
         'source: serum index 0514.xlsx. '
         'Serum Pro: not measured — Exp 1.1 pending. '
         'Error bars = SEM. p-values from Welch\'s t-test.',
         ha='center', fontsize=5.5, color='#999999', fontstyle='italic')

outpath = './fig1_liver_serum_validation.png'
fig.savefig(outpath, dpi=300, facecolor='white', edgecolor='none')
plt.close(fig)
print(f'Saved: {outpath}')
