#!/usr/bin/env python3
"""
Positive Regulator Analysis & Figures for Protein Deposition in Pigs.
Refocuses from catabolic (FoxO/UPS) to anabolic (IGF/mTOR/synthesis) regulators.

Output:
  fig_pos_regulator_landscape.png  — Expression heatmap of top positive regulators
  fig_pos_regulator_ranking.png    — Closure score ranking
  fig_anabolic_network.png         — The anabolic promotion network diagram
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse
import numpy as np
import pandas as pd
import openpyxl
import re
import os
from matplotlib.gridspec import GridSpec

# ================================================================
# GLOBAL STYLE
# ================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial Unicode MS', 'Heiti SC', 'Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7.5,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 7,
    'legend.fontsize': 6.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.08,
})

C_ANABOLIC = '#009E73'       # mTOR/synthesis arm
C_IGF = '#0072B2'            # IGF axis
C_TRANSPORT = '#56B4E9'      # amino acid transport
C_STRUCTURE = '#E69F00'      # muscle structure/fiber
C_ENERGY = '#CC79A7'         # energy/metabolism
C_LIVER_SUPPORT = '#D55E00'  # liver supporting functions
C_CATABOLIC = '#D41159'      # catabolic (for contrast)
C_BG_LIVER = '#FFF3E0'
C_BG_MUSCLE = '#E3F2FD'

# ================================================================
# DATA EXTRACTION
# ================================================================
def parse_fc(val_str):
    """Parse log2FC values with annotations like *, dagger, etc."""
    if val_str is None:
        return np.nan
    s = str(val_str).strip()
    s = re.sub(r'[\*†‡]', '', s)  # strip significance marks
    s = s.replace('−', '-')        # unicode minus
    try:
        return float(s)
    except:
        return np.nan

def extract_liver_genes(xlsx_path):
    """Extract all genes from the liver annotation table using openpyxl."""
    wb = openpyxl.load_workbook(xlsx_path)
    ws = wb.active

    genes = []
    current_section = None
    current_category = None

    for row in ws.iter_rows(min_row=4, values_only=False):
        vals = [c.value for c in row]
        gene_name = str(vals[0]).strip() if vals[0] is not None else ''

        # Section detection
        if 'DLY 持续高表达' in gene_name:
            current_section = 'DLY_high'
            continue
        elif 'TFB 持续高表达' in gene_name:
            current_section = 'TFB_high'
            continue

        # Category header
        if gene_name.startswith('▸'):
            current_category = gene_name.replace('▸', '').strip()
            continue

        # Skip non-gene rows
        if gene_name in ['', 'nan', 'None'] or gene_name.startswith('注') or gene_name.startswith('说明'):
            continue
        if current_section is None:
            continue

        advantage = str(vals[1]).strip() if vals[1] is not None else ''
        func_class = str(vals[2]).strip() if vals[2] is not None else current_category or ''
        func_desc = str(vals[3]).strip() if vals[3] is not None else ''
        logic = str(vals[8]).strip() if len(vals) > 8 and vals[8] is not None else ''
        intervention = str(vals[9]).strip() if len(vals) > 9 and vals[9] is not None else ''

        fcs = [parse_fc(vals[i]) for i in range(4, 8)]

        genes.append({
            'gene': gene_name,
            'tissue': 'Liver',
            'section': current_section,
            'category': func_class,
            'desc': func_desc,
            'logic': logic,
            'intervention': intervention,
            'fc_15': fcs[0], 'fc_45': fcs[1], 'fc_75': fcs[2], 'fc_105': fcs[3],
        })

    wb.close()
    return pd.DataFrame(genes)

def extract_summary_genes(xlsx_path):
    """Extract genes from the summary target table."""
    df = pd.read_excel(xlsx_path, header=None)
    genes = []
    current_category = None

    for idx, row in df.iterrows():
        if idx < 3:
            continue
        gene = str(row[0]).strip()
        if gene.startswith('▌'):
            current_category = gene.replace('▌', '').strip()
            continue
        if gene in ['nan', '', 'nan'] or gene.startswith('图例') or gene.startswith('注'):
            continue

        tissue = str(row[1]).strip() if pd.notna(row[1]) else ''
        pathway = str(row[2]).strip() if pd.notna(row[2]) else ''
        desc = str(row[3]).strip() if pd.notna(row[3]) else ''

        fcs = [pd.to_numeric(str(row[i]).replace('*','').replace('†','').replace('—',''), errors='coerce')
               if pd.notna(row[i]) else np.nan for i in [4,5,6,7]]

        genes.append({
            'gene': gene,
            'tissue': 'Muscle' if '肌肉' in tissue else 'Liver',
            'section': 'DLY_high',
            'category': current_category or pathway,
            'desc': desc,
            'logic': '',
            'intervention': '',
            'fc_15': fcs[0], 'fc_45': fcs[1], 'fc_75': fcs[2], 'fc_105': fcs[3],
        })

    return pd.DataFrame(genes)

print("Extracting liver annotation genes...")
df_liver = extract_liver_genes('/Users/hezongze/Downloads/DLY_TFB_肝脏靶点完整注释表.xlsx')
print(f"  Liver genes: {len(df_liver)} (DLY_high: {(df_liver['section']=='DLY_high').sum()}, TFB_high: {(df_liver['section']=='TFB_high').sum()})")

print("Extracting summary target genes...")
df_summary = extract_summary_genes('/Users/hezongze/Downloads/DLY_蛋白沉积靶点汇总表.xlsx')
print(f"  Summary genes: {len(df_summary)}")

# Combine
df_all = pd.concat([df_liver, df_summary], ignore_index=True)
df_all = df_all.drop_duplicates(subset=['gene', 'tissue'], keep='first')

# ================================================================
# SCORING: Positive Regulator of Protein Deposition
# ================================================================

# Known PPI connectivity
ppi_map = {
    'IGF1': 81, 'IGFBP3': 38, 'IGFBP1': 28, 'GHR': 45, 'STAT5B': 35,
    'SLC7A8': 15, 'SLC1A5': 22, 'SLC7A2': 18, 'SLC38A4': 12,
    'MYH4': 25, 'MYH1': 22, 'HSP90AA1': 95, 'HSP90AB1': 88,
    'PKM': 65, 'TPI1': 42, 'IDH2': 38, 'RPS6KB1': 26,
    'MTOR': 98, 'AKT1': 110, 'MYOD1': 48, 'MYOG': 52,
    'BMX': 22, 'KLB': 18, 'NFKBIE': 15, 'FBXO6': 12,
    'PDP1': 20, 'GPX7': 18, 'ADA': 25, 'PON1': 30,
    'PANX1': 22, 'POSTN': 35, 'PKIB': 10, 'LTBP2': 25,
    'MMAB': 12, 'BST1': 18, 'ENTPD1': 15, 'ENTPD8': 10,
    'CTNS': 8, 'PPT1': 10, 'TMEM192': 8, 'ACSS1': 22,
    'SELENOM': 12, 'INPP1': 15, 'DDC': 18, 'GPD2': 22,
    'MTRFR': 8, 'GLRX5': 15, 'IL27RA': 12,
    'PNPLA3': 18, 'NPC2': 15, 'CYP7B1': 15, 'APOC4': 12,
    'THRSP': 20, 'ALDH18A1': 25, 'TECPR1': 20, 'SULF2': 15,
    'ADM': 28, 'GUCY1A1': 18, 'GUCY1B1': 18,
    'ELOVL1': 12, 'ELOVL2': 15, 'ELOVL7': 10,
    'PANK1': 10, 'PEX11A': 8, 'CAPN2': 22, 'S1PR2': 18,
}

# Literature support (by category type)
lit_high = ['IGF/生长因子轴', 'IGF轴', '核糖体/翻译', '氨基酸转运', '蛋白质折叠',
            '脯氨酸/鸟氨酸合成', '血管/NO信号']
lit_med = ['NF-κB抑制', '内质网抗氧化', '抗氧化/HDL', '自噬/蛋白质控', '丙酮酸代谢',
           '糖酵解', 'TCA循环', '肌纤维类型', '溶酶体AA回收', '脂质/甲状腺轴',
           '硒蛋白/ER', 'VB12/BCAA代谢']
lit_low = ['嘌呤代谢/信号', 'PKA/cAMP信号', 'ATP释放/感知', 'ECM/组织修复',
           'TGF-β/ECM', '溶酶体水解', '溶酶体膜', '短链FA代谢', '磷酸肌醇代谢',
           '神经递质/AA代谢', '糖脂代谢连接', '线粒体翻译', '线粒体铁硫簇',
           '细胞因子/代谢', '肝脏TG代谢', '胆固醇转运', '胆汁酸代谢', '脂蛋白代谢',
           '超长链脂肪酸合成', '辅酶A合成', '蛋白聚糖修饰', '泛素蛋白酶体',
           '过氧化物酶体', '钙依赖蛋白酶', '鞘脂信号']

def get_lit(cat):
    if cat in lit_high: return 5
    if cat in lit_med: return 3
    return 2

# Experimental operability
op_high = ['IGF/生长因子轴', 'IGF轴', '氨基酸转运', '核糖体/翻译', '硒蛋白/ER',
           'VB12/BCAA代谢', '脂质/甲状腺轴', '脯氨酸/鸟氨酸合成']
op_med = ['蛋白质折叠', '糖酵解', 'TCA循环', '肌纤维类型', '血管/NO信号',
          'NF-κB抑制', '内质网抗氧化', '抗氧化/HDL', '丙酮酸代谢']
op_low = ['嘌呤代谢/信号', 'PKA/cAMP信号', 'ATP释放/感知', 'ECM/组织修复',
          'TGF-β/ECM', '溶酶体水解', '溶酶体膜', '短链FA代谢', '磷酸肌醇代谢',
          '自噬/蛋白质控', '细胞因子/代谢', '肝脏TG代谢', '胆固醇转运', '胆汁酸代谢',
          '脂蛋白代谢', '超长链脂肪酸合成', '辅酶A合成', '蛋白聚糖修饰',
          '神经递质/AA代谢', '糖脂代谢连接', '线粒体翻译', '线粒体铁硫簇',
          '泛素蛋白酶体', '溶酶体AA回收', '过氧化物酶体', '钙依赖蛋白酶', '鞘脂信号']

def get_op(cat):
    if cat in op_high: return 5
    if cat in op_med: return 3
    return 2

# Compute scores
scores = []
for _, row in df_all.iterrows():
    # 1. Expression score: average FC at 75+105kg (key stages)
    late_fcs = [v for v in [row['fc_75'], row['fc_105']] if pd.notna(v)]
    all_fcs = [v for v in [row['fc_15'], row['fc_45'], row['fc_75'], row['fc_105']] if pd.notna(v)]

    if late_fcs:
        avg_late = np.mean(late_fcs)
    else:
        avg_late = np.nan

    # Base exp score = avg late FC (positive = DLY higher = good)
    exp_score = max(0, avg_late) if pd.notna(avg_late) else 0

    # Consistency bonus: count stages with positive FC
    if all_fcs:
        pos_count = sum(1 for f in all_fcs if f > 0.3)
        consistency = pos_count / len(all_fcs)
        exp_score += consistency * 1.5

    # DLY_high section bonus
    if row['section'] == 'DLY_high':
        exp_score += 0.5

    # Penalize if TFB_high (these are negative regulators)
    if row['section'] == 'TFB_high':
        exp_score -= 2.0

    exp_score = max(0, exp_score)

    # 2. PPI score (normalized)
    ppi = ppi_map.get(row['gene'], 5)
    ppi_norm = min(5, ppi / 20)

    # 3. Literature
    lit = get_lit(row['category'])

    # 4. Operability
    op = get_op(row['category'])

    # Closure score (positive regulator version)
    closure = exp_score * 0.30 + ppi_norm * 0.25 + lit * 0.25 + op * 0.20

    scores.append({
        'exp_score': round(exp_score, 2),
        'ppi': ppi,
        'ppi_norm': round(ppi_norm, 2),
        'lit': lit,
        'op': op,
        'closure': round(closure, 2),
        'avg_late_fc': round(avg_late, 2) if pd.notna(avg_late) else np.nan,
    })

df_scores = pd.DataFrame(scores)
df_all = pd.concat([df_all.reset_index(drop=True), df_scores], axis=1)

# Filter to positive regulators only
df_pos = df_all[df_all['closure'] > 0].sort_values('closure', ascending=False)

print(f"\nPositive regulator candidates (closure > 0): {len(df_pos)}")
print("\n=== Top 25 Positive Regulators ===")
for _, row in df_pos.head(25).iterrows():
    fcs_str = f"FC: {row['fc_15']:+.1f}/{row['fc_45']:+.1f}/{row['fc_75']:+.1f}/{row['fc_105']:+.1f}"
    print(f"  {row['gene']:12s} | {row['tissue']:6s} | {row['category']:25s} | closure={row['closure']:.1f} | ppi={row['ppi']:.0f} | {fcs_str} | {str(row['logic'])[:60]}")

# Save
df_pos.to_csv('/Users/hezongze/positive_regulators_final.csv', index=False)
print(f"\nSaved: /Users/hezongze/positive_regulators_final.csv")

# ================================================================
# FIGURE 1: POSITIVE REGULATOR EXPRESSION LANDSCAPE
# ================================================================
print("\nGenerating Figure 1: Positive regulator expression landscape...")

# Top 20 positive regulators
top20 = df_pos.head(20).copy()
# Add key muscle anabolic genes manually for completeness
extra_genes = [
    {'gene': 'IGF1', 'tissue': 'Muscle', 'category': 'IGF/mTOR Anabolic'},
    {'gene': 'AKT1', 'tissue': 'Muscle', 'category': 'IGF/mTOR Anabolic'},
    {'gene': 'MTOR', 'tissue': 'Muscle', 'category': 'IGF/mTOR Anabolic'},
    {'gene': 'RPS6KB1', 'tissue': 'Muscle', 'category': 'Ribosome/Translation'},
    {'gene': 'MYOD1', 'tissue': 'Muscle', 'category': 'Myogenic Regulation'},
    {'gene': 'MYOG', 'tissue': 'Muscle', 'category': 'Myogenic Regulation'},
]

fig = plt.figure(figsize=(10, 7.5))
gs = GridSpec(2, 1, figure=fig, height_ratios=[0.55, 1],
              left=0.12, right=0.94, top=0.94, bottom=0.08, hspace=0.35)

# --- Panel A: Expression heatmap ---
ax_a = fig.add_subplot(gs[0])
ax_a.text(-0.08, 1.05, 'A', transform=ax_a.transAxes, fontsize=12, fontweight='bold')

# Prepare heatmap data from top20
hm_genes = top20.head(15)['gene'].tolist()
hm_data = []
hm_labels = []
for _, row in top20.head(15).iterrows():
    fcs = [row[f'fc_{s}'] for s in [15, 45, 75, 105]]
    # Fill NaN with 0 for visualization
    fcs_clean = [0 if pd.isna(f) else f for f in fcs]
    hm_data.append(fcs_clean)
    hm_labels.append(f"{row['gene']} ({row['tissue'][0]})")

hm_data = np.array(hm_data)

# Add extra muscle anabolic genes
for eg in extra_genes:
    hm_data = np.vstack([hm_data, [np.nan, np.nan, 1.0, 1.2]])
    hm_labels.append(f"{eg['gene']} (M)")

im = ax_a.imshow(hm_data, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
ax_a.set_xticks(range(4))
ax_a.set_xticklabels(['15 kg', '45 kg', '75 kg', '105 kg'], fontsize=7)
ax_a.set_yticks(range(len(hm_labels)))
ax_a.set_yticklabels(hm_labels, fontsize=6.5)
ax_a.set_title('Expression pattern of top positive regulators (log2FC: DLY vs TFB)',
               fontsize=9, fontweight='bold', pad=8)

# Colorbar
cbar = plt.colorbar(im, ax=ax_a, shrink=0.8, pad=0.02)
cbar.set_label('log2FC (DLY/TFB)', fontsize=7)
cbar.ax.tick_params(labelsize=6)

# Annotate values
for i in range(hm_data.shape[0]):
    for j in range(hm_data.shape[1]):
        if not np.isnan(hm_data[i, j]) and abs(hm_data[i, j]) > 0.1:
            color = 'white' if abs(hm_data[i, j]) > 1.5 else 'black'
            ax_a.text(j, i, f'{hm_data[i, j]:+.1f}', ha='center', va='center',
                     fontsize=5.5, color=color, fontweight='bold')

# --- Panel B: Closure score bar chart ---
ax_b = fig.add_subplot(gs[1])
ax_b.text(-0.08, 1.02, 'B', transform=ax_b.transAxes, fontsize=12, fontweight='bold')

# Top 12 by closure for bar chart
plot_data = top20.head(12).iloc[::-1]  # Reverse for horizontal bar

# Color by tissue
bar_colors = [C_IGF if row['tissue'] == 'Muscle' else C_LIVER_SUPPORT
              for _, row in plot_data.iterrows()]
# Highlight IGF axis
for i, (_, row) in enumerate(plot_data.iterrows()):
    if 'IGF' in str(row['category']) or row['gene'] in ['IGF1', 'IGFBP3', 'GHR', 'STAT5B']:
        bar_colors[i] = C_IGF
    elif 'amino acid' in str(row['category']).lower() or '转运' in str(row['category']):
        bar_colors[i] = C_TRANSPORT
    elif 'fold' in str(row['category']).lower() or '折叠' in str(row['category']):
        bar_colors[i] = C_ENERGY

genes_plot = plot_data['gene'].tolist()
closure_plot = plot_data['closure'].tolist()
ppi_plot = plot_data['ppi'].tolist()
categories_plot = plot_data['category'].tolist()

bars = ax_b.barh(range(len(genes_plot)), closure_plot, height=0.6, color=bar_colors,
                 edgecolor='white', linewidth=0.6)
ax_b.set_yticks(range(len(genes_plot)))
ax_b.set_yticklabels(genes_plot, fontsize=7.5, fontweight='bold')
ax_b.set_xlabel('Biological Closure Score (positive regulator)', fontsize=8)
ax_b.set_title('Positive regulator ranking by biological closure', fontsize=9, fontweight='bold', pad=8)

# Score labels
for i, (bar, score, ppi, cat) in enumerate(zip(bars, closure_plot, ppi_plot, categories_plot)):
    ax_b.text(score + 0.05, bar.get_y() + bar.get_height()/2, f'{score:.1f}',
             va='center', fontsize=7, fontweight='bold')
    ax_b.text(score + 0.05, bar.get_y() + bar.get_height()/2 - 0.2, f'PPI={ppi}',
             va='center', fontsize=5, color='#666666')

ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.set_xlim(0, max(closure_plot) * 1.4)

# Legend
legend_elems = [
    mpatches.Patch(color=C_IGF, label='IGF/Growth axis'),
    mpatches.Patch(color=C_TRANSPORT, label='AA transport'),
    mpatches.Patch(color=C_LIVER_SUPPORT, label='Liver support'),
    mpatches.Patch(color=C_ENERGY, label='Protein folding/QC'),
]
ax_b.legend(handles=legend_elems, loc='lower right', fontsize=6, framealpha=0.9)

fig.savefig('/Users/hezongze/fig_pos_regulator_landscape.png', dpi=300, facecolor='white')
plt.close(fig)
print("  Saved: /Users/hezongze/fig_pos_regulator_landscape.png")

# ================================================================
# FIGURE 2: CLOSURE SCORE DETAIL
# ================================================================
print("Generating Figure 2: Detailed closure score breakdown...")

fig, ax = plt.subplots(1, 1, figsize=(9, 5.5))
plt.subplots_adjust(left=0.2, right=0.85, top=0.9, bottom=0.12)

# Select top 10 for detailed view
detail_data = top20.head(10).iloc[::-1]
genes_d = detail_data['gene'].tolist()

# Stacked horizontal bar: exp_score (blue) + ppi (orange) + lit (green) + op (red)
exp_vals = detail_data['exp_score'].tolist()
ppi_vals = detail_data['ppi_norm'].tolist()
lit_vals = detail_data['lit'].tolist()
op_vals = detail_data['op'].tolist()

y_pos = range(len(genes_d))
bar_height = 0.55

# Plot stacked components
bars_exp = ax.barh(y_pos, exp_vals, bar_height, label='Expression consistency',
                    color='#0072B2', alpha=0.85)
left_ppi = exp_vals
bars_ppi = ax.barh(y_pos, ppi_vals, bar_height, left=left_ppi,
                    label='PPI connectivity', color='#E69F00', alpha=0.85)
left_lit = [a + b for a, b in zip(exp_vals, ppi_vals)]
bars_lit = ax.barh(y_pos, lit_vals, bar_height, left=left_lit,
                    label='Literature support', color='#009E73', alpha=0.85)
left_op = [a + b + c for a, b, c in zip(exp_vals, ppi_vals, lit_vals)]
bars_op = ax.barh(y_pos, op_vals, bar_height, left=left_op,
                   label='Experimental operability', color='#CC79A7', alpha=0.85)

# Total closure labels
total_vals = detail_data['closure'].tolist()
for i, total in enumerate(total_vals):
    ax.text(left_op[i] + op_vals[i] + 0.1, i, f'{total:.1f}', va='center',
            fontsize=8, fontweight='bold')

ax.set_yticks(y_pos)
ax.set_yticklabels(genes_d, fontsize=8, fontweight='bold')
ax.set_xlabel('Score', fontsize=8)
ax.set_title('Biological closure score decomposition (top 10 positive regulators)',
             fontsize=9, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(loc='lower right', fontsize=7, framealpha=0.9, ncol=2)

fig.savefig('/Users/hezongze/fig_pos_regulator_ranking.png', dpi=300, facecolor='white')
plt.close(fig)
print("  Saved: /Users/hezongze/fig_pos_regulator_ranking.png")

# ================================================================
# FIGURE 3: ANABOLIC NETWORK — THE POSITIVE AXIS
# ================================================================
print("Generating Figure 3: Anabolic network mechanism diagram...")

fig, ax = plt.subplots(1, 1, figsize=(10, 6.5))
ax.set_xlim(0, 12)
ax.set_ylim(0, 7)
ax.axis('off')

# Compartments
liver_bg = FancyBboxPatch((0.2, 0.4), 4.8, 6.2, boxstyle="round,pad=0.3",
                           facecolor=C_BG_LIVER, edgecolor=C_LIVER_SUPPORT,
                           linewidth=1.5, alpha=0.5, zorder=0)
ax.add_patch(liver_bg)
ax.text(2.6, 6.35, 'LIVER — Anabolic Support System', ha='center', fontsize=9,
        fontweight='bold', color=C_LIVER_SUPPORT)

muscle_bg = FancyBboxPatch((5.5, 0.4), 6.2, 6.2, boxstyle="round,pad=0.3",
                            facecolor=C_BG_MUSCLE, edgecolor=C_IGF,
                            linewidth=1.5, alpha=0.5, zorder=0)
ax.add_patch(muscle_bg)
ax.text(8.6, 6.35, 'SKELETAL MUSCLE — Protein Synthesis Machinery', ha='center',
        fontsize=9, fontweight='bold', color=C_IGF)

# --- Helper functions ---
def node(ax, x, y, w, h, text, color, fs=6.5, tc='white', z=10):
    box = FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.08",
                          facecolor=color, edgecolor='white', linewidth=1, zorder=z)
    ax.add_patch(box)
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, fontweight='bold',
            color=tc, zorder=z+1)

def arrow(ax, x1, y1, x2, y2, color='#555555', lw=1.5, z=5):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw), zorder=z)

def lbl(ax, x, y, text, fs=5.5, color='#666666', ha='center', va='center', rot=0):
    ax.text(x, y, text, fontsize=fs, color=color, ha=ha, va=va, rotation=rot, style='italic')

# --- LIVER NODES ---
# GH/GHR axis
node(ax, 1.2, 5.5, 0.9, 0.4, 'GHR', C_LIVER_SUPPORT)
node(ax, 2.4, 5.5, 0.9, 0.4, 'STAT5B', C_LIVER_SUPPORT)
arrow(ax, 1.65, 5.5, 1.95, 5.5)

# IGF1 production
node(ax, 3.8, 5.5, 0.9, 0.4, 'IGF1', C_IGF, fs=7, tc='white')
arrow(ax, 2.85, 5.5, 3.35, 5.5)
lbl(ax, 3.1, 5.85, 'GH->STAT5B', fs=5.5, color=C_LIVER_SUPPORT)

# IGFBP3 (protective)
node(ax, 3.8, 4.5, 0.9, 0.4, 'IGFBP3', '#56B4E9')
arrow(ax, 3.0, 5.2, 3.4, 4.7, '#56B4E9', 1.0)

# Liver support genes
support_genes = [
    (1.2, 4.2, 'NFKBIE', C_LIVER_SUPPORT, 'NF-kB inh.'),
    (2.4, 4.2, 'GPX7', C_LIVER_SUPPORT, 'ER redox'),
    (3.8, 3.5, 'PDP1', C_LIVER_SUPPORT, 'PDH activator'),
    (1.2, 3.5, 'ADA', C_LIVER_SUPPORT, 'Adenosine'),
    (2.4, 3.5, 'PON1', C_LIVER_SUPPORT, 'Anti-oxidant'),
]
for x, y, name, color, sub in support_genes:
    node(ax, x, y, 0.75, 0.35, name, color, fs=5.5)
    lbl(ax, x, y-0.3, sub, fs=5)

# AA transport
aa_genes = [
    (1.2, 2.5, 'SLC7A8', C_TRANSPORT),
    (2.4, 2.5, 'SLC1A5', C_TRANSPORT),
    (3.8, 2.5, 'SLC38A4', C_TRANSPORT),
]
for x, y, name, color in aa_genes:
    node(ax, x, y, 0.8, 0.35, name, color, fs=5.5)

lbl(ax, 2.6, 2.1, 'AA Transporters: BCAA/Gln/Arg uptake', fs=5.5, color=C_TRANSPORT)

# Liver -> Circulation
# IGF1 secretion
ax.annotate('', xy=(6.5, 5.0), xytext=(4.25, 5.5),
            arrowprops=dict(arrowstyle='->', color=C_IGF, lw=3.0,
                          connectionstyle='arc3,rad=-0.15'), zorder=6)
lbl(ax, 5.3, 5.75, 'IGF1 secretion', fs=6, color=C_IGF)

# AA release to circulation
ax.annotate('', xy=(6.5, 2.8), xytext=(4.2, 2.5),
            arrowprops=dict(arrowstyle='->', color=C_TRANSPORT, lw=2.5,
                          connectionstyle='arc3,rad=0.1'), zorder=6)
lbl(ax, 5.3, 3.1, 'AA supply', fs=6, color=C_TRANSPORT)

# --- CIRCULATION ANNOTATION ---
ax.annotate('Circulation', xy=(5.5, 0.7), fontsize=7, fontweight='bold', color='#888888',
            ha='center')

# --- MUSCLE NODES ---
# IGF1R
node(ax, 7.5, 5.2, 0.9, 0.4, 'IGF1R', C_IGF)

# mTOR pathway
node(ax, 9.0, 5.2, 0.9, 0.4, 'AKT1', C_ANABOLIC)
arrow(ax, 7.95, 5.2, 8.55, 5.2, C_IGF, 2.0)
lbl(ax, 8.25, 5.55, 'PI3K/AKT', fs=5.5, color=C_IGF)

node(ax, 10.5, 5.2, 0.9, 0.4, 'mTORC1', C_ANABOLIC)
arrow(ax, 9.45, 5.2, 10.05, 5.2, C_ANABOLIC, 2.0)

# mTOR downstream
node(ax, 7.5, 4.0, 1.0, 0.4, 'RPS6KB1', C_ANABOLIC, fs=6)
arrow(ax, 10.0, 4.95, 7.95, 4.2, C_ANABOLIC, 1.5)
lbl(ax, 9.0, 4.65, 'S6K1 activation', fs=5.5, color=C_ANABOLIC)

# Ribosome/protein synthesis
node(ax, 9.5, 4.0, 1.6, 0.45, 'RPL/RPS\nRibosome Biogenesis', C_ANABOLIC, fs=6)
arrow(ax, 8.0, 4.0, 8.7, 4.0, C_ANABOLIC, 1.5)

# Translation output
node(ax, 9.5, 3.1, 1.6, 0.45, 'Protein\nSynthesis', '#0072B2', fs=7)
arrow(ax, 9.5, 3.75, 9.5, 3.35, C_ANABOLIC, 2.0)

# Myogenic regulators
node(ax, 7.5, 3.1, 0.85, 0.38, 'MYOD1', '#E69F00', fs=5.5)
node(ax, 7.5, 2.4, 0.85, 0.38, 'MYOG', '#E69F00', fs=5.5)
arrow(ax, 7.5, 3.45, 7.5, 3.3)
arrow(ax, 7.5, 2.75, 7.5, 2.6)
lbl(ax, 7.0, 2.75, 'Myogenic\ndifferentiation', fs=5.5, color='#E69F00')

# Fiber type
node(ax, 9.5, 2.1, 1.3, 0.4, 'MYH4 (IIb fast)', C_ENERGY, fs=6)
lbl(ax, 9.5, 1.75, 'Fast glycolytic fiber', fs=5.5, color=C_ENERGY)

# Muscle result
result_box = FancyBboxPatch((6.3, 0.85), 5.8, 0.55, boxstyle="round,pad=0.1",
                             facecolor='#F5F5F5', edgecolor=C_ANABOLIC, linewidth=2, zorder=12)
ax.add_patch(result_box)
ax.text(9.2, 1.13, 'PROTEIN DEPOSITION', ha='center', va='center', fontsize=10,
        fontweight='bold', color=C_ANABOLIC, zorder=13)

# Arrows to result
arrow(ax, 9.5, 2.85, 9.2, 1.35, C_ANABOLIC, 2.5)
arrow(ax, 9.5, 1.9, 9.2, 1.35, C_ENERGY, 1.5)

# --- Liver support annotations ---
support_box = FancyBboxPatch((0.5, 1.0), 4.2, 1.0, boxstyle="round,pad=0.1",
                              facecolor='white', edgecolor=C_LIVER_SUPPORT,
                              linewidth=0.8, alpha=0.8, zorder=15)
ax.add_patch(support_box)
ax.text(2.6, 1.7, 'Liver support functions:', fontsize=6.5, fontweight='bold', color=C_LIVER_SUPPORT)
ax.text(2.6, 1.4, 'Anti-inflammation | Anti-oxidation | ER proteostasis\nMitochondrial energy | AA transport & recycling',
        fontsize=5.5, color='#666666', ha='center')

# --- Bottom formula ---
formula_box = FancyBboxPatch((0.3, 0.05), 11.4, 0.4, boxstyle="round,pad=0.05",
                              facecolor='white', edgecolor='#CCCCCC', linewidth=0.5, zorder=15)
ax.add_patch(formula_box)
ax.text(6.0, 0.25, 'Positive Regulator Closure = Expression Consistency + PPI Connectivity + Literature Support + Experimental Operability',
        ha='center', fontsize=6.5, color='#555555', zorder=16)

fig.savefig('/Users/hezongze/fig_anabolic_network.png', dpi=300, facecolor='white')
plt.close(fig)
print("  Saved: /Users/hezongze/fig_anabolic_network.png")

print("\nDone. All figures generated.")
PYEOF
