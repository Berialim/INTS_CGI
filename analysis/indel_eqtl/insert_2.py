import pandas as pd
import pyarrow.parquet as pq
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # 使 PDF 中字体可编辑
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"
# --- 0. 设定前缀变量 ---
prefix = 'u1tocgi'
# 'u1tocgi'  # 您可以根据需要修改此变量
# ‘1stexon'
# --- 1. 配置与数据提取 ---
bed_path = '/Users/berial/ljy/lailab2024/features/u1_to_cgi_boundaries_both_strands_filtered_sorted.bed6'
# '/Users/berial/ljy/lailab2024/features/1stexon.bed' 
# '/Users/berial/ljy/lailab2024/features/u1_to_cgi_boundaries_both_strands_filtered_sorted.bed6'

data_dir = './' 
output_file = f'{prefix}_GTEx_v11_Ins_Del_Summary_Filtered.csv' # 加前缀

bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
regions = pd.read_csv(bed_path, sep='\t', names=bed_cols)
regions['chrom'] = regions['chrom'].astype(str)

def process_tissue(file_path, regions_df, data_dir):
    tissue_name = os.path.basename(file_path).split('.')[0]
    print(f"[{prefix}] 正在处理组织: {tissue_name}...") # 加前缀
    columns = ['variant_id', 'phenotype_id', 'slope', 'pval_nominal', 'af']
    
    try:
        table = pq.read_table(file_path, columns=columns)
        df = table.to_pandas().rename(columns={'phenotype_id': 'gene_id'})
        
        egenes_path = os.path.join(data_dir, f"{tissue_name}.v11.eGenes.txt.gz")
        if not os.path.exists(egenes_path):
            return pd.DataFrame()
            
        egenes_df = pd.read_csv(egenes_path, sep='\t', usecols=['gene_id', 'biotype', 'gene_start', 'gene_end'])
        egenes_pc = egenes_df[egenes_df['biotype'] == 'protein_coding']
        
        df = pd.merge(df, egenes_pc, on='gene_id', how='inner')
        if df.empty: return pd.DataFrame()
        
        var_info = df['variant_id'].str.split('_', expand=True)
        df['chrom'] = var_info[0]
        df['pos'] = var_info[1].astype(int)
        df['ref'] = var_info[2]
        df['alt'] = var_info[3]
        
        df = df[(df['pos'] >= df['gene_start']) & (df['pos'] <= df['gene_end'])]

        is_ins = df['alt'].str.len() > df['ref'].str.len() + 5
        is_del = df['ref'].str.len() > df['alt'].str.len() + 5
        
        df_ins = df[is_ins].copy()
        df_ins['sv_type'] = 'Insertion'
        df_del = df[is_del].copy()
        df_del['sv_type'] = 'Deletion'
        
        df_combined = pd.concat([df_ins, df_del])
        
        matched = pd.merge(df_combined, regions_df, on='chrom')
        
        mask_pos_strand = (matched['strand'] == '+') & (matched['pos'] > matched['start'] + 5) & (matched['pos'] <= matched['end'])
        mask_neg_strand = (matched['strand'] == '-') & (matched['pos'] >= matched['start']) & (matched['pos'] < matched['end'] - 5)
        
        final_hits = matched[mask_pos_strand | mask_neg_strand].copy()
        final_hits['tissue'] = tissue_name
        
        return final_hits
        
    except Exception as e:
        print(f"[{prefix}] 处理 {tissue_name} 失败: {e}")
        return pd.DataFrame()

# --- 2. 执行处理与去重 ---
all_files = glob.glob(os.path.join(data_dir, "*.signif_pairs.parquet"))
results = [process_tissue(f, regions, data_dir) for f in all_files]
df_raw = pd.concat(results, ignore_index=True)

if not df_raw.empty:
    df_main = df_raw.sort_values(by=['variant_id', 'gene_id', 'slope'], ascending=[True, True, False])
    
    print(f"[{prefix}] 原始记录数: {len(df_raw)}")
    print(f"[{prefix}] 去重后记录数: {len(df_main)}")
    
    df_main.to_csv(output_file, index=False)
else:
    print(f"[{prefix}] 最终未发现符合条件的记录。")

# --- 3. 数据分类与预处理 ---
COLOR_INS    = '#2980b9'
COLOR_DEL    = '#e74c3c'
MIN_N_BUBBLE = 5
is_ins   = (df_main['sv_type'] == 'Insertion').values
seq_arr  = np.where(is_ins, df_main['alt'].str.upper(), df_main['ref'].str.upper())
seq_len  = np.array([len(s) for s in seq_arr])
cg_cnt   = np.array([s.count('CG') for s in seq_arr])
a_cnt    = np.array([s.count('A')  for s in seq_arr])
t_cnt    = np.array([s.count('T')  for s in seq_arr])

safe_len = np.where(seq_len > 0, seq_len, 1)
cpg_pct  = (cg_cnt * 2) / safe_len * 100
at_pct   = (a_cnt + t_cnt) / safe_len * 100

is_cgi     = cpg_pct >= 20.0
is_at_rich = (~is_cgi) & (at_pct >= 65.0)
sv_label   = np.where(is_ins, 'Insertion', 'Deletion')
sv4_label  = np.where(is_ins, 'Ins', 'Del')

df_main['Group']     = np.where(is_cgi,      sv_label + '_CGI',
                                             sv_label + '_Others')
df_main['Group4']    = np.where(is_cgi, sv4_label + '_CGI', sv4_label + '_Others')
df_main['CGI_Group'] = np.where(is_cgi, sv_label + '_CGI', None)
df_final    = df_main[df_main['CGI_Group'].notna()].copy()
df_plot_all = df_main[df_main['Group4'].notna()].copy()
print(f"\n[{prefix}] CGI子集: {len(df_final):,} 行 | 四分组: {len(df_plot_all):,} 行")

# --- 4. 图A：全局箱线图 ---
sns.set_style("white")
target_order = ['Insertion_CGI', 'Deletion_CGI']
palette      = {"Insertion_CGI": COLOR_INS, "Deletion_CGI": COLOR_DEL}

plt.figure(figsize=(5, 6))
ax = sns.boxplot(data=df_final, x='CGI_Group', y='slope',
                 order=target_order, palette=palette,
                 showfliers=False, width=0.8, linewidth=2, notch=True)
g1 = df_final[df_final['CGI_Group'] == 'Insertion_CGI']['slope']
g2 = df_final[df_final['CGI_Group'] == 'Deletion_CGI']['slope']
if len(g1) >= 3 and len(g2) >= 3:
    _, p_val = stats.mannwhitneyu(g1, g2)
    plt.text(0.5, df_final['slope'].quantile(0.95) * 0.8,
             f'p = {p_val:.2e}', ha='center', va='bottom',
             fontsize=12, fontweight='bold')
plt.axhline(0, color='black', linestyle='--', alpha=0.4)
plt.title(f'{prefix}: CGI Insertion vs Deletion\n(variant inside gene body)', fontsize=13) # 标题加前缀
plt.ylabel('Effect Size (Slope)')
plt.tight_layout()
plt.savefig(f'{prefix}_Global_CGI_InGene_Ins_vs_Del.pdf') # 文件名加前缀
plt.savefig(f'{prefix}_Global_CGI_InGene_Ins_vs_Del.png', dpi=300)
plt.show()

# --- 5. 气泡图：组织分类汇总 ---
TISSUE_CATEGORY = {
    'Brain_Amygdala': 'Brain', 'Brain_Anterior_cingulate_cortex_BA24': 'Brain',
    'Brain_Caudate_basal_ganglia': 'Brain', 'Brain_Cerebellar_Hemisphere': 'Brain',
    'Brain_Cerebellum': 'Brain', 'Brain_Cortex': 'Brain',
    'Brain_Frontal_Cortex_BA9': 'Brain', 'Brain_Hippocampus': 'Brain',
    'Brain_Hypothalamus': 'Brain', 'Brain_Nucleus_accumbens_basal_ganglia': 'Brain',
    'Brain_Putamen_basal_ganglia': 'Brain', 'Brain_Spinal_cord_cervical_c-1': 'Brain',
    'Brain_Substantia_nigra': 'Brain', 'Heart_Atrial_Appendage': 'Heart',
    'Heart_Left_Ventricle': 'Heart', 'Artery_Aorta': 'Artery',
    'Artery_Coronary': 'Artery', 'Artery_Tibial': 'Artery',
    'Adipose_Subcutaneous': 'Adipose', 'Adipose_Visceral_Omentum': 'Adipose',
    'Skin_Not_Sun_Exposed_Suprapubic': 'Skin', 'Skin_Sun_Exposed_Lower_leg': 'Skin',
    'Colon_Sigmoid': 'Colon', 'Colon_Transverse': 'Colon',
    'Esophagus_Gastroesophageal_Junction': 'Esophagus', 'Esophagus_Mucosa': 'Esophagus',
    'Esophagus_Muscularis': 'Esophagus', 'Ovary': 'Reproductive',
    'Uterus': 'Reproductive', 'Vagina': 'Reproductive',
    'Prostate': 'Reproductive', 'Testis': 'Reproductive',
    'Adrenal_Gland': 'Adrenal Gland', 'Bladder': 'Bladder',
    'Breast_Mammary_Tissue': 'Breast', 'Kidney_Cortex': 'Kidney',
    'Liver': 'Liver', 'Lung': 'Lung', 'Minor_Salivary_Gland': 'Salivary Gland',
    'Muscle_Skeletal': 'Muscle', 'Pancreas': 'Pancreas',
    'Pituitary': 'Pituitary', 'Small_Intestine_Terminal_Ileum': 'Small Intestine',
    'Spleen': 'Spleen', 'Stomach': 'Stomach', 'Thyroid': 'Thyroid',
    'Whole_Blood': 'Whole Blood',
}

df_bubble = df_final.copy()
df_bubble['category'] = df_bubble['tissue'].map(TISSUE_CATEGORY)
df_bubble = df_bubble.dropna(subset=['category'])

agg = (df_bubble.groupby(['category', 'CGI_Group'])['slope']
       .agg(n='count', median='median', mean='mean', sem='sem')
       .reset_index())

n_pivot  = agg.pivot_table(index='category', columns='CGI_Group', values='n', fill_value=0)
all_cats = n_pivot.index
ins_col  = n_pivot.get('Insertion_CGI', pd.Series(0, index=all_cats)).reindex(all_cats, fill_value=0)
del_col  = n_pivot.get('Deletion_CGI',  pd.Series(0, index=all_cats)).reindex(all_cats, fill_value=0)
valid_cats = set(all_cats[(ins_col >= MIN_N_BUBBLE) & (del_col >= MIN_N_BUBBLE)])

agg_filtered = agg[agg['category'].isin(valid_cats)].copy()
total_n      = agg_filtered.groupby('category')['n'].sum().sort_values(ascending=False)
top_cats     = list(total_n.head(15).index)
agg_top      = agg_filtered[agg_filtered['category'].isin(top_cats)].copy()
present_cats = top_cats
y_pos        = {cat: i for i, cat in enumerate(reversed(present_cats))}

print(f"\n[{prefix}] 气泡图大分类（各组 n ≥ {MIN_N_BUBBLE}）：")
for cat in present_cats:
    ni = int(ins_col.get(cat, 0))
    nd = int(del_col.get(cat, 0))
    print(f"  {cat:<20} Ins={ni:>3}  Del={nd:>3}  total={ni+nd:>3}")

pval_dict = {}
for cat in present_cats:
    sub = df_bubble[df_bubble['category'] == cat]
    g1  = sub[sub['CGI_Group'] == 'Insertion_CGI']['slope']
    g2  = sub[sub['CGI_Group'] == 'Deletion_CGI']['slope']
    if len(g1) >= 3 and len(g2) >= 3:
        _, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
        pval_dict[cat] = p
    else:
        pval_dict[cat] = np.nan
agg_top['pval'] = agg_top['category'].map(pval_dict)

# 气泡图绘制
fig, axes = plt.subplots(1, 2, figsize=(16, max(7, len(present_cats) * 0.65 + 2)), gridspec_kw={'wspace': 0.06})
for ax_idx, metric in enumerate(['median', 'mean']):
    ax = axes[ax_idx]
    sig_annotations = []
    for i in range(len(present_cats)):
        if i % 2 == 0: ax.axhspan(i - 0.5, i + 0.5, color='#f4f6f8', zorder=0)
    ax.axvline(0, color='#7f8c8d', lw=1, linestyle='--', alpha=0.7, zorder=1)

    for cat in present_cats:
        y = y_pos[cat]
        for grp, color, marker in [('Insertion_CGI', COLOR_INS, 'o'), ('Deletion_CGI', COLOR_DEL, 's')]:
            row = agg_top[(agg_top['category'] == cat) & (agg_top['CGI_Group'] == grp)]
            if not row.empty:
                xv = row[metric].values[0]
                sz = np.clip(np.log10(row['n'].values[0] + 1) * 150, 30, 700)
                ax.scatter(xv, y, s=sz, c=color, marker=marker, edgecolors='white', lw=0.8, alpha=0.88, zorder=4)
                ax.errorbar(xv, y, xerr=row['sem'].values[0], fmt='none', ecolor=color, elinewidth=1.5, capsize=4, alpha=0.75, zorder=3)
        
        p = pval_dict.get(cat, np.nan)
        sig = ('***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns')
        sig_annotations.append((y, sig))

    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(list(reversed(present_cats)), fontsize=10.5)
    if ax_idx == 1: ax.set_yticklabels([])
    ax.set_xlabel(f'Slope {metric.capitalize()} (±SEM)')
    ax.set_title(f'{"Median" if metric == "median" else "Mean"} Effect Size')
    ax.set_ylim(-0.8, len(present_cats) - 0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    xlim = ax.get_xlim()
    x_sig = xlim[1] + (xlim[1] - xlim[0]) * 0.04
    for y_s, sig in sig_annotations:
        ax.text(x_sig, y_s, sig, va='center', ha='left', fontsize=9, color='#c0392b' if sig != 'ns' else '#95a5a6', fontweight='bold' if sig != 'ns' else 'normal')

fig.suptitle(f'[{prefix}] CGI Insertion (blue ●) vs Deletion (red ■): Slope by Tissue Category\nGTEx v11 | >10bp Indels', fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f'{prefix}_Bubble_Plot_InGene_Category.pdf', bbox_inches='tight', dpi=300) # 文件名加前缀
plt.savefig(f'{prefix}_Bubble_Plot_InGene_Category.png', bbox_inches='tight', dpi=300)
plt.show()

print(f"✅ [{prefix}] 气泡图已保存（{len(present_cats)} 个大类）")
agg_top.to_csv(f'{prefix}_Bubble_InGene_Category_Summary.csv', index=False) # 文件名加前缀
print(f"✅ [{prefix}] 汇总表: {prefix}_Bubble_InGene_Category_Summary.csv")

# %%

# --- 修改后的横置气泡图绘制代码 ---

# 计算绘图所需的宽度（根据组织数量动态调整）
plot_width = max(10, len(present_cats) * 0.8)
fig, axes = plt.subplots(
    2, 1, # 改为 2 行 1 列
    figsize=(plot_width, 12), 
    gridspec_kw={'hspace': 0.3} # 调整行间距
)

# 重新定义 x 轴位置
x_indices = np.arange(len(present_cats))
x_pos = {cat: i for i, cat in enumerate(present_cats)}

for ax_idx, metric in enumerate(['median', 'mean']):
    ax = axes[ax_idx]
    sig_annotations = []

    # 背景条纹 (改为垂直方向)
    for i in range(len(present_cats)):
        if i % 2 == 0:
            ax.axvspan(i - 0.5, i + 0.5, color='#f4f6f8', zorder=0)
    
    # 水平基准线 (Slope=0)
    ax.axhline(0, color='#7f8c8d', lw=1, linestyle='--', alpha=0.7, zorder=1)

    for cat in present_cats:
        x = x_pos[cat]
        # 绘制 Insertion
        ins_row = agg_top[(agg_top['category'] == cat) & (agg_top['CGI_Group'] == 'Insertion_CGI')]
        if not ins_row.empty:
            yv = ins_row[metric].values[0]
            sz = np.clip(np.log10(ins_row['n'].values[0] + 1) * 150, 30, 700)
            ax.scatter(x, yv, s=sz, c=COLOR_INS, marker='o', edgecolors='white', lw=0.8, alpha=0.88, zorder=4)
            ax.errorbar(x, yv, yerr=ins_row['sem'].values[0], fmt='none', ecolor=COLOR_INS, elinewidth=1.5, capsize=4, alpha=0.75, zorder=3)

        # 绘制 Deletion
        del_row = agg_top[(agg_top['category'] == cat) & (agg_top['CGI_Group'] == 'Deletion_CGI')]
        if not del_row.empty:
            yv = del_row[metric].values[0]
            sz = np.clip(np.log10(del_row['n'].values[0] + 1) * 150, 30, 700)
            ax.scatter(x, yv, s=sz, c=COLOR_DEL, marker='s', edgecolors='white', lw=0.8, alpha=0.88, zorder=4)
            ax.errorbar(x, yv, yerr=del_row['sem'].values[0], fmt='none', ecolor=COLOR_DEL, elinewidth=1.5, capsize=4, alpha=0.75, zorder=3)

        # 显著性标注
        p = pval_dict.get(cat, np.nan)
        sig = ('***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns')
        sig_annotations.append((x, sig))

    # 设置轴标签
    ax.set_xticks(x_indices)
    if ax_idx == 1: # 仅在底部子图显示组织名称
        ax.set_xticklabels(present_cats, fontsize=10, rotation=45, ha='right')
    else:
        ax.set_xticklabels([])
        
    ax.set_ylabel(f'Slope {metric.capitalize()} (±SEM)', fontsize=11)
    ax.set_title(f'{"Median" if metric == "median" else "Mean"} Effect Size', fontsize=12, fontweight='bold')
    
    # 设置显著性文本位置（放在 y 轴顶部）
    ylim = ax.get_ylim()
    y_sig = ylim[1] + (ylim[1] - ylim[0]) * 0.05
    for x_s, sig in sig_annotations:
        ax.text(x_s, y_sig, sig, va='bottom', ha='center', fontsize=9, 
                color='#c0392b' if sig != 'ns' else '#95a5a6',
                fontweight='bold' if sig != 'ns' else 'normal')
    ax.set_ylim(ylim[0], y_sig + (ylim[1] - ylim[0]) * 0.1)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# 图例（挪到右侧或下方）
legend_elements = [
    plt.scatter([], [], s=100, c=COLOR_INS, marker='o', edgecolors='white', label='Insertion CGI'),
    plt.scatter([], [], s=100, c=COLOR_DEL, marker='s', edgecolors='white', label='Deletion CGI'),
    mpatches.Patch(color='white', label=''), # 占位
    plt.scatter([], [], s=np.log10(51)*150, c='gray', marker='o', alpha=0.6, label='n = 50'),
]
fig.legend(handles=legend_elements, loc='center right', bbox_to_anchor=(1.1, 0.5), frameon=True)

fig.suptitle(f'[{prefix}] CGI Ins vs Del: Slope by Tissue Category (Horizontal View)', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 0.9, 1]) # 为右侧图例留出空间

# 保存
plt.savefig(f'{prefix}_Bubble_Plot_Horizontal.pdf', bbox_inches='tight', dpi=300)
plt.savefig(f'{prefix}_Bubble_Plot_Horizontal.png', bbox_inches='tight', dpi=300)
plt.show()