import pandas as pd
import pyarrow.parquet as pq
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# --- 1. 配置与数据提取 ---
bed_path =  '/Users/berial/ljy/lailab2024/features/u1_to_cgi_boundaries_both_strands_filtered_sorted.bed6'
# '/Users/berial/ljy/lailab2024/features/1stexon.bed' # '/Users/berial/ljy/lailab2024/features/u1_to_cgi_boundaries_both_strands_filtered_sorted.bed6'
data_dir = './' 
output_file = 'GTEx_v11_Ins_Del_Summary.csv'

bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
regions = pd.read_csv(bed_path, sep='\t', names=bed_cols)
regions['chrom'] = regions['chrom'].astype(str)

def process_tissue(file_path, regions_df):
    tissue_name = os.path.basename(file_path).split('.')[0]
    print(f"正在处理组织: {tissue_name}...")
    columns = ['variant_id', 'phenotype_id', 'slope', 'pval_nominal', 'af']
    
    try:
        table = pq.read_table(file_path, columns=columns)
        df = table.to_pandas().rename(columns={'phenotype_id': 'gene_id'})
        
        var_info = df['variant_id'].str.split('_', expand=True)
        if var_info.shape[1] < 4: return pd.DataFrame()

        # 修改核心逻辑：提取 >10bp 的插入 OR 缺失
        is_ins = var_info[3].str.len() > var_info[2].str.len() + 10
        is_del = var_info[2].str.len() > var_info[3].str.len() + 10
        
        df_ins = df[is_ins].copy()
        df_ins['sv_type'] = 'Insertion'
        
        df_del = df[is_del].copy()
        df_del['sv_type'] = 'Deletion'
        
        df_combined = pd.concat([df_ins, df_del])
        
        df_combined['chrom'] = var_info.loc[df_combined.index, 0]
        df_combined['pos'] = var_info.loc[df_combined.index, 1].astype(int)
        
        matched = pd.merge(df_combined, regions_df, on='chrom')
        final_hits = matched[(matched['pos'] >= matched['start']) & (matched['pos'] <= matched['end'])].copy()
        final_hits['tissue'] = tissue_name
        return final_hits
    except Exception as e:
        print(f"失败: {e}")
        return pd.DataFrame()

# 循环提取
all_files = glob.glob(os.path.join(data_dir, "*.signif_pairs.parquet"))
results = [process_tissue(f, regions) for f in all_files]
df_main = pd.concat(results, ignore_index=True)

# --- 2. 增强分类逻辑 (Ins vs Del & CGI) ---
def classify_sv(row):
    try:
        parts = row['variant_id'].split('_')
        ref, alt = parts[2].upper(), parts[3].upper()
        
        # 如果是插入，分析插入的序列 (Alt)；如果是缺失，分析被删掉的序列 (Ref)
        seq_to_analyze = alt if row['sv_type'] == 'Insertion' else ref
        seq_len = len(seq_to_analyze)
        
        cg_count = seq_to_analyze.count('CG')
        cpg_pct = (cg_count * 2) / seq_len * 100
        at_pct = ((seq_to_analyze.count('A') + seq_to_analyze.count('T')) / seq_len) * 100
        
        if cpg_pct >= 30.0:
            return f"{row['sv_type']}_CGI"
        elif at_pct >= 65.0:
            return f"{row['sv_type']}_AT-rich"
        else:
            return f"{row['sv_type']}_Others"
    except:
        return None

df_main['Group'] = df_main.apply(classify_sv, axis=1)
df_main.to_csv(output_file, index=False)





# %%

def filter_cgi_only(row):
    try:
        parts = row['variant_id'].split('_')
        # 插入看 Alt，缺失看 Ref
        seq = parts[3].upper() if row['sv_type'] == 'Insertion' else parts[2].upper()
        cg_count = seq.count('CG')
        cpg_pct = (cg_count * 2) / len(seq) * 100
        # 只要 CpG% >= 20% 就标记为对应的 CGI 类
        if cpg_pct >= 20.0:
            return f"{row['sv_type']}_CGI"
        return None
    except:
        return None

df_main['CGI_Group'] = df_main.apply(filter_cgi_only, axis=1)
df_final = df_main.dropna(subset=['CGI_Group']).copy()
# --- 3. 统计展示与绘图 ---
sns.set_style("white")
target_order = ['Insertion_CGI', 'Deletion_CGI']
palette = {"Insertion_CGI": "#e74c3c", "Deletion_CGI": "#3498db"}

# A. 全局对比图 (Combined Tissues)
plt.figure(figsize=(5, 6))
ax = sns.boxplot(data=df_final, x='CGI_Group', y='slope', order=target_order, 
                 palette=palette, showfliers=False, width=0.8, 
                 linewidth=2, notch=True) # 使用 notch 展示中位数置信区间

# 计算全局 P 值 (Mann-Whitney U)
g1 = df_final[df_final['CGI_Group'] == 'Insertion_CGI']['slope']
g2 = df_final[df_final['CGI_Group'] == 'Deletion_CGI']['slope']
stat, p_val = stats.mannwhitneyu(g1, g2)

# 标注 P 值
y_max = df_final['slope'].max() * 0.8
plt.text(0.5, y_max, f'p = {p_val:.2e}', ha='center', va='bottom', 
         fontsize=12, fontweight='bold', color='black')
plt.axhline(0, color='black', linestyle='--', alpha=0.4)
plt.title('CGI Insertion vs Deletion (All Tissues)', fontsize=14)
plt.ylabel('Effect Size (Slope)')
plt.savefig('Global_CGI_Ins_vs_Del.pdf')
plt.show()

# B. 分组织展示 (FacetPlot)
# 选取样本量较多的前 12 个组织
top_tissues = df_final['tissue'].value_counts()
display_tissues = top_tissues[top_tissues >= 10].index[:12]
df_facet = df_final[df_final['tissue'].isin(display_tissues)]

g = sns.FacetGrid(df_facet, col="tissue", col_wrap=4, height=4, aspect=0.8, sharey=False)
g.map_dataframe(sns.boxplot, x='CGI_Group', y='slope', order=target_order, 
                palette=palette, showfliers=False, width=0.6)

# 为每个子图添加 P 值
for ax in g.axes.flat:
    t_name = ax.get_title().split('=')[-1].strip()
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)
    t_data = df_facet[df_facet['tissue'] == t_name]
    ins_s = t_data[t_data['CGI_Group'] == 'Insertion_CGI']['slope']
    del_s = t_data[t_data['CGI_Group'] == 'Deletion_CGI']['slope']
    
    if len(ins_s) > 3 and len(del_s) > 3:
        _, p = stats.mannwhitneyu(ins_s, del_s)
        ax.set_title(f"{t_name}\n(p={p:.1e})", fontsize=10)
    else:
        ax.set_title(t_name, fontsize=10)

g.set_axis_labels("", "Slope")
plt.subplots_adjust(top=0.9, hspace=0.4)
g.fig.suptitle('Tissue-specific CGI Insertion vs Deletion (No Scatter)', fontsize=16)
plt.savefig('Tissue_CGI_Ins_vs_Del_Box.png', dpi=300)
plt.show()


# %%

# ================================================================
# 以下代码直接追加到原脚本末尾（B图之后）
# 依赖变量：df_final（已有 CGI_Group 和 tissue 列）
# ================================================================
import matplotlib.patches as mpatches

COLOR_INS = '#2980b9'   # 插入=蓝
COLOR_DEL = '#e74c3c'   # 删除=红

# ================================================================
# --- C. 气泡图：按大分类汇总，展示 Top 15（Ins_CGI+Del_CGI 总数排名）
#         过滤：每个大类 Ins_CGI>=5 且 Del_CGI>=5
#         ● 蓝圆=Insertion CGI   ■ 红方=Deletion CGI
#         X轴=slope median（左）/ mean（右）±SEM
#         气泡大小 ∝ log10(n)
# ================================================================

# ── 大分类映射 ────────────────────────────────────────────────
TISSUE_CATEGORY = {
    'Brain_Amygdala':                        'Brain',
    'Brain_Anterior_cingulate_cortex_BA24':  'Brain',
    'Brain_Caudate_basal_ganglia':           'Brain',
    'Brain_Cerebellar_Hemisphere':           'Brain',
    'Brain_Cerebellum':                      'Brain',
    'Brain_Cortex':                          'Brain',
    'Brain_Frontal_Cortex_BA9':              'Brain',
    'Brain_Hippocampus':                     'Brain',
    'Brain_Hypothalamus':                    'Brain',
    'Brain_Nucleus_accumbens_basal_ganglia': 'Brain',
    'Brain_Putamen_basal_ganglia':           'Brain',
    'Brain_Spinal_cord_cervical_c-1':        'Brain',
    'Brain_Substantia_nigra':                'Brain',
    'Heart_Atrial_Appendage':                'Heart',
    'Heart_Left_Ventricle':                  'Heart',
    'Artery_Aorta':                          'Artery',
    'Artery_Coronary':                       'Artery',
    'Artery_Tibial':                         'Artery',
    'Adipose_Subcutaneous':                  'Adipose',
    'Adipose_Visceral_Omentum':              'Adipose',
    'Skin_Not_Sun_Exposed_Suprapubic':       'Skin',
    'Skin_Sun_Exposed_Lower_leg':            'Skin',
    'Colon_Sigmoid':                         'Colon',
    'Colon_Transverse':                      'Colon',
    'Esophagus_Gastroesophageal_Junction':   'Esophagus',
    'Esophagus_Mucosa':                      'Esophagus',
    'Esophagus_Muscularis':                  'Esophagus',
    'Cells_Cultured_fibroblasts':            'Cells',
    'Cells_EBV-transformed_lymphocytes':     'Cells',
    'Ovary':                                 'Reproductive',
    'Uterus':                                'Reproductive',
    'Vagina':                                'Reproductive',
    'Prostate':                              'Reproductive',
    'Testis':                                'Reproductive',
    'Adrenal_Gland':                         'Adrenal Gland',
    'Bladder':                               'Bladder',
    'Breast_Mammary_Tissue':                 'Breast',
    'Kidney_Cortex':                         'Kidney',
    'Liver':                                 'Liver',
    'Lung':                                  'Lung',
    'Minor_Salivary_Gland':                  'Salivary Gland',
    'Muscle_Skeletal':                       'Muscle',
    'Nerve_Tibial':                          'Nerve',
    'Pancreas':                              'Pancreas',
    'Pituitary':                             'Pituitary',
    'Small_Intestine_Terminal_Ileum':        'Small Intestine',
    'Spleen':                                'Spleen',
    'Stomach':                               'Stomach',
    'Thyroid':                               'Thyroid',
    'Whole_Blood':                           'Whole Blood',
}

MIN_N_BUBBLE = 5   # 每个大类每组最小样本量

# ── Step 1: 映射大分类 ────────────────────────────────────────
df_bubble = df_final.copy()
df_bubble['category'] = df_bubble['tissue'].map(TISSUE_CATEGORY)
df_bubble = df_bubble.dropna(subset=['category'])

# ── Step 2: 按 (category, CGI_Group) 聚合 ────────────────────
agg = (df_bubble.groupby(['category', 'CGI_Group'])['slope']
       .agg(n='count', median='median', mean='mean', sem='sem')
       .reset_index())

# ── Step 3: 过滤 Ins_CGI>=MIN_N 且 Del_CGI>=MIN_N ────────────
n_pivot = agg.pivot_table(index='category', columns='CGI_Group',
                           values='n', fill_value=0)
all_cats = n_pivot.index
ins_col  = (n_pivot['Insertion_CGI'] if 'Insertion_CGI' in n_pivot.columns
            else pd.Series(0, index=all_cats))
del_col  = (n_pivot['Deletion_CGI']  if 'Deletion_CGI'  in n_pivot.columns
            else pd.Series(0, index=all_cats))
ins_col  = ins_col.reindex(all_cats, fill_value=0)
del_col  = del_col.reindex(all_cats, fill_value=0)
valid_cats = set(all_cats[(ins_col >= MIN_N_BUBBLE) & (del_col >= MIN_N_BUBBLE)])

agg_filtered = agg[agg['category'].isin(valid_cats)].copy()

# ── Step 4: 按 Ins_CGI+Del_CGI 总数排名，取 Top 15 ───────────
total_n = (agg_filtered.groupby('category')['n'].sum()
           .sort_values(ascending=False))
top_cats = list(total_n.head(15).index)
agg_top  = agg_filtered[agg_filtered['category'].isin(top_cats)].copy()

# Y轴顺序：从上到下 = 总数最多 → 最少
# matplotlib Y轴从下到上，所以 enumerate(reversed) 让最多的在顶部
present_cats = list(total_n.head(15).index)           # 降序
y_pos = {cat: i for i, cat in enumerate(reversed(present_cats))}

print(f"\n气泡图大分类（各组 n ≥ {MIN_N_BUBBLE}，按总数排名）：")
for cat in present_cats:
    n_ins = int(n_pivot.loc[cat, 'Insertion_CGI']) if cat in n_pivot.index else 0
    n_del = int(n_pivot.loc[cat, 'Deletion_CGI'])  if cat in n_pivot.index else 0
    print(f"  {cat:<20} Ins={n_ins:>3}  Del={n_del:>3}  total={n_ins+n_del:>3}")

# ── Step 5: Mann-Whitney p值 ──────────────────────────────────
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

# ── Step 6: 绘图 ──────────────────────────────────────────────
fig, axes = plt.subplots(
    1, 2,
    figsize=(16, max(7, len(present_cats) * 0.65 + 2)),
    gridspec_kw={'wspace': 0.06}
)

for ax_idx, metric in enumerate(['median', 'mean']):
    ax = axes[ax_idx]
    sig_annotations = []

    # 斑马纹背景
    for i in range(len(present_cats)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color='#f4f6f8', zorder=0)
    ax.axvline(0, color='#7f8c8d', lw=1, linestyle='--', alpha=0.7, zorder=1)

    for cat in present_cats:
        y = y_pos[cat]
        ins_row = agg_top[(agg_top['category'] == cat) &
                          (agg_top['CGI_Group'] == 'Insertion_CGI')]
        del_row = agg_top[(agg_top['category'] == cat) &
                          (agg_top['CGI_Group'] == 'Deletion_CGI')]

        # 插入：蓝色圆形
        if not ins_row.empty:
            xv  = ins_row[metric].values[0]
            n   = ins_row['n'].values[0]
            sem = ins_row['sem'].values[0]
            sz  = np.clip(np.log10(n + 1) * 150, 30, 700)
            ax.scatter(xv, y, s=sz, c=COLOR_INS, marker='o',
                       edgecolors='white', lw=0.8, alpha=0.88, zorder=4)
            ax.errorbar(xv, y, xerr=sem, fmt='none',
                        ecolor=COLOR_INS, elinewidth=1.5,
                        capsize=4, alpha=0.75, zorder=3)

        # 删除：红色方形
        if not del_row.empty:
            xv  = del_row[metric].values[0]
            n   = del_row['n'].values[0]
            sem = del_row['sem'].values[0]
            sz  = np.clip(np.log10(n + 1) * 150, 30, 700)
            ax.scatter(xv, y, s=sz, c=COLOR_DEL, marker='s',
                       edgecolors='white', lw=0.8, alpha=0.88, zorder=4)
            ax.errorbar(xv, y, xerr=sem, fmt='none',
                        ecolor=COLOR_DEL, elinewidth=1.5,
                        capsize=4, alpha=0.75, zorder=3)

        # 连线
        if not ins_row.empty and not del_row.empty:
            ax.plot([ins_row[metric].values[0], del_row[metric].values[0]],
                    [y, y], color='#95a5a6', lw=1.2, alpha=0.5, zorder=2)

        # 显著性标记
        p   = pval_dict.get(cat, np.nan)
        sig = ('***' if p < 0.001 else
               '**'  if p < 0.01  else
               '*'   if p < 0.05  else 'ns')
        sig_annotations.append((y, sig))

    # Y轴
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(list(reversed(present_cats)), fontsize=10.5)
    if ax_idx == 1:
        ax.set_yticklabels([])

    ax.set_xlabel(f'Slope {metric.capitalize()} (±SEM)', fontsize=12)
    ax.set_title(
        f'{"Median" if metric == "median" else "Mean"} Effect Size',
        fontsize=13, fontweight='bold', pad=8
    )
    ax.set_ylim(-0.8, len(present_cats) - 0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x', labelsize=9)

    # 显著性标注（图右侧）
    xlim    = ax.get_xlim()
    x_range = xlim[1] - xlim[0]
    x_sig   = xlim[1] + x_range * 0.04
    for y_s, sig in sig_annotations:
        ax.text(x_sig, y_s, sig,
                va='center', ha='left', fontsize=9,
                color='#c0392b' if sig != 'ns' else '#95a5a6',
                fontweight='bold' if sig != 'ns' else 'normal')
    ax.set_xlim(xlim[0], x_sig + x_range * 0.15)

# ── 图例 ──────────────────────────────────────────────────────
legend_elements = [
    plt.scatter([], [], s=130, c=COLOR_INS, marker='o',
                edgecolors='white', label=f'Insertion CGI  (n ≥ {MIN_N_BUBBLE})'),
    plt.scatter([], [], s=130, c=COLOR_DEL, marker='s',
                edgecolors='white', label=f'Deletion CGI   (n ≥ {MIN_N_BUBBLE})'),
    mpatches.Patch(color='white', label=''),
    plt.scatter([], [], s=np.log10(6)*150,   c='gray', marker='o',
                alpha=0.6, label='n = 5'),
    plt.scatter([], [], s=np.log10(51)*150,  c='gray', marker='o',
                alpha=0.6, label='n = 50'),
    plt.scatter([], [], s=np.log10(201)*150, c='gray', marker='o',
                alpha=0.6, label='n = 200'),
    mpatches.Patch(color='white', label='bubble ∝ log₁₀(n)'),
    mpatches.Patch(color='white', label=''),
    mpatches.Patch(color='#c0392b', label='*** p < 0.001', alpha=0.8),
    mpatches.Patch(color='#c0392b', label='**  p < 0.01',  alpha=0.6),
    mpatches.Patch(color='#c0392b', label='*   p < 0.05',  alpha=0.4),
    mpatches.Patch(color='#95a5a6', label='ns  p ≥ 0.05',  alpha=0.6),
]
fig.legend(
    handles=legend_elements,
    loc='center right', bbox_to_anchor=(1.20, 0.5),
    fontsize=9, frameon=True, framealpha=0.95,
    edgecolor='#bdc3c7', title='Legend', title_fontsize=10
)
fig.suptitle(
    f'CGI Insertion (blue ●) vs Deletion (red ■): Slope by Tissue Category\n'
    f'GTEx v11 | >10bp Indels in BED regions | Top {len(present_cats)} categories '
    f'(Ins_CGI & Del_CGI ≥ {MIN_N_BUBBLE})',
    fontsize=12, fontweight='bold', y=1.02
)
plt.tight_layout()
plt.savefig('Bubble_Plot_Category_Slope.pdf', bbox_inches='tight', dpi=300)
plt.savefig('Bubble_Plot_Category_Slope.png', bbox_inches='tight', dpi=300)
plt.show()
print(f"✅ 气泡图已保存（{len(present_cats)} 个大类）")

# ── 汇总表 ────────────────────────────────────────────────────
agg_top.to_csv('Bubble_Category_Summary.csv', index=False)
print("✅ 汇总表已保存: Bubble_Category_Summary.csv")
print(agg_top[['category', 'CGI_Group', 'n', 'median', 'mean', 'sem', 'pval']]
      .sort_values(['category', 'CGI_Group']).to_string(index=False))


# %%



'''





# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statannotations.Annotator import Annotator

# --- 1. 数据分类与预处理 ---
def classify_four_groups(row):
    try:
        parts = row['variant_id'].split('_')
        # 插入看 Alt (索引3)，缺失看 Ref (索引2)
        seq = parts[3].upper() if row['sv_type'] == 'Insertion' else parts[2].upper()
        
        seq_len = len(seq)
        if seq_len < 2: return None
        
        cg_count = seq.count('CG')
        cpg_pct = (cg_count * 2) / seq_len * 100
        
        # 定义 CGI 阈值 (20%)
        is_cgi = cpg_pct >= 20.0
        
        if row['sv_type'] == 'Insertion':
            return "Ins_CGI" if is_cgi else "Ins_Others"
        else:
            return "Del_CGI" if is_cgi else "Del_Others"
    except:
        return None

# 应用分类逻辑并清洗数据
df_main['Group'] = df_main.apply(classify_four_groups, axis=1)
df_plot_all = df_main.dropna(subset=['Group']).copy()

# --- 2. 绘图参数设置 ---
order = ['Ins_CGI', 'Ins_Others', 'Del_CGI']
palette = {
    "Ins_CGI": "#d63031",    # 深红
    "Ins_Others": "#fab1a0", # 浅橙
    "Del_CGI": "#0984e3",    # 深蓝
    "Del_Others": "#74b9ff"  # 浅蓝
}

# --- 3. 绘制全局对比图 (使用 Annotator) ---
def draw_global_with_annotator(data, order, palette):
    plt.figure(figsize=(9, 7))
    # sns.set_style("ticks")
    
    # 绘制基础箱线图
    ax = sns.boxplot(data=data, x='Group', y='slope', order=order, 
                     palette=palette, showfliers=False, width=0.6, linewidth=1.5,  notch=True)
    
    # 添加基准线
    plt.axhline(0, color='black', linestyle='--', alpha=0.5, lw=1)

    # 定义需要对比的组对
    pairs = [
        ("Ins_CGI", "Del_CGI"),
        ("Ins_CGI", "Ins_Others")
    ]

    # 初始化 Annotator
    annotator = Annotator(ax, pairs, data=data, x='Group', y='slope', order=order)
    
    # 配置统计检验 (Mann-Whitney) 与 标注格式
    annotator.configure(
        test='Mann-Whitney', 
        text_format='full', 
        loc='outside', 
        show_test_name=False
    )
    
    # 执行计算并应用到图中
    annotator.apply_and_annotate()

    # plt.title('Expression Change (Slope) Across Four SV Categories', fontsize=14, pad=20)
    plt.ylabel('Effect Size (Slope)')
    plt.xlabel('')
    sns.despine()
    plt.tight_layout()
    plt.savefig('Four_Groups_Global_Annotated.pdf')
    plt.show()

# 执行绘图
draw_global_with_annotator(df_plot_all, order, palette)

# --- 4. 分组织展示 (FacetPlot 保持原有逻辑) ---
# 选取样本量较多的组织
top_tissues = df_plot_all['tissue'].value_counts()
display_tissues = top_tissues[top_tissues >= 15].index[:12]
df_facet = df_plot_all[df_plot_all['tissue'].isin(display_tissues)]

g = sns.FacetGrid(df_facet, col="tissue", col_wrap=4, height=4, aspect=0.9, sharey=False)
g.map_dataframe(sns.boxplot, x='Group', y='slope', order=order, 
                palette=palette, showfliers=False, width=0.7)

for ax in g.axes.flat:
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)
    ax.set_xticklabels(order, rotation=45)
    t_name = ax.get_title().split('=')[-1].strip()
    t_data = df_facet[df_facet['tissue'] == t_name]
    
    # 对每个组织内部进行简单的显著性标注 (Ins_CGI vs Del_CGI)
    i_cgi = t_data[t_data['Group'] == 'Ins_CGI']['slope']
    d_cgi = t_data[t_data['Group'] == 'Del_CGI']['slope']
    if len(i_cgi) > 3 and len(d_cgi) > 3:
        _, p = stats.mannwhitneyu(i_cgi, d_cgi)
        ax.set_title(f"{t_name}\n(p={p:.1e})", fontsize=10)

g.set_axis_labels("", "Slope")
plt.subplots_adjust(top=0.9, hspace=0.6)
plt.savefig('Four_Groups_Tissue_Facet.png', dpi=300)
plt.show()
# %%

# %%

################# FOR down stream region as control

import pandas as pd
import pyarrow.parquet as pq
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statannotations.Annotator import Annotator

# --- 1. 配置路径 ---
bed_path = '/Users/berial/ljy/lailab2024/features/u1_to_cgi_boundaries_both_strands_filtered_sorted.bed6'
data_dir = './' 
output_file = 'GTEx_v11_Downstream_2kb_Summary.csv'

# --- 2. 加载并偏移 BED 区域 (计算下游 2kb) ---
print("正在计算下游 2kb 区域...")
bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
regions = pd.read_csv(bed_path, sep='\t', names=bed_cols)
regions['chrom'] = regions['chrom'].astype(str)

def shift_to_downstream(row, window=2000):
    """根据正负链计算下游 2000bp 的新坐标"""
    if row['strand'] == '+':
        new_start = row['end']
        new_end = row['end'] + window
    else:
        new_start = row['start'] - window
        new_end = row['start']
    return pd.Series([new_start, new_end])

# 更新坐标为下游 2kb
regions[['start', 'end']] = regions.apply(shift_to_downstream, axis=1)

# --- 3. 定义处理函数 (提取插入与缺失) ---
def process_tissue_downstream(file_path, regions_df):
    tissue_name = os.path.basename(file_path).split('.')[0]
    columns = ['variant_id', 'phenotype_id', 'slope', 'pval_nominal']
    
    try:
        table = pq.read_table(file_path, columns=columns)
        df = table.to_pandas().rename(columns={'phenotype_id': 'gene_id'})
        
        var_info = df['variant_id'].str.split('_', expand=True)
        if var_info.shape[1] < 4: return pd.DataFrame()

        # 提取 >10bp 的变异
        is_ins = var_info[3].str.len() > var_info[2].str.len() + 10
        is_del = var_info[2].str.len() > var_info[3].str.len() + 10
        
        df_ins = df[is_ins].copy()
        df_ins['sv_type'] = 'Insertion'
        df_del = df[is_del].copy()
        df_del['sv_type'] = 'Deletion'
        
        df_combined = pd.concat([df_ins, df_del])
        df_combined['chrom'] = var_info.loc[df_combined.index, 0]
        df_combined['pos'] = var_info.loc[df_combined.index, 1].astype(int)
        
        # 空间合并 (针对下游 2kb 区域)
        matched = pd.merge(df_combined, regions_df, on='chrom')
        final_hits = matched[(matched['pos'] >= matched['start']) & (matched['pos'] <= matched['end'])].copy()
        final_hits['tissue'] = tissue_name
        return final_hits
    except Exception as e:
        print(f"处理 {tissue_name} 失败: {e}")
        return pd.DataFrame()

# 循环提取
all_files = glob.glob(os.path.join(data_dir, "*.signif_pairs.parquet"))
results = [process_tissue_downstream(f, regions) for f in all_files]
df_downstream = pd.concat(results, ignore_index=True)
def classify_four_groups(row):
    try:
        parts = row['variant_id'].split('_')
        seq = parts[3].upper() if row['sv_type'] == 'Insertion' else parts[2].upper()
        seq_len = len(seq)
        if seq_len < 2: return None
        cg_count = seq.count('CG')
        cpg_pct = (cg_count * 2) / seq_len * 100
        is_cgi = cpg_pct >= 20.0
        return f"{'Ins' if row['sv_type'] == 'Insertion' else 'Del'}_{'CGI' if is_cgi else 'Others'}"
    except: return None

df_downstream['Group'] = df_downstream.apply(classify_four_groups, axis=1)
df_plot = df_downstream.dropna(subset=['Group']).copy()
# %%


# 绘图配置
order = ['Ins_CGI', 'Ins_Others', 'Del_CGI']
palette = {"Ins_CGI": "#d63031", "Ins_Others": "#fab1a0", "Del_CGI": "#0984e3", "Del_Others": "#74b9ff"}

plt.figure(figsize=(9, 7))
ax = sns.boxplot(data=df_plot, x='Group', y='slope', order=order, palette=palette, showfliers=False, width=0.8, notch=True)
plt.axhline(0, color='black', linestyle='--', alpha=0.5)

# 使用 Annotator 标注 p 值
pairs = [("Ins_CGI", "Del_CGI"), ("Ins_CGI", "Ins_Others"), ("Del_CGI", "Del_Others")]
annotator = Annotator(ax, pairs, data=df_plot, x='Group', y='slope', order=order)
annotator.configure(test='Mann-Whitney', text_format='full', loc='outside')
annotator.apply_and_annotate()

plt.title('Downstream 2kb: Effect of SV Composition on Slope', fontsize=14)
plt.ylabel('Effect Size (Slope)')
plt.savefig('Downstream_2kb_Four_Groups.pdf')
plt.show()


'''