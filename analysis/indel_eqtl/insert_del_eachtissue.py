import pandas as pd
import pyarrow.parquet as pq
import glob
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from concurrent.futures import ThreadPoolExecutor, as_completed
from statannotations.Annotator import Annotator

# ================================================================
# --- 1. 配置
# ================================================================
bed_path    = '/Users/berial/ljy/lailab2024/features/u1_to_cgi_boundaries_both_strands_filtered_sorted.bed6'
data_dir    = './'
output_file = 'GTEx_v11_Ins_Del_Summary.csv'
N_WORKERS   = 8      # ★ 并行进程数，建议 CPU核心数-2，Apple M系列可设8~10
MIN_N       = 10     # 气泡图每组最小样本量
COLOR_INS   = '#2980b9'   # 插入=蓝
COLOR_DEL   = '#e74c3c'   # 删除=红

t_start = time.time()

# ================================================================
# --- 2. 读取BED，构建区间索引
# ================================================================
bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
regions  = pd.read_csv(bed_path, sep='\t', names=bed_cols)
regions['chrom'] = regions['chrom'].astype(str)
print(f"✅ BED区域加载完成，共 {len(regions):,} 个区域")

# 按染色体分组，按start排序，供二分查找
bed_by_chrom = {}
for chrom, grp in regions.groupby('chrom'):
    arr = grp[['start', 'end']].values
    bed_by_chrom[chrom] = arr[arr[:, 0].argsort()]


# ================================================================
# --- 3. 向量化区间查找（无需pyranges，兼容Apple Silicon）
# ================================================================
def variants_in_bed(chrom_arr, pos_arr, bed_by_chrom):
    """
    BED 0-based [start,end) ↔ GTEx 1-based pos
    条件: start < pos <= end
    """
    result = np.zeros(len(pos_arr), dtype=bool)
    for chrom in np.unique(chrom_arr):
        if chrom not in bed_by_chrom:
            continue
        intervals = bed_by_chrom[chrom]
        starts    = intervals[:, 0]
        ends      = intervals[:, 1]
        mask      = (chrom_arr == chrom)
        pos0      = pos_arr[mask] - 1
        idx       = np.searchsorted(starts, pos0, side='right') - 1
        valid     = (idx >= 0) & (pos0 < ends[np.clip(idx, 0, len(ends)-1)])
        valid[idx < 0] = False
        result[mask] = valid
    return result


# ================================================================
# --- 4. 单组织处理函数（供多进程调用）
# ================================================================
def process_tissue(args):
    file_path, bed_by_chrom = args
    tissue_name = os.path.basename(file_path).split('.')[0]

    columns = ['variant_id', 'phenotype_id', 'slope', 'pval_nominal', 'af']
    try:
        df = pq.read_table(file_path, columns=columns).to_pandas()
        df = df.rename(columns={'phenotype_id': 'gene_id'})

        # 向量化解析 variant_id
        split_vals = df['variant_id'].str.split('_', n=4, expand=True)
        if split_vals.shape[1] < 4:
            return pd.DataFrame(), tissue_name, 0, 0, 0

        df['chrom'] = split_vals[0].astype(str)
        df['pos']   = pd.to_numeric(split_vals[1], errors='coerce')
        df['ref']   = split_vals[2]
        df['alt']   = split_vals[3]
        df = df.dropna(subset=['pos'])
        df['pos'] = df['pos'].astype(int)

        # 向量化筛选 >10bp indel
        ref_len  = df['ref'].str.len().values
        alt_len  = df['alt'].str.len().values
        is_ins   = alt_len > ref_len + 10
        is_del   = ref_len > alt_len + 10
        df_indel = df[is_ins | is_del].copy()
        df_indel['sv_type'] = np.where(
            df_indel['alt'].str.len().values > df_indel['ref'].str.len().values + 10,
            'Insertion', 'Deletion'
        )

        if df_indel.empty:
            return pd.DataFrame(), tissue_name, 0, 0, 0

        # 区间查找
        in_bed      = variants_in_bed(df_indel['chrom'].values,
                                      df_indel['pos'].values, bed_by_chrom)
        final_hits  = df_indel[in_bed].copy()
        if final_hits.empty:
            return pd.DataFrame(), tissue_name, 0, 0, 0

        final_hits['tissue'] = tissue_name
        n_ins = (final_hits['sv_type'] == 'Insertion').sum()
        n_del = (final_hits['sv_type'] == 'Deletion').sum()
        return final_hits, tissue_name, len(final_hits), n_ins, n_del

    except Exception as e:
        return pd.DataFrame(), tissue_name, -1, 0, 0


# ================================================================
# --- 5. 多进程并行处理所有组织 ★ 核心提速
# ================================================================
# ================================================================
# --- 5. 串行处理所有组织（兼容所有环境）
#        向量化已大幅提速，串行也足够快
# ================================================================
all_files = sorted(glob.glob(os.path.join(data_dir, "*.signif_pairs.parquet")))
print(f"\n共找到 {len(all_files)} 个parquet文件\n")

args_list   = [(f, bed_by_chrom) for f in all_files]
results_raw = []

for i, arg in enumerate(args_list, 1):
    df_res, tname, n_hits, n_ins, n_del = process_tissue(arg)
    if n_hits > 0:
        results_raw.append(df_res)
        print(f"  [{i:>2}/{len(all_files)}] ✅ {tname:<50} {n_hits:>5,} hits "
              f"(ins={n_ins}, del={n_del})")
    elif n_hits == 0:
        print(f"  [{i:>2}/{len(all_files)}] ○  {tname:<50} 0 hits")
    else:
        print(f"  [{i:>2}/{len(all_files)}] ✗  {tname:<50} 处理失败")

print(f"\n并行处理完成，耗时 {time.time()-t_start:.1f}s")
df_main = pd.concat([r for r in results_raw if not r.empty], ignore_index=True)
print(f"✅ 共 {len(df_main):,} 行")


# ================================================================
# --- 6. 向量化分类（替代逐行apply，速度提升10-50x）★
# ================================================================
print("向量化分类中...")

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

df_main['Group']     = np.where(is_cgi, sv_label + '_CGI',
                        np.where(is_at_rich, sv_label + '_AT-rich',
                                             sv_label + '_Others'))
df_main['Group4']    = np.where(is_cgi, sv4_label + '_CGI', sv4_label + '_Others')
df_main['CGI_Group'] = np.where(is_cgi, sv_label + '_CGI', None)

df_main.to_csv(output_file, index=False)
print(f"✅ 主结果已保存: {output_file}")
print(df_main['Group'].value_counts())

df_final    = df_main[df_main['CGI_Group'].notna()].copy()
df_plot_all = df_main[df_main['Group4'].notna()].copy()
print(f"\nCGI子集: {len(df_final):,} 行  |  四分组: {len(df_plot_all):,} 行")


# ================================================================
# --- 7. 大分类映射
# ================================================================
TISSUE_CATEGORY = {
    'Brain_Amygdala':                        'Brain',
    'Brain_Anterior_cingulate_cortex_BA24':  'Brain',
    'Brain_Caudate_basal_ganglia':           'Brain',
    'Brain_Cerebellar_Hemisphere':           'Brain',
    'Brain_Cerebellum':                      'Brain',
    'Brain_Cortex':                          'Brain',
    'Brain_Frontal_Cortex_BA9':             'Brain',
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
    'Adipose_Visceral_Omentum':             'Adipose',
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

CATEGORY_ORDER = [
    'Brain', 'Heart', 'Artery', 'Adipose', 'Skin',
    'Colon', 'Esophagus', 'Cells', 'Reproductive',
    'Adrenal Gland', 'Bladder', 'Breast', 'Kidney',
    'Liver', 'Lung', 'Salivary Gland', 'Muscle',
    'Nerve', 'Pancreas', 'Pituitary', 'Small Intestine',
    'Spleen', 'Stomach', 'Thyroid', 'Whole Blood',
]

df_final['category']    = df_final['tissue'].map(TISSUE_CATEGORY)
df_plot_all['category'] = df_plot_all['tissue'].map(TISSUE_CATEGORY)


# ================================================================
# --- 8. 图A：全局箱线图
# ================================================================
sns.set_style("white")
target_order = ['Insertion_CGI', 'Deletion_CGI']
palette      = {"Insertion_CGI": COLOR_INS, "Deletion_CGI": COLOR_DEL}

plt.figure(figsize=(5, 6))
ax = sns.boxplot(data=df_final, x='CGI_Group', y='slope',
                 order=target_order, palette=palette,
                 showfliers=False, width=0.8, linewidth=2, notch=True)
g1 = df_final[df_final['CGI_Group'] == 'Insertion_CGI']['slope']
g2 = df_final[df_final['CGI_Group'] == 'Deletion_CGI']['slope']
_, p_val = stats.mannwhitneyu(g1, g2)
plt.text(0.5, df_final['slope'].quantile(0.95) * 0.8,
         f'p = {p_val:.2e}', ha='center', va='bottom', fontsize=12, fontweight='bold')
plt.axhline(0, color='black', linestyle='--', alpha=0.4)
plt.title('CGI Insertion vs Deletion (All Tissues)', fontsize=14)
plt.ylabel('Effect Size (Slope)')
plt.tight_layout()
plt.savefig('Global_CGI_Ins_vs_Del.pdf')
plt.show()


# ================================================================
# --- 9. 图B：大分类 FacetGrid（CGI Ins vs Del）
# ================================================================
df_facet_cat = df_final.dropna(subset=['category']).copy()
# 筛选两组各至少5个的大类
cat_ok = (df_facet_cat.groupby('category')['CGI_Group']
          .apply(lambda x: ((x=='Insertion_CGI').sum() >= 5) &
                            ((x=='Deletion_CGI').sum() >= 5)))
display_cats = [c for c in CATEGORY_ORDER if cat_ok.get(c, False)]
df_facet     = df_facet_cat[df_facet_cat['category'].isin(display_cats)]

g = sns.FacetGrid(df_facet, col='category', col_wrap=4,
                  height=4, aspect=0.8, sharey=False)
g.map_dataframe(sns.boxplot, x='CGI_Group', y='slope',
                order=target_order, palette=palette,
                showfliers=False, width=0.6)
for ax in g.axes.flat:
    cat_name = ax.get_title().split('=')[-1].strip()
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)
    t_data = df_facet[df_facet['category'] == cat_name]
    ins_s  = t_data[t_data['CGI_Group'] == 'Insertion_CGI']['slope']
    del_s  = t_data[t_data['CGI_Group'] == 'Deletion_CGI']['slope']
    if len(ins_s) > 3 and len(del_s) > 3:
        _, p = stats.mannwhitneyu(ins_s, del_s)
        ax.set_title(f"{cat_name}\n(p={p:.1e})", fontsize=10)
    else:
        ax.set_title(cat_name, fontsize=10)
    ax.set_xticklabels(['Ins\nCGI', 'Del\nCGI'], fontsize=8)
g.set_axis_labels("", "Slope")
plt.subplots_adjust(top=0.9, hspace=0.5)
g.fig.suptitle('Category-level CGI Insertion vs Deletion', fontsize=14)
plt.savefig('Category_CGI_Ins_vs_Del_Box.png', dpi=300, bbox_inches='tight')
plt.show()


# ================================================================
# --- 10. 图C：四分组全局对比
# ================================================================
order4   = ['Ins_CGI', 'Ins_Others', 'Del_CGI', 'Del_Others']
palette4 = {"Ins_CGI": COLOR_INS, "Ins_Others": "#85c1e9",
            "Del_CGI": COLOR_DEL, "Del_Others": "#f1948a"}

plt.figure(figsize=(9, 7))
ax = sns.boxplot(data=df_plot_all, x='Group4', y='slope',
                 order=order4, palette=palette4,
                 showfliers=False, width=0.6, linewidth=1.5, notch=True)
plt.axhline(0, color='black', linestyle='--', alpha=0.5, lw=1)
pairs = [("Ins_CGI", "Del_CGI"), ("Ins_CGI", "Ins_Others")]
annotator = Annotator(ax, pairs, data=df_plot_all,
                      x='Group4', y='slope', order=order4)
annotator.configure(test='Mann-Whitney', text_format='full',
                    loc='outside', show_test_name=False)
annotator.apply_and_annotate()
plt.ylabel('Effect Size (Slope)')
plt.xlabel('')
sns.despine()
plt.tight_layout()
plt.savefig('Four_Groups_Global_Annotated.pdf')
plt.show()


# ================================================================
# --- 11. 图D：四分组 × 大分类 FacetGrid
# ================================================================
df_plot_cat  = df_plot_all.dropna(subset=['category']).copy()
cat4_counts  = df_plot_cat.groupby('category').size()
display_cats4 = [c for c in CATEGORY_ORDER if cat4_counts.get(c, 0) >= 15][:12]
df_facet4     = df_plot_cat[df_plot_cat['category'].isin(display_cats4)]

g4 = sns.FacetGrid(df_facet4, col='category', col_wrap=4,
                   height=4, aspect=0.9, sharey=False)
g4.map_dataframe(sns.boxplot, x='Group4', y='slope',
                 order=order4, palette=palette4,
                 showfliers=False, width=0.7)
for ax in g4.axes.flat:
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)
    ax.set_xticklabels(order4, rotation=45, ha='right', fontsize=8)
    cat_name = ax.get_title().split('=')[-1].strip()
    t_data   = df_facet4[df_facet4['category'] == cat_name]
    i_cgi    = t_data[t_data['Group4'] == 'Ins_CGI']['slope']
    d_cgi    = t_data[t_data['Group4'] == 'Del_CGI']['slope']
    if len(i_cgi) > 3 and len(d_cgi) > 3:
        _, p = stats.mannwhitneyu(i_cgi, d_cgi)
        ax.set_title(f"{cat_name}\n(p={p:.1e})", fontsize=10)
g4.set_axis_labels("", "Slope")
plt.subplots_adjust(top=0.9, hspace=0.6)
plt.savefig('Four_Groups_Category_Facet.png', dpi=300, bbox_inches='tight')
plt.show()


# ================================================================
# --- 12. 气泡图：大分类 slope 中位数 & 均值
#         ● 插入=蓝  ■ 删除=红  气泡大小∝log10(n)  各组n≥MIN_N
# ================================================================
df_bubble = df_plot_all.dropna(subset=['category']).copy()

agg     = (df_bubble.groupby(['category', 'Group4'])
           .apply(lambda g: pd.Series({'n': len(g), 'median': g['slope'].median(),
                                        'mean': g['slope'].mean(), 'sem': g['slope'].sem()}))
           .reset_index())
agg_cgi = agg[agg['Group4'].isin(['Ins_CGI', 'Del_CGI'])].copy()

n_pivot    = agg_cgi.pivot_table(index='category', columns='Group4', values='n', fill_value=0)
valid_mask = (n_pivot.get('Ins_CGI', pd.Series(0, index=n_pivot.index)) >= MIN_N) & \
             (n_pivot.get('Del_CGI', pd.Series(0, index=n_pivot.index)) >= MIN_N)
valid_cats = set(valid_mask[valid_mask].index)

agg_filtered = agg_cgi[agg_cgi['category'].isin(valid_cats)].copy()
print(f"\n气泡图：{agg_cgi['category'].nunique()} 个大类 → 过滤后 {len(valid_cats)} 个（各组 n ≥ {MIN_N}）")
print(f"  保留: {sorted(valid_cats)}")

pval_dict = {}
for cat in valid_cats:
    sub = df_bubble[df_bubble['category'] == cat]
    _, p = stats.mannwhitneyu(sub[sub['Group4']=='Ins_CGI']['slope'],
                               sub[sub['Group4']=='Del_CGI']['slope'],
                               alternative='two-sided')
    pval_dict[cat] = p

agg_filtered['pval'] = agg_filtered['category'].map(pval_dict)
present_cats = [c for c in CATEGORY_ORDER if c in valid_cats]
y_pos        = {cat: i for i, cat in enumerate(reversed(present_cats))}

fig, axes = plt.subplots(1, 2,
                          figsize=(15, max(6, len(present_cats) * 0.55 + 2)),
                          gridspec_kw={'wspace': 0.06})

for ax_idx, metric in enumerate(['median', 'mean']):
    ax  = axes[ax_idx]
    sig_annotations = []

    for i in range(len(present_cats)):
        if i % 2 == 0:
            ax.axhspan(i-0.5, i+0.5, color='#f4f6f8', zorder=0)
    ax.axvline(0, color='#7f8c8d', lw=1, linestyle='--', alpha=0.7, zorder=1)

    for cat in present_cats:
        y       = y_pos[cat]
        ins_row = agg_filtered[(agg_filtered['category']==cat) & (agg_filtered['Group4']=='Ins_CGI')]
        del_row = agg_filtered[(agg_filtered['category']==cat) & (agg_filtered['Group4']=='Del_CGI')]

        if not ins_row.empty:
            xv = ins_row[metric].values[0]
            sz = np.clip(np.log10(ins_row['n'].values[0]+1)*130, 50, 700)
            ax.scatter(xv, y, s=sz, c=COLOR_INS, marker='o',
                       edgecolors='white', lw=0.8, alpha=0.88, zorder=4)
            ax.errorbar(xv, y, xerr=ins_row['sem'].values[0],
                        fmt='none', ecolor=COLOR_INS, elinewidth=1.5, capsize=4, alpha=0.75, zorder=3)
        if not del_row.empty:
            xv = del_row[metric].values[0]
            sz = np.clip(np.log10(del_row['n'].values[0]+1)*130, 50, 700)
            ax.scatter(xv, y, s=sz, c=COLOR_DEL, marker='s',
                       edgecolors='white', lw=0.8, alpha=0.88, zorder=4)
            ax.errorbar(xv, y, xerr=del_row['sem'].values[0],
                        fmt='none', ecolor=COLOR_DEL, elinewidth=1.5, capsize=4, alpha=0.75, zorder=3)
        if not ins_row.empty and not del_row.empty:
            ax.plot([ins_row[metric].values[0], del_row[metric].values[0]], [y, y],
                    color='#95a5a6', lw=1.2, alpha=0.5, zorder=2)

        p   = pval_dict.get(cat, np.nan)
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        sig_annotations.append((y, sig))

    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(list(reversed(present_cats)), fontsize=10.5)
    if ax_idx == 1:
        ax.set_yticklabels([])
    ax.set_xlabel(f'Slope {metric.capitalize()} (±SEM)', fontsize=12)
    ax.set_title(f'{"Median" if metric=="median" else "Mean"} Effect Size',
                 fontsize=13, fontweight='bold', pad=8)
    ax.set_ylim(-0.8, len(present_cats)-0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    xlim    = ax.get_xlim()
    x_range = xlim[1] - xlim[0]
    x_sig   = xlim[1] + x_range * 0.04
    for y_s, sig in sig_annotations:
        ax.text(x_sig, y_s, sig, va='center', ha='left', fontsize=9,
                color='#c0392b' if sig != 'ns' else '#95a5a6',
                fontweight='bold' if sig != 'ns' else 'normal')
    ax.set_xlim(xlim[0], x_sig + x_range * 0.12)

legend_elements = [
    plt.scatter([], [], s=120, c=COLOR_INS, marker='o',
                edgecolors='white', label=f'Insertion CGI  (n ≥ {MIN_N})'),
    plt.scatter([], [], s=120, c=COLOR_DEL, marker='s',
                edgecolors='white', label=f'Deletion CGI   (n ≥ {MIN_N})'),
    mpatches.Patch(color='white', label=''),
    plt.scatter([], [], s=np.log10(11)*130,  c='gray', marker='o', alpha=0.6, label='n = 10'),
    plt.scatter([], [], s=np.log10(51)*130,  c='gray', marker='o', alpha=0.6, label='n = 50'),
    plt.scatter([], [], s=np.log10(201)*130, c='gray', marker='o', alpha=0.6, label='n = 200'),
    mpatches.Patch(color='white', label='bubble ∝ log₁₀(n)'),
    mpatches.Patch(color='white', label=''),
    mpatches.Patch(color='#c0392b', label='*** p < 0.001', alpha=0.8),
    mpatches.Patch(color='#c0392b', label='**  p < 0.01',  alpha=0.6),
    mpatches.Patch(color='#c0392b', label='*   p < 0.05',  alpha=0.4),
    mpatches.Patch(color='#95a5a6', label='ns  p ≥ 0.05',  alpha=0.6),
]
fig.legend(handles=legend_elements, loc='center right', bbox_to_anchor=(1.20, 0.5),
           fontsize=9, frameon=True, framealpha=0.95, edgecolor='#bdc3c7',
           title='Legend', title_fontsize=10)
fig.suptitle(
    f'CGI Insertion (blue ●) vs Deletion (red ■): Slope by Tissue Category\n'
    f'GTEx v11 | >10bp Indels in BED regions | Min n = {MIN_N} per group',
    fontsize=12, fontweight='bold', y=1.02
)
plt.tight_layout()
plt.savefig('Bubble_Plot_Category_Slope.pdf', bbox_inches='tight', dpi=300)
plt.savefig('Bubble_Plot_Category_Slope.png', bbox_inches='tight', dpi=300)
plt.show()
print(f"✅ 气泡图已保存（{len(present_cats)} 个大类）")

agg_filtered.to_csv('Bubble_Category_Summary.csv', index=False)
print("✅ 汇总表: Bubble_Category_Summary.csv")
print(agg_filtered[['category','Group4','n','median','mean','sem','pval']]
      .sort_values(['category','Group4']).to_string(index=False))

print(f"\n🎉 全部完成，总耗时 {time.time()-t_start:.1f}s")
