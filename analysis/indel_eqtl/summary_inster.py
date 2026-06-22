import pandas as pd
import pyarrow.parquet as pq
import glob
import os

# --- 1. 配置路径 ---
bed_path = '/Users/berial/ljy/lailab2024/features/u1_to_cgi_boundaries_both_strands_filtered_sorted.bed6'
data_dir = './' 
output_file = 'GTEx_v11_Insertions_Summary.csv'

# --- 2. 加载 BED 区域 ---
print("正在加载 BED 区域...")
bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
regions = pd.read_csv(bed_path, sep='\t', names=bed_cols)
regions['chrom'] = regions['chrom'].astype(str)

# --- 3. 定义处理函数 ---
def process_tissue(file_path, regions_df):
    tissue_name = os.path.basename(file_path).split('.')[0]
    print(f"正在处理组织: {tissue_name}...")
    
    # 根据你报错的信息，修正列名：将 gene_id 改为 phenotype_id
    columns = ['variant_id', 'phenotype_id', 'slope', 'pval_nominal', 'af']
    
    try:
        table = pq.read_table(file_path, columns=columns)
        df = table.to_pandas()
        
        # 重命名以便于后续处理
        df = df.rename(columns={'phenotype_id': 'gene_id'})
        
        # A. 提取插入片段 (Ref 长度 < Alt 长度)
        # GTEx variant_id: chr1_155205634_G_A_b38
        var_info = df['variant_id'].str.split('_', expand=True)
        
        # 确保有足够的列（处理可能的异常 ID）
        if var_info.shape[1] < 4:
            return pd.DataFrame()

        # 核心逻辑：Alt(3) 比 Ref(2) 长即为插入
        is_insertion = var_info[3].str.len() > var_info[2].str.len() + 10
        df = df[is_insertion].copy()
        
        # B. 转换坐标用于比对
        df['chrom'] = var_info[0]
        df['pos'] = var_info[1].astype(int)
        
        # C. 空间比对：判断插入位置是否在 BED 区域内
        # 1. 先通过染色体合并缩小范围
        matched = pd.merge(df, regions_df, on='chrom')
        # 2. 物理坐标过滤
        final_hits = matched[(matched['pos'] >= matched['start']) & (matched['pos'] <= matched['end'])].copy()
        
        final_hits['tissue'] = tissue_name
        
        return final_hits[['tissue', 'variant_id', 'gene_id', 'slope', 'pval_nominal', 'af', 'name', 'start', 'end']]
    
    except Exception as e:
        print(f"处理 {tissue_name} 失败: {e}")
        return pd.DataFrame()

# --- 4. 循环处理 ---
all_parquet_files = glob.glob(os.path.join(data_dir, "*.signif_pairs.parquet"))
results_list = []

for f in all_parquet_files:
    res = process_tissue(f, regions)
    if not res.empty:
        results_list.append(res)

# --- 5. 汇总保存 ---
if results_list:
    final_summary = pd.concat(results_list, ignore_index=True)
    # 按照显著性排序
    final_summary = final_summary.sort_values('pval_nominal')
    final_summary.to_csv(output_file, index=False)
    print(f"\n成功！共在 BED 区域内发现 {len(final_summary)} 个插入相关 eQTL 对。")
    print(final_summary.head())
else:
    print("\n未发现符合条件的插入片段。请检查：1. BED 文件染色体格式是否带 'chr'；2. 坐标版本是否同为 hg38。")
    
    
# %%

###################################### Plot FOR correlation between instert CpG% and expression change   

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# 1. 读取结果文件
df = pd.read_csv('GTEx_v11_Insertions_Summary.csv')

def classify_cgi(variant_id, cpg_pct_threshold=30.0):
    """提取插入序列并分类为 CGI 或 Non-CGI"""
    try:
        parts = variant_id.split('_')
        alt_seq = parts[3].upper()
        
        if len(alt_seq) < 2:
            return "Non-CGI"
        
        # 计算 CpG 百分比
        cg_count = alt_seq.count('CG')
        cpg_pct = (cg_count * 2) / len(alt_seq) * 100
        
        # 分类逻辑：根据 CpG 占比是否超过阈值
        # 你也可以改为逻辑：如果 cg_count >= 1 则为 CGI，否则为 Non-CGI
        if cpg_pct >= cpg_pct_threshold:
            return "CGI"
        else:
            return "Non-CGI"
    except:
        return None

# 2. 数据处理与分类
df['CGI_Status'] = df['variant_id'].apply(classify_cgi)
df = df.dropna(subset=['CGI_Status'])

# 3. 统计计算 (Mann-Whitney U Test)
# 提取两组数据
cgi_group = df[df['CGI_Status'] == 'CGI']['slope']
non_cgi_group = df[df['CGI_Status'] == 'Non-CGI']['slope']

# 计算 p-value
u_stat, p_val = stats.mannwhitneyu(cgi_group, non_cgi_group, alternative='two-sided')

print(f"CGI 组样本量: {len(cgi_group)}")
print(f"Non-CGI 组样本量: {len(non_cgi_group)}")
print(f"Mann-Whitney U 检验 p-value: {p_val:.2e}")

# 4. 绘图
plt.figure(figsize=(4, 5))
sns.set_style("whitegrid")

# 使用小提琴图结合箱线图
ax = sns.boxplot(data=df, x='CGI_Status', y='slope', palette="muted", showfliers=False)
# 或者使用简单的箱线图
# ax = sns.boxplot(data=df, x='CGI_Status', y='slope', palette="Set2")

# 在图上标注 p-value
plt.text(0.5, 1, f'p-value = {p_val:.2e}', 
         ha='center', va='bottom', fontsize=12, fontweight='bold', color='red')

plt.title('Effect of Insertion Type (CGI vs Non-CGI) on Slope', fontsize=14)
plt.xlabel('Insertion Category', fontsize=12)
plt.ylabel('Effect Size (Slope)', fontsize=12)

plt.axhline(0, color='black', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('CGI_vs_NonCGI_Slope.pdf')
plt.show()


# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# 1. 读取数据
df = pd.read_csv('GTEx_v11_Insertions_Summary.csv')

def classify_insertion(variant_id):
    """三分类逻辑：CGI, AT-rich, Others"""
    try:
        parts = variant_id.split('_')
        alt_seq = parts[3].upper()
        seq_len = len(alt_seq)
        if seq_len < 2: return "Others"
        
        # 计算 CpG%
        cg_count = alt_seq.count('CG')
        cpg_pct = (cg_count * 2) / seq_len * 100
        
        # 计算 AT%
        at_count = alt_seq.count('A') + alt_seq.count('T')
        at_pct = (at_count / seq_len) * 100
        
        # 分类阈值可根据需求微调
        if cpg_pct >= 20.0:
            return "CGI"
        elif at_pct >= 65.0:
            return "AT-rich"
        else:
            return "Others"
    except:
        return None

# 应用分类
df['Group'] = df['variant_id'].apply(classify_insertion)
df = df.dropna(subset=['Group', 'tissue'])

# 2. 筛选组织：只保留样本量足够的组织（例如至少有15个位点）
# 否则统计检验和绘图没有意义
min_sites = 15
tissue_counts = df['tissue'].value_counts()
valid_tissues = tissue_counts[tissue_counts >= min_sites].index.tolist()
df_filtered = df[df['tissue'].isin(valid_tissues)].copy()

# 为了美观，我们选取前 12 个组织进行展示（或根据需要调整）
display_tissues = valid_tissues[:12]
df_plot = df_filtered[df_filtered['tissue'].isin(display_tissues)]

# 3. 统计计算函数（用于在子图中显示各组差异的 p-value）
def get_p_value(data):
    groups = [data[data['Group'] == g]['slope'] for g in ['CGI', 'AT-rich', 'Others']]
    # 至少两组有数据才能做 Kruskal-Wallis 检验
    if sum(len(g) > 0 for g in groups) > 1:
        _, p = stats.kruskal(*[g for g in groups if len(g) > 0])
        return p
    return None

# 4. 绘图：使用 FacetGrid 按组织分面
plt.figure(figsize=(16, 12))
sns.set_style("white")

# 调色盘
palette = {"CGI": "#e74c3c", "AT-rich": "#3498db", "Others": "#95a5a6"}

g = sns.FacetGrid(df_plot, col="tissue", col_wrap=4, height=4, aspect=1, 
                  sharey=False) # sharey=False 因为不同组织的 slope 范围可能不同

# 映射箱线图
g.map_dataframe(sns.boxplot, x='Group', y='slope', order=['CGI', 'AT-rich', 'Others'],
                palette=palette, showfliers=False, width=0.6)

# 映射散点（增加透明度防止重叠）
g.map_dataframe(sns.stripplot, x='Group', y='slope', order=['CGI', 'AT-rich', 'Others'],
                palette=palette, alpha=0.4, size=3, jitter=True)

# 5. 装饰子图（添加 p-value 和 0 线）
for ax in g.axes.flat:
    tissue_name = ax.get_title().split('=')[-1].strip()
    ax.set_title(tissue_name, fontweight='bold', fontsize=12)
    ax.axhline(0, color='black', linestyle='--', alpha=0.3)
    
    # 计算并显示该组织的全局 p-value
    p = get_p_value(df_plot[df_plot['tissue'] == tissue_name])
    if p is not None:
        ax.text(0.5, 0.9, f'P={p:.2e}', transform=ax.transAxes, 
                ha='center', color='red', fontsize=10)

g.set_axis_labels("Category", "Slope")
plt.subplots_adjust(top=0.92, hspace=0.4)
g.fig.suptitle('Tissue-specific comparison of Insertion Composition on Expression', fontsize=16)

plt.savefig('Tissue_Comparison_CGI_AT_Others.png', dpi=300)
plt.show()