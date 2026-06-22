import pandas as pd
import pybedtools

# 路径设置
u1_bed = "/Users/berial/ljy/lailab2025/u1site/sm.bed"
cpg_bed = "/Users/berial/ljy/lailab2024/features/CpG/CpGtogene.bed"
output_tsv = "cpg_u1_distance.tsv"

# 读取 U1 和 CpG BED 文件
u1_binding_sites = pybedtools.BedTool(u1_bed)
cpg_islands = pybedtools.BedTool(cpg_bed)

# 计算 U1 与 CpG 的交集（保留链特异性）
intersections = u1_binding_sites.intersect(cpg_islands, wa=True, wb=True, s=True)

# 提取交集信息
data = []
for feature in intersections:
    u1_start = int(feature[1])
    u1_end = int(feature[2])
    u1_strand = feature[5]
    cpg_start = int(feature[7])
    cpg_end = int(feature[8])
    cpg_name = feature[9]

    u1_center = (u1_start + u1_end) / 2
    if u1_strand == "+":
        distance = cpg_end - u1_center
    else:
        distance = u1_center - cpg_start

    data.append((cpg_name, distance))

# 转为 DataFrame，取每个 CpG name 最大距离
df = pd.DataFrame(data, columns=["name", "distance_tocpg"])
df = df.groupby("name")["distance_tocpg"].max().reset_index()

# 加入 CpG 岛长度
cpg_df = pd.read_csv(cpg_bed, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand"])
cpg_df["cpg_length"] = cpg_df["end"] - cpg_df["start"]
df = df.merge(cpg_df[["name", "cpg_length"]], on="name", how="left")

# 保存为 TSV 文件
df.to_csv(output_tsv, sep="\t", index=False)
print(f"✅ Done. Output saved to {output_tsv}")
