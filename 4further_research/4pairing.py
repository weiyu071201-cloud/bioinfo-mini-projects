import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# -------------------------- 1. 读取与清洗数据 --------------------------
file_path = 'D:/biocode/4further_research/IL6_INS_summary.csv'
df = pd.read_csv(file_path)

# 清理夹杂的表头或空行
df = df.dropna(subset=['Species', 'GC_Content', 'Gene'])
df = df[df['Gene'].isin(['INS', 'IL6'])]
df['GC_Num'] = df['GC_Content'].astype(str).str.strip('%').astype(float)

# 重组为按物种对齐的宽表 (索引为物种，列为 IL6 和 INS)
pivot_df = df.pivot_table(index='Species', columns='Gene', values='GC_Num').dropna()

# -------------------------- 2. 物种分类 --------------------------
# 定义哺乳动物列表 (基于你之前数据中的物种)
mammals_list = [
    'Bos taurus', 'Canis lupus familiaris', 'Homo sapiens', 
    'Macaca mulatta', 'Mus musculus', 'Rattus norvegicus', 'Sus scrofa'
]

# 添加分组标签
pivot_df['Class'] = pivot_df.index.map(lambda x: 'Mammal' if x in mammals_list else 'Non-Mammal')

# 分离两组数据
mam_df = pivot_df[pivot_df['Class'] == 'Mammal']
non_mam_df = pivot_df[pivot_df['Class'] == 'Non-Mammal']

# -------------------------- 3. 统计学检验 (配对 T 检验) --------------------------
# 计算哺乳动物的统计学差异
t_mam, p_mam = stats.ttest_rel(mam_df['INS'], mam_df['IL6'])
mam_diff = mam_df['INS'].mean() - mam_df['IL6'].mean()

# 计算非哺乳动物的统计学差异
t_non, p_non = stats.ttest_rel(non_mam_df['INS'], non_mam_df['IL6'])
non_diff = non_mam_df['INS'].mean() - non_mam_df['IL6'].mean()

# -------------------------- 4. 可视化绘图 (配对连线图) --------------------------
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(8, 6))

# 遍历数据绘制配对连线
for species, row in pivot_df.iterrows():
    color = '#d62728' if row['Class'] == 'Mammal' else '#1f77b4' # 哺乳动物红色，非哺乳蓝色
    alpha = 0.8
    # 绘制线段与散点
    ax.plot(['IL-6 (免疫)', 'INS (代谢)'], [row['IL6'], row['INS']], 
            marker='o', color=color, alpha=alpha, linewidth=2, markersize=8)
    
    # 在 INS 右侧添加物种名称标签，避免图例拥挤
    ax.text(1.05, row['INS'], species, color=color, va='center', fontsize=9)

# 设置图形属性
ax.set_title('哺乳类 vs 非哺乳类: 免疫与代谢基因 GC 含量配对比较', fontsize=14, pad=15)
ax.set_ylabel('GC 含量 (%)', fontsize=12)
ax.set_xlim(-0.2, 1.5) # 留出右侧空间显示文字
ax.grid(axis='y', alpha=0.3, linestyle='--')

# -------------------------- 5. 组装检验结果框 --------------------------
stats_text = (
    f"【哺乳动物 (n={len(mam_df)}) - 红色】\n"
    f"INS 比 IL6 平均高: {mam_diff:.2f}%\n"
    f"P-value: {p_mam:.6f}\n\n"
    f"【非哺乳动物 (n={len(non_mam_df)}) - 蓝色】\n"
    f"INS 比 IL6 平均高: {non_diff:.2f}%\n"
    f"P-value: {p_non:.6f}"
)

ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))

plt.tight_layout()
plt.savefig('GC_paired_comparison.png', dpi=300)
plt.show()