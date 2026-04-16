import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

# -------------------------- 1. 数据读取与清洗 --------------------------
file_path = 'D:/biocode/4further_research/IL6_INS_summary.csv'
df = pd.read_csv(file_path)

df = df[df['Species'] != 'Species']

df = df.dropna(subset=['CDS_Length','GC_Content']).reset_index(drop=True)# 这里假设CDS_Length和GC_Content列是关键指标，剔除这些列有空值的行

df['GC_Num'] = df['GC_Content'].astype(str).str.strip('%').astype(float)# 处理GC含量：去掉%转为数值，方便绘图和统计分析

# -------------------------- 2. 数据对齐 --------------------------
# 将数据拆分为两组，确保物种顺序是一致的
ins_group = df[df['Gene'] == 'INS'].sort_values('Species')
il6_group = df[df['Gene'] == 'IL6'].sort_values('Species')

# 提取数值
species = ins_group['Species'].tolist()
ins_vals = ins_group['GC_Num'].values
il6_vals = il6_group['GC_Num'].values

# -------------------------- 3. 统计计算 --------------------------
# 计算均值和标准差 (SD)
ins_mean, ins_std = np.mean(ins_vals), np.std(ins_vals, ddof=1)
il6_mean, il6_std = np.mean(il6_vals), np.std(il6_vals, ddof=1)

t_stat, p_val = stats.ttest_ind(ins_vals, il6_vals, equal_var=False)# 进行独立样本 T 检验，比较两组数据的均值差异是否显著
# -------------------------- 4. 可视化绘图 --------------------------
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS'] 
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(figsize=(12, 7))

# 绘制折线图
x = range(len(species))
ax.plot(x, ins_vals, marker='o', label=f'INS (代谢) Mean:{ins_mean:.1f} std:{ins_std:.2f}', 
        color='#ff7f0e', linewidth=2)
ax.plot(x, il6_vals, marker='s', label=f'IL-6 (免疫) Mean:{il6_mean:.1f} std:{il6_std:.2f}', 
        color='#1f77b4', linewidth=2)

# 装饰图表
ax.set_xticks(x)
ax.set_xticklabels(species, rotation=45, ha='right')
ax.set_title('免疫基因(IL-6)与代谢基因(INS) GC 含量在不同物种间的波动对比', fontsize=14)
ax.set_ylabel('CDS(coding sequence)GC 含量 (%)', fontsize=12)
ax.grid(alpha=0.3, linestyle='--')
ax.legend()

# 将 T 检验结果和 P 值标注在图中
stats_box = (f"统计检验结果:\n"
             f"T-statistic: {t_stat:.3f}\n"
             f"P-value: {p_val:.6f}\n"
             f"结论: {'显著差异' if p_val < 0.05 else '无显著差异'}")

# 在图左下角添加文本框
plt.text(0.05, 0.05, stats_box, transform=ax.transAxes, fontsize=10,
         bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

plt.tight_layout()
plt.show()

