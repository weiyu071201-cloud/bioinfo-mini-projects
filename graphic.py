import pandas as pd
import matplotlib.pyplot as plt

# -------------------------- 1. 读取并清洗数据（关键：剔除标题行、处理空值） --------------------------
df = pd.read_csv('D:/biocode/3cytokine/cytokine_summary.csv')
# 剔除空值行和标题行，重置索引
df = df.dropna(subset=['CDS_Length', 'GC_Content']).reset_index(drop=True)# CDS_Length和GC_Content列是关键指标，剔除这些列有空值的行
# 处理GC含量：去掉%转为数值
df['GC_Num'] = df['GC_Content'].str.strip('%').astype(float)

# -------------------------- 2. 分组（按你的：前4条=IL-6，后4条=INS） --------------------------
group1 = df.iloc[:4]  # IL-6序列组
group2 = df.iloc[4:]  # INS序列组
# 生成x轴标签
x_labels = ['Human', 'Mouse','Rat','cow']

# -------------------------- 3. 计算GC标准差 --------------------------
group1_std = group1['GC_Num'].std()
group2_std = group2['GC_Num'].std()

# -------------------------- 4. 绘图设置（中文字体+画布） --------------------------
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']  # 解决中文乱码
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
fig, axes = plt.subplots(1, 2, figsize=(14, 5))  # 1行2列子图

# -------------------------- 5. 子图1：CDS长度对比折线图 --------------------------
axes[0].plot(x_labels, group1['CDS_Length'], marker='o', color='#1f77b4', 
             linewidth=2, markersize=8, label='IL-6序列')
axes[0].plot(x_labels, group2['CDS_Length'], marker='s', color="#1f76b4", 
             linewidth=1, markersize=6, label='INS序列')
axes[0].set_title('IL-6 与 INS序列 CDS长度 对比', fontsize=14, pad=10)
axes[0].set_xlabel('基因编号', fontsize=12)
axes[0].set_ylabel('CDS长度 (bp)', fontsize=12)
axes[0].legend(loc='best')
axes[0].grid(alpha=0.3, linestyle='--')  # 加网格，更易看

# -------------------------- 6. 子图2：GC含量对比折线图（核心对比） --------------------------
axes[1].plot(x_labels, group1['GC_Num'], marker='o', color='#ff7f0e', 
             linewidth=2, markersize=8, label='IL-6序列')
axes[1].plot(x_labels, group2['GC_Num'], marker='s', color='#2ca02c', 
             linewidth=2, markersize=8, label='INS序列')
axes[1].set_title(f'IL-6 与 INS序列 GC含量 对比 (std:IL-6: {group1_std:.2f}, INS: {group2_std:.2f})', fontsize=14, pad=10)
axes[1].set_xlabel('基因编号', fontsize=12)
axes[1].set_ylabel('GC含量 (%)', fontsize=12)
axes[1].legend(loc='best')
axes[1].grid(alpha=0.3, linestyle='--')

# -------------------------- 7. 调整布局+显示/保存 --------------------------
plt.tight_layout()  # 自动调整间距，避免标签截断
# 可选：保存高清图片到指定路径，dpi=300适合报告/论文
plt.savefig('D:/biocode/3cytokine/IL-6_INS_CDS_mRNA_compare.png', dpi=300, bbox_inches='tight')
plt.show()