import pandas as pd
import matplotlib.pyplot as plt

# -------------------------- 1. 读取并清洗数据（关键：剔除标题行、处理空值） --------------------------
df = pd.read_csv('D:/biocode/2project_insulin_analysis/insulin_summary.csv')
# 剔除空值行和标题行（CDS、record.seq这两行），重置索引
df = df.dropna(subset=['CDS_Length', 'GC_Content']).reset_index(drop=True)
# 处理GC含量：去掉%转为数值，方便绘图
df['GC_Num'] = df['GC_Content'].str.strip('%').astype(float)

# -------------------------- 2. 精准分组（按你的CSV实际结构：前8条=CDS序列，后8条=整个mRNA序列） --------------------------
group1 = df.iloc[:8]  # CDS序列组（标题为CDS的那组）
group2 = df.iloc[8:]  # 整个mRNA序列组（标题为record.seq的那组）
# 生成x轴标签（基因1-8，简化显示，避免文件名过长）
x_labels = [f'Gene {i+1}' for i in range(8)]

# -------------------------- 3. 绘图设置（中文字体+画布） --------------------------
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']  # 解决中文乱码
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
fig, axes = plt.subplots(1, 2, figsize=(14, 5))  # 1行2列子图，适配长度+GC两个指标

# -------------------------- 4. 子图1：CDS长度对比折线图 --------------------------
# 两组CDS长度一致，可只画一组，或叠加显示（这里叠加，颜色一致）
axes[0].plot(x_labels, group1['CDS_Length'], marker='o', color='#1f77b4', 
             linewidth=2, markersize=8, label='CDS序列')
axes[0].plot(x_labels, group2['CDS_Length'], marker='s', color='#1f77b4', 
             linewidth=1, markersize=6, label='mRNA全序列(同长度)')
axes[0].set_title('胰岛素基因 CDS长度 对比', fontsize=14, pad=10)
axes[0].set_xlabel('基因编号', fontsize=12)
axes[0].set_ylabel('CDS长度 (bp)', fontsize=12)
axes[0].legend(loc='best')
axes[0].grid(alpha=0.3, linestyle='--')  # 加网格，更易看

# -------------------------- 5. 子图2：GC含量对比折线图（核心对比） --------------------------
axes[1].plot(x_labels, group1['GC_Num'], marker='o', color='#ff7f0e', 
             linewidth=2, markersize=8, label='CDS序列')
axes[1].plot(x_labels, group2['GC_Num'], marker='s', color='#2ca02c', 
             linewidth=2, markersize=8, label='整个mRNA序列')
axes[1].set_title('胰岛素基因 GC含量 对比', fontsize=14, pad=10)
axes[1].set_xlabel('基因编号', fontsize=12)
axes[1].set_ylabel('GC含量 (%)', fontsize=12)
axes[1].legend(loc='best')
axes[1].grid(alpha=0.3, linestyle='--')

# -------------------------- 6. 调整布局+显示/保存 --------------------------
plt.tight_layout()  # 自动调整间距，避免标签截断
# 可选：保存高清图片到指定路径，dpi=300适合报告/论文
plt.savefig('D:/biocode/2project_insulin_analysis/insulin_CDS_mRNA_compare.png', dpi=300, bbox_inches='tight')
plt.show()