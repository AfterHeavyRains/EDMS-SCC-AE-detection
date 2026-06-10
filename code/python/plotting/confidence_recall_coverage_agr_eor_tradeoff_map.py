import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.ticker import AutoMinorLocator
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.title_fontsize": 10,
    "axes.linewidth": 0.9,
    "lines.linewidth": 0.9,
    "grid.linewidth": 0.6,
    "grid.alpha": 0.5,
    "axes.edgecolor": "black",
    "xtick.direction": "in",
    "ytick.direction": "in",
    "figure.dpi": 300,
     "savefig.dpi": 650,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})
matplotlib.use('Agg')# ==== 数据 ====
thickness = np.array(["2mm", "4mm", "6mm", "8mm"])
AGR_Orig = [0.795, 0.812, 0.597, 0.725]
AGR_Layered = [0.96, 0.966, 0.91, 0.94]
EOR_Orig = [0.569, 0.655, 0.781, 0.623]
EOR_Layered = [0.071, 0.133, 0.070, 0.073]
Coverage_Orig = [0.562, 0.532, 0.411, 0.514]
Coverage_Layered = [0.929, 0.867, 0.930, 0.927]

x = np.arange(len(thickness))
width = 0.34

# ==== 绘图 ====
fig, ax1 = plt.subplots(figsize=(5, 3), dpi=320)

# ===== 分组绘制：针对每个厚度单独控制堆叠顺序 =====
for i in range(len(thickness)):
    # ---- Original ----
    if i == 2:  # 特例：6mm，让 AGR 在上
        ax1.bar(x[i] - width/2, EOR_Orig[i], width, color='#c38323', alpha=1, edgecolor='none', label=None if i>0 else 'EOR (Base)')
        ax1.bar(x[i] - width/2, AGR_Orig[i], width, bottom=0, color='#75b9e6', alpha=1, edgecolor='none', label=None if i>0 else 'AGR (Base)')
    else:       # 其他厚度正常（AGR 在下）
        ax1.bar(x[i] - width/2, AGR_Orig[i], width, color='#75b9e6', alpha=1, edgecolor='none', label=None if i>0 else 'AGR (Base)')
        ax1.bar(x[i] - width/2, EOR_Orig[i], width, bottom=0, color='#c38323', alpha=1, edgecolor='none', label=None if i>0 else 'EOR (Base)')

    # ---- Layered ----
    if i == 2:  # 特例：6mm，让 AGR 在上
        ax1.bar(x[i] + width/2, EOR_Layered[i], width, color='#8c7357', alpha=1, edgecolor='none', label=None if i>0 else 'EOR (EDMS)')
        ax1.bar(x[i] + width/2, AGR_Layered[i], width, bottom=EOR_Layered[i], color='#4889b3', alpha=1, edgecolor='none', label=None if i>0 else 'AGR (EDMS)')
    else:
        ax1.bar(x[i] + width/2, AGR_Layered[i], width, color='#4889b3', alpha=1, edgecolor='none', label=None if i>0 else 'AGR (EDMS)')
        ax1.bar(x[i] + width/2, EOR_Layered[i], width, bottom=0, color='#8c7357', alpha=1, edgecolor='none', label=None if i>0 else 'EOR (EDMS)')

# ==== 主轴设置 ====
ax1.set_ylabel('AGR / EOR', fontsize=10)
ax1.set_xlabel('Specimen Thickness', fontsize=10)
ax1.set_xticks(x)
ax1.set_xticklabels(thickness, fontsize=9)
ax1.set_ylim(0, 1.2)
ax1.grid(axis='y', linestyle=':', alpha=0.5)

# ==== Coverage 折线图 ====
ax2 = ax1.twinx()
ax2.plot(x, Coverage_Orig, 'o--', color='#6a3d9a', label='Cov (Base)', linewidth=1.0, markersize=4)
ax2.plot(x, Coverage_Layered, 'o-', color='#00c9ad', label='Cov (EDMS)', linewidth=1.0, markersize=4)
ax2.set_ylabel('Coverage', fontsize=8)
ax2.set_ylim(0, 1.3)

# ==== 自定义图例（透明匹配） ====
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    # --- Baseline 三个 ---
    Patch(facecolor='#75b9e6', alpha=1, label='AGR (Baseline)'),
    Patch(facecolor='#4889b3', alpha=1, label='AGR (Proposed)'),
    Patch(facecolor='#c38323', alpha=0.8, label='EOR (Baseline)'),
    Patch(facecolor='#8c7357', alpha=0.8, label='EOR (Proposed)'),
    Line2D([0], [0], color='#6a3d9a', marker='o', linestyle='--', label='Cov (Baseline)'),
    Line2D([0], [0], color='#00c9ad', marker='o', linestyle='-', label='Cov (Proposed)')
]
ax1.legend(handles=legend_elements, loc='upper left', ncol=3, fontsize=8.5, frameon=False, handlelength=1.6, columnspacing=1.2).set_zorder(10)

# 主坐标轴刻度朝内
ax1.tick_params(axis='both', direction='in', length=4, width=0.8, top=False, right=False)
# 右轴也要朝内
ax2.tick_params(axis='both', direction='in', length=4, width=0.8,  labelsize=9,top=False, right=True)

# ==== 标题 ====
# plt.title('Confidence–Recall and Coverage Comparison Across Thickness', fontsize=13, pad=10)
plt.tight_layout()
plt.show()


# === 映射颜色和尺寸 ===
colors = ['#5dade2', '#2874a6', '#1b4f72', '#154360']
sizes = [60, 70, 80, 90]  # 不同厚度的点大小

# === 绘图 ===
fig, ax = plt.subplots(figsize=(5.2, 3.8), dpi=300)

for i in range(len(thickness)):
    # 连接线
    ax.plot([EOR_Orig[i], EOR_Layered[i]],
            [AGR_Orig[i], AGR_Layered[i]],
            linestyle='--', color=colors[i], alpha=0.75, linewidth=1.0)
    # 原始框架：圆点
    ax.scatter(EOR_Orig[i], AGR_Orig[i],
               s=sizes[i], color=colors[i], edgecolor='k', linewidth=0.5,
               marker='o', label=None if i>0 else 'Baseline')
    # 分层框架：三角形
    ax.scatter(EOR_Layered[i], AGR_Layered[i],
               s=sizes[i], color=colors[i], edgecolor='k', linewidth=0.5,
               marker='^', label=None if i>0 else 'Proposed (EDMS)')

# === 坐标轴与网格 ===
ax.set_xlabel('Error Occurrence Rate (EOR)', fontsize=13)
ax.set_ylabel('Accuracy Gain Ratio (AGR)', fontsize=13)
ax.set_xlim(0, 1.0)
ax.set_ylim(0, 1.0)
ax.grid(True, linestyle=':', alpha=0.4)
ax.tick_params(direction='in', length=4, width=0.8, top=True, right=True, labelsize=12)
from matplotlib.lines import Line2D

# === 自定义整合图例 ===
legend_elements = [
    # 空心部分（方法）
    Line2D([0], [0], marker='o', color='none', markerfacecolor='none',
           markeredgecolor='k', markersize=6, label='Baseline'),
    Line2D([0], [0], marker='^', color='none', markerfacecolor='none',
           markeredgecolor='k', markersize=6, label='Proposed (EDMS)')
]

# === 厚度部分 ===
for i in range(len(thickness)):
    legend_elements.append(
        Patch(
            facecolor=colors[i],   # 纯色块
            edgecolor='none',      # 去掉边框
            label=thickness[i],
            alpha=1.0
        )
    )


# === 添加图例 ===
ax.legend(handles=legend_elements,
          loc='lower right',
          fontsize=11,
          frameon=False,
          title='Framework & Thickness',
          title_fontsize=11,
          handletextpad=0.8,  # 图例间距
          labelspacing=0.4,   # 行距
          borderaxespad=0.6)


# === 标题 ===
# ax.set_title('AGR–EOR Trade-off Map: Original vs Layered Framework', fontsize=13, pad=8)

plt.tight_layout()
plt.show()
